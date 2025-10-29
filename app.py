#!/usr/bin/env python3



import argparse
import os
import time
import json
import shutil
import tempfile
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

def log(msg):
    print(f"[star-id] {msg}", flush=True)

def has_wcs(header):
    try:
        w = WCS(header)
        return w.has_celestial
    except Exception:
        return False

def ensure_dir(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

def estimate_background(data, box_size=64, filter_size=3, sigma_clip=3.0):
    sigma_clip = SigmaClip(sigma=sigma_clip)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, box_size=box_size, filter_size=filter_size,
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return bkg.background, bkg.background_rms

def detect_sources(data, fwhm=3.0, threshold_sigma=5.0, bkg=None, rms=None,
                   sharplo=0.2, sharphi=1.5):
    if bkg is None or rms is None:
        bkg, rms = estimate_background(data)
    data_sub = data - bkg
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * float(np.nanmedian(rms)),
                            sharplo=sharplo, sharphi=sharphi)
    tbl = daofind(data_sub)
    if tbl is None or len(tbl) == 0:
        from astropy.table import Table
        return Table(names=["xcentroid", "ycentroid", "flux"], dtype=[float, float, float])
    tbl.sort("flux")
    tbl.reverse()
    return tbl[["xcentroid", "ycentroid", "flux"]]

def image_center_and_radius(wcs, shape):
    ny, nx = shape
    corners = np.array([[0, 0], [0, ny-1], [nx-1, ny-1], [nx-1, 0]])
    world = wcs.pixel_to_world(corners[:, 0], corners[:, 1])
    center = wcs.pixel_to_world(nx/2.0, ny/2.0)
    center = SkyCoord(center.ra.deg*u.deg, center.dec.deg*u.deg, frame="icrs")
    sep = center.separation(world)
    radius = np.max(sep)
    return center, radius

def query_gaia(center, radius):
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    r = (radius * 1.05).to(u.deg).value
    query = f"""
        SELECT source_id, ra, dec, phot_g_mean_mag
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {r})
        )=1
    """
    job = Gaia.launch_job(query)
    t = job.get_results()
    t.rename_columns(["ra", "dec", "phot_g_mean_mag"], ["RA_ICRS", "DE_ICRS", "Gmag"])
    from astropy.table import Table
    return Table(t)

def crossmatch(detections, wcs, catalog, match_radius_arcsec=1.0):
    from astropy.table import Table
    if len(detections) == 0 or len(catalog) == 0:
        return Table()

    ra, dec = wcs.pixel_to_world(detections["xcentroid"], detections["ycentroid"])
    det_coords = SkyCoord(ra.deg*u.deg, dec.deg*u.deg, frame="icrs")

    if "RA_ICRS" not in catalog.colnames or "DE_ICRS" not in catalog.colnames:
        raise ValueError("Catalog must contain RA_ICRS, DE_ICRS in degrees.")

    cat_coords = SkyCoord(catalog["RA_ICRS"]*u.deg, catalog["DE_ICRS"]*u.deg, frame="icrs")

    idx, sep2d, _ = det_coords.match_to_catalog_sky(cat_coords)
    mask = sep2d <= (match_radius_arcsec * u.arcsec)

    out = Table()
    out["x"] = detections["xcentroid"]
    out["y"] = detections["ycentroid"]
    out["RA_deg"] = det_coords.ra.deg
    out["Dec_deg"] = det_coords.dec.deg
    out["sep_arcsec"] = sep2d.to(u.arcsec)
    out["matched"] = mask

    for col in catalog.colnames:
        is_str = catalog[col].dtype.kind in "OUS"
        out[col] = np.full(len(out), "" if is_str else np.nan, dtype=object if is_str else float)

    for i, ok in enumerate(mask):
        if ok:
            j = idx[i]
            for col in catalog.colnames:
                out[col][i] = catalog[col][j]

    return out

def annotate_plot(data, wcs, det_tbl, matches, out_png):
    img = data.copy().astype(float)
    p1 = np.nanpercentile(img, 1)
    p99 = np.nanpercentile(img, 99)
    img = (img - p1) / (p99 - p1 + 1e-9)
    img = np.clip(img, 0, 1)

    fig = plt.figure(figsize=(9, 7), dpi=150)
    ax = plt.subplot(projection=wcs)
    ax.imshow(img, origin="lower", cmap="gray")
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    ax.set_title("Detections and Catalog Matches")

    ax.scatter(det_tbl["xcentroid"], det_tbl["ycentroid"], transform=ax.get_transform('pixel'),
               s=20, marker='o', facecolors='none', edgecolors='C0', linewidths=1.2, label="Detections")

    if len(matches) > 0:
        matched_px = np.array(matches["matched"], dtype=bool)
        ax.scatter(det_tbl["xcentroid"][matched_px], det_tbl["ycentroid"][matched_px],
                   transform=ax.get_transform('pixel'),
                   s=30, marker='x', color='C1', linewidths=1.2, label="Matched")
        for i in np.where(matched_px)[0]:
            x = det_tbl["xcentroid"][i]
            y = det_tbl["ycentroid"][i]
            label = None
            if "source_id" in matches.colnames and str(matches["source_id"][i]) not in ("", "nan"):
                label = str(matches["source_id"][i])
            elif "Gmag" in matches.colnames:
                try:
                    g = float(matches["Gmag"][i])
                    label = f"G={g:.2f}"
                except Exception:
                    label = "match"
            else:
                label = "match"
            ax.text(x+6, y+6, label, transform=ax.get_transform('pixel'), fontsize=7, color='w',
                    bbox=dict(facecolor='k', alpha=0.35, pad=1.0, edgecolor='none'))

    ax.legend(loc="upper right", frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

def solve_with_local_astrometry(input_path: Path, out_dir: Path) -> Path:
    if shutil.which("solve-field") is None:
        raise RuntimeError("solve-field not found on PATH.")
    out_base = out_dir / (input_path.stem + "_wcs")
    cmd = [
        "solve-field",
        str(input_path),
        "--overwrite",
        "--no-plots",
        "--downsample", "2",
        "--dir", str(out_dir),
        "--new-fits", str(out_base.with_suffix(".fits")),
    ]
    log(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        log(res.stdout)
        log(res.stderr)
        raise RuntimeError("solve-field failed.")
    wcs_fits = out_base.with_suffix(".fits")
    if not wcs_fits.exists():
        raise RuntimeError("Expected WCS FITS not found after solve-field.")
    return wcs_fits

def solve_with_astrometry_web(input_path: Path, out_dir: Path, api_key: str) -> Path:
    import requests
    log("Logging in to astrometry.net API...")
    r = requests.post("https://nova.astrometry.net/api/login", data={"request-json": json.dumps({"apikey": api_key})})
    r.raise_for_status()
    sess = r.json()["session"]
    log("Uploading image...")
    files = {"file": open(input_path, "rb")}
    payload = {"request-json": json.dumps({"session": sess})}
    r = requests.post("https://nova.astrometry.net/api/upload", data=payload, files=files)
    r.raise_for_status()
    upload_resp = r.json()
    if "status" in upload_resp and upload_resp["status"] == "error":
        raise RuntimeError(f"Upload error: {upload_resp}")
    subid = upload_resp.get("subid", None)
    if subid is None:
        raise RuntimeError(f"No submission id: {upload_resp}")
    log(f"Submission id: {subid}")
    job_id = None
    for _ in range(120):
        time.sleep(5)
        r = requests.get(f"https://nova.astrometry.net/api/submissions/{subid}")
        r.raise_for_status()
        jobs = r.json().get("jobs", [])
        jobs = [j for j in jobs if j is not None]
        if jobs:
            job_id = jobs[0]
            break
        log("Waiting for job assignment...")
    if job_id is None:
        raise RuntimeError("Timed out waiting for job id.")
    while True:
        time.sleep(5)
        r = requests.get(f"https://nova.astrometry.net/api/jobs/{job_id}")
        r.raise_for_status()
        status = r.json().get("status", "")
        log(f"Job {job_id} status: {status}")
        if status in ("success", "failure"):
            break
    if status != "success":
        raise RuntimeError("Plate solve failed on server.")
    wcs_url = f"https://nova.astrometry.net/wcs_file/{job_id}"
    log(f"Downloading WCS FITS: {wcs_url}")
    r = requests.get(wcs_url)
    r.raise_for_status()
    wcs_fits = out_dir / (input_path.stem + "_wcs.fits")
    with open(wcs_fits, "wb") as f:
        f.write(r.content)
    return wcs_fits

def ensure_wcs(input_path: Path, out_dir: Path) -> Path:
    if input_path.suffix.lower() in (".fits", ".fit", ".fts"):
        with fits.open(input_path) as hdul:
            for h in hdul:
                if h.data is not None and h.data.ndim == 2 and has_wcs(h.header):
                    return input_path
        raise RuntimeError("FITS provided but no celestial WCS found.")
    api_key = os.environ.get("ASTROMETRY_API_KEY", "").strip()
    if shutil.which("solve-field") is not None:
        log("Using local astrometry.net 'solve-field'...")
        return solve_with_local_astrometry(input_path, out_dir)
    if api_key:
        log("Using astrometry.net web API...")
        return solve_with_astrometry_web(input_path, out_dir, api_key)
    raise RuntimeError("Need WCS. Provide FITS+WCS, install 'solve-field', or set ASTROMETRY_API_KEY.")

def main():
    ap = argparse.ArgumentParser(description="Detect and identify stars (FITS+WCS or JPEG/PNG).")
    ap.add_argument("image_path", help="Path to input image (FITS, JPG, PNG).")
    ap.add_argument("--out-prefix", default="stars", help="Output prefix for CSV/PNG")
    ap.add_argument("--fwhm", type=float, default=3.0, help="FWHM (pixels) for DAOStarFinder")
    ap.add_argument("--threshold", type=float, default=5.0, help="Detection threshold (sigmas)")
    ap.add_argument("--match-radius", type=float, default=1.5, help="Match radius (arcsec)")
    args = ap.parse_args()

    in_path = Path(args.image_path).expanduser().resolve()
    out_prefix = Path(args.out_prefix).resolve()
    ensure_dir(out_prefix)

    log(f"Input: {in_path}")
    wcs_fits = ensure_wcs(in_path, out_prefix.parent)
    log(f"WCS FITS: {wcs_fits}")

    with fits.open(wcs_fits) as hdul:
        hdu = None
        for h in hdul:
            if h.data is not None and h.data.ndim == 2 and has_wcs(h.header):
                hdu = h
                break
        if hdu is None:
            raise RuntimeError("No 2D WCS image found in FITS.")
        data = np.array(hdu.data, dtype=float)
        wcs = WCS(hdu.header)

    bkg, rms = estimate_background(data)
    det = detect_sources(data, fwhm=args.fwhm, threshold_sigma=args.threshold, bkg=bkg, rms=rms)
    log(f"Detected {len(det)} sources.")

    center, radius = image_center_and_radius(wcs, data.shape)
    log(f"Field center: {center.to_string('hmsdms')}, radius ≈ {radius.to(u.arcmin):.1f}")

    log("Querying Gaia DR3...")
    cat = query_gaia(center, radius)
    log(f"Gaia rows: {len(cat)}")

    matches = crossmatch(det, wcs, cat, match_radius_arcsec=args.match_radius)
    n_match = int(np.count_nonzero(matches['matched'])) if len(matches) else 0
    log(f"Matched {n_match} / {len(det)} within {args.match_radius} arcsec.")

    out_csv = f"{args.out_prefix}_matches.csv"
    matches.write(out_csv, overwrite=True)
    log(f"Wrote: {out_csv}")

    out_png = f"{args.out_prefix}_annotated.png"
    annotate_plot(data, wcs, det, matches, out_png)
    log(f"Wrote: {out_png}")

if __name__ == "__main__":
    main()
