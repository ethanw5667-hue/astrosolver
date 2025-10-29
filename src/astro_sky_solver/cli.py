import argparse
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .solver import ensure_wcs
from .detect import estimate_background, detect_sources
from .catalog import image_center_and_radius, query_gaia, crossmatch
from .annotate import annotate_plot
from .extinction import estimate_zero_point, detection_flux_threshold, limiting_magnitude_from_threshold, classify_bortle_like
from .utils import log, ensure_dir_for

def main():
    ap = argparse.ArgumentParser(description="Detect + identify stars and estimate limiting magnitude / sky quality.")
    ap.add_argument("image_path", help="Path to input image (FITS, JPG, PNG).")
    ap.add_argument("--out-prefix", default="stars", help="Output prefix for CSV/PNG")
    ap.add_argument("--fwhm", type=float, default=3.0, help="Stellar FWHM in pixels (for detection)")
    ap.add_argument("--threshold", type=float, default=5.0, help="Detection threshold in sigmas")
    ap.add_argument("--match-radius", type=float, default=1.5, help="Cross-match radius (arcsec)")
    args = ap.parse_args()

    in_path = Path(args.image_path).expanduser().resolve()
    out_prefix = Path(args.out_prefix).expanduser().resolve()
    ensure_dir_for(out_prefix)

    log(f"Input: {in_path}")
    wcs_fits = ensure_wcs(in_path, out_prefix.parent)
    log(f"WCS FITS: {wcs_fits}")

    with fits.open(wcs_fits) as hdul:
        hdu = None
        for h in hdul:
            if h.data is not None and getattr(h.data, 'ndim', 0) == 2:
                try:
                    w = WCS(h.header)
                    if w.has_celestial:
                        hdu = h
                        break
                except Exception:
                    pass
        if hdu is None:
            raise RuntimeError("No 2D WCS image found in FITS.")
        data = np.array(hdu.data, dtype=float)
        wcs = WCS(hdu.header)

    # Detect
    bkg, rms = estimate_background(data)
    det = detect_sources(data, fwhm=args.fwhm, threshold_sigma=args.threshold, bkg=bkg, rms=rms)
    log(f"Detected {len(det)} sources.")

    # Field & query
    center, radius = image_center_and_radius(wcs, data.shape)
    log(f"Field center: {center.to_string('hmsdms')}, radius ≈ {radius.to('arcmin'):.1f}")
    log("Querying Gaia DR3...")
    cat = query_gaia(center, radius)
    log(f"Gaia rows: {len(cat)}")

    # Cross-match
    matches = crossmatch(det, wcs, cat, match_radius_arcsec=args.match_radius)
    n_match = int(np.count_nonzero(matches['matched'])) if len(matches) else 0
    log(f"Matched {n_match} / {len(det)} within {args.match_radius} arcsec.")

    # Photometry + limits
    zp = estimate_zero_point(det, matches)
    glim = None
    if zp is not None and len(det) > 0:
        flux_thresh = detection_flux_threshold(rms=float(np.nanmedian(rms)), fwhm=args.fwhm, snr_thresh=args.threshold)
        glim = limiting_magnitude_from_threshold(zp, flux_thresh)
        log(f"Estimated camera limiting magnitude (Gaia G): {glim:.2f} mag")
    sky_class = classify_bortle_like(glim) if glim is not None else "unknown"
    log(f"Rough sky quality (Bortle-like): {sky_class}")

    # Save outputs
    out_csv = f"{args.out_prefix}_matches.csv"
    matches.write(out_csv, overwrite=True)
    log(f"Wrote: {out_csv}")

    out_png = f"{args.out_prefix}_annotated.png"
    annotate_plot(data, wcs, det, matches, out_png,
                  subtitle=(f"G_lim≈{glim:.2f} | Sky≈{sky_class}" if glim is not None else None))
    log(f"Wrote: {out_png}")

if __name__ == "__main__":
    main()
