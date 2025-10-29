import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import Table

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
    query = f"""        SELECT source_id, ra, dec, phot_g_mean_mag
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {r})
        )=1
    """
    job = Gaia.launch_job(query)
    t = job.get_results()
    t.rename_columns(["ra","dec","phot_g_mean_mag"],["RA_ICRS","DE_ICRS","Gmag"])
    return Table(t)

def crossmatch(detections, wcs, catalog, match_radius_arcsec=1.5):
    if len(detections) == 0 or len(catalog) == 0:
        return Table()
    ra, dec = wcs.pixel_to_world(detections["xcentroid"], detections["ycentroid"])
    det_coords = SkyCoord(ra.deg*u.deg, dec.deg*u.deg, frame="icrs")
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
