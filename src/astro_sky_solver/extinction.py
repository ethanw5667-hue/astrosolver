import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def estimate_zero_point(det_tbl, matches):
    if len(det_tbl) == 0 or len(matches) == 0:
        return None
    vals = []
    for i in range(len(matches)):
        if not matches['matched'][i]: continue
        flux = det_tbl['flux'][i]
        try: gmag = float(matches['Gmag'][i])
        except Exception: continue
        if not np.isfinite(flux) or flux <= 0 or not np.isfinite(gmag): continue
        minst = -2.5 * np.log10(flux)
        vals.append(gmag - minst)
    if len(vals) < 3: return None
    return float(np.median(vals))

def detection_flux_threshold(rms: float, fwhm: float, snr_thresh: float) -> float:
    npix = np.pi * (max(fwhm, 1.0) / 2.0)**2
    return float(max(snr_thresh, 1.0) * np.sqrt(npix) * rms)

def limiting_magnitude_from_threshold(zp: float, flux_thresh: float) -> float:
    return float(zp - 2.5 * np.log10(max(flux_thresh, 1e-12)))

def airmass_from_center(center: SkyCoord, obstime: str, lat: float, lon: float, elev_m: float = 0.0) -> float:
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=elev_m*u.m)
    altaz = center.transform_to(AltAz(obstime=Time(obstime), location=loc))
    z = (90*u.deg - altaz.alt).to(u.deg).value
    if z >= 90: return np.inf
    return float(1/np.cos(np.deg2rad(z)))

def apply_extinction_correction(mag_lim: float, airmass: float, k_g: float = 0.15) -> float:
    if not np.isfinite(airmass) or airmass <= 0: return mag_lim
    return float(mag_lim + k_g*(airmass - 1.0))

def classify_bortle_like(glim: float) -> str:
    if glim is None or not np.isfinite(glim): return 'unknown'
    if glim >= 20.5: return 'Bortle ~2 (excellent)'
    if glim >= 19.5: return 'Bortle ~3'
    if glim >= 18.5: return 'Bortle ~4'
    if glim >= 17.5: return 'Bortle ~5'
    if glim >= 16.5: return 'Bortle ~6'
    if glim >= 15.5: return 'Bortle ~7'
    return 'Bortle ~8–9 (urban)'
