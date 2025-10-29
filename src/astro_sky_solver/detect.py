import numpy as np
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from astropy.table import Table

def estimate_background(data, box_size=64, filter_size=3, sigma_clip=3.0):
    sc = SigmaClip(sigma=sigma_clip)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, box_size=box_size, filter_size=filter_size,
                       sigma_clip=sc, bkg_estimator=bkg_estimator)
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
        return Table(names=["xcentroid","ycentroid","flux"], dtype=[float,float,float])
    tbl.sort("flux"); tbl.reverse()
    return tbl[["xcentroid","ycentroid","flux"]]
