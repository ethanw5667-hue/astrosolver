# Example usage for astro-sky-solver
# Ensure dependencies are installed and ASTROMETRY_API_KEY is set (for web solving).
# Then run:
#   python -m astro_sky_solver.cli samples/astropic.jpg --out-prefix out/result
# Or for a solved FITS:
#   python -m astro_sky_solver.cli samples/astropic_wcs.fits --out-prefix out/result
from astro_sky_solver.cli import main

if __name__ == "__main__":
    main()
