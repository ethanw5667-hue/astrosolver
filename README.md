# astro-sky-solver

Detect stars in an image, plate-solve to figure out **which part of the sky** you're looking at,
identify **which stars**, and estimate **limiting magnitude** + a rough **light pollution / sky quality** index
from your detection depth.

## Quick Start

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Recommended: web plate solver (free key at nova.astrometry.net)
export ASTROMETRY_API_KEY=YOUR_KEY
```

### Run (JPG/PNG -> auto-solve on web)
```bash
python -m astro_sky_solver.cli samples/astropic.jpg --out-prefix out/result
```

### Run (already-solved FITS)
```bash
python -m astro_sky_solver.cli samples/astropic_wcs.fits --out-prefix out/result
```

Outputs: `*_matches.csv`, `*_annotated.png`, and console with limiting magnitude + rough sky class.

See `src/astro_sky_solver/extinction.py` for how limiting magnitude is derived.


## CLI

Run via module:

```bash
python -m astro_sky_solver.cli <image_path> --out-prefix out/result
```
