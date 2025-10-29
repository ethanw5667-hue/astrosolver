# run_example.py — simple demo runner for astro-sky-solver
# Make sure your environment is set up:
#   export ASTROMETRY_API_KEY=YOUR_KEY
# Then run:
#   python run_example.py

import subprocess
import os

# Example usage
cmd = [
    "python", "app.py",
    "samples/astropic.jpg",
    "--out-prefix", "out/example",
    "--lat", "43.65",
    "--lon", "-79.38",
    "--time", "2025-07-01 02:30:00"
]

print("Running:", " ".join(cmd))
subprocess.run(cmd)
