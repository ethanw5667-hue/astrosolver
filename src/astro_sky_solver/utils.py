from pathlib import Path

def log(m): print(f"[astro-sky-solver] {m}")

def ensure_dir_for(p: Path): p.parent.mkdir(parents=True, exist_ok=True)
