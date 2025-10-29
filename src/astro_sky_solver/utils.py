from pathlib import Path

def log(msg: str) -> None:
    print(f"[astro-sky-solver] {msg}", flush=True)

def ensure_dir_for(prefix_path: Path) -> None:
    prefix_path.parent.mkdir(parents=True, exist_ok=True)
