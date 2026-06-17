#!/usr/bin/env python3
"""
download_geo.py — fetch the WT 10x matrices declared in params/samples.yml.

Same job as download_geo.sh, but pure Python: it spawns no shells, so it
sidesteps the SHLVL problem entirely. URLs live only in samples.yml.

Idempotent (skips files already present, so re-running resumes a partial pull).

Run from the project root:
    python3 scripts/download_geo.py
"""
import sys, urllib.request, shutil
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("Need PyYAML:  conda install -c conda-forge pyyaml   (or pip install pyyaml)")

here    = Path(__file__).resolve().parent
roster  = here.parent /"phase1" /"params" / "samples.yml"
raw_dir = here.parent / "raw_data"

if not roster.exists():
    sys.exit(f"Roster not found: {roster}")

samples = yaml.safe_load(roster.read_text())["samples"]


def fetch(url, target, tries=3):
    if target.exists() and target.stat().st_size > 0:
        print(f"  [skip] {target.name}")
        return
    print(f"  [get ] {target.name}")
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r, open(target, "wb") as out:
                shutil.copyfileobj(r, out)
            return
        except Exception as e:
            if target.exists():
                target.unlink()
            if attempt == tries:
                sys.exit(f"  FAILED after {tries} tries: {url}\n  {e}")
            print(f"  retry {attempt}/{tries - 1} ...")


for s in samples:
    sid = s["id"]
    print(f"== {sid} ==")
    dest = raw_dir / sid
    dest.mkdir(parents=True, exist_ok=True)
    fetch(s["barcodes_url"], dest / "barcodes.tsv.gz")
    fetch(s["features_url"], dest / "features.tsv.gz")   # GEO names this *_genes.tsv.gz
    fetch(s["matrix_url"],   dest / "matrix.mtx.gz")

print(f"\nDone. Verify each sample has 3 files:\n  ls -lh {raw_dir}/*/")
