"""Regenerate the README sample figures using the server's own tools.

Usage:  MPLBACKEND=Agg .venv/bin/python scripts/make_readme_figures.py [--skip-network-heavy]
Writes PNGs to output/readme_figs/; copy the ones used in the README to docs/figures/.
"""
from __future__ import annotations

import json
import re
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.bao import (  # noqa: E402
    fetch_des_y6_bao_likelihood,
    fetch_des_y6_clustering,
    fit_bao_cosmology,
    plot_bao_contours,
    plot_bao_distances,
)
from tools.euclid import fetch_euclid_catalog, fetch_euclid_strong_lenses  # noqa: E402

OUT = ROOT / "output" / "readme_figs"



def _pngs(res) -> list[Path]:
    """Existing .png paths mentioned anywhere in a tool result."""
    d = res.model_dump() if hasattr(res, "model_dump") else res
    blob = json.dumps(d, default=str)
    return [Path(m) for m in re.findall(r'"([^"]+\.png)"', blob) if Path(m).exists()]


def run(name, fn, **kw):
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    try:
        res = fn(output_dir=str(out), **kw)
        pngs = _pngs(res) or sorted(out.glob("*.png"))
        print(f"[ok] {name}: {[p.name for p in pngs]}")
        return res, pngs
    except Exception:
        print(f"[fail] {name}")
        traceback.print_exc()
        return None, []


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    made = {}
    _, made["bao_distances"] = run("bao_distances", plot_bao_distances,
                                   datasets="desi_dr2,desi_dr1,boss_dr12_bao,sixdf,sdss_mgs")
    for key, ds in [("desi", ["desi_dr2"]), ("desi_desy6", ["desi_dr2", "des_y6_nodesi"])]:
        res, _ = run(f"fit_{key}", fit_bao_cosmology, datasets=ds, n_grid=80)
        grid = sorted((OUT / f"fit_{key}").glob("*grid*.npz"))
        if grid:
            _, made[f"contour_{key}"] = run(f"contour_{key}", plot_bao_contours, grid_file=str(grid[0]),
                                            dataset_label=" + ".join(ds))
    _, made["des_y6_alpha"] = run("des_y6_alpha", fetch_des_y6_bao_likelihood, variant="fid", plot=True)
    _, made["des_y6_acf"] = run("des_y6_acf", fetch_des_y6_clustering, statistic="acf", plot=True)
    if "--skip-network-heavy" not in sys.argv:
        _, made["euclid_spec"] = run("euclid_spec", fetch_euclid_catalog, sample="spectroscopic",
                                     ra=80.3, dec=-51.25, radius_deg=0.3)
        _, made["euclid_lenses"] = run("euclid_lenses", fetch_euclid_strong_lenses, plot=True)
    print(json.dumps({k: [str(p) for p in v] for k, v in made.items()}, indent=1))


if __name__ == "__main__":
    main()
