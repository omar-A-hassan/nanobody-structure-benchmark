"""Recreate paper Figures 3–5 from saved data; never run the superseded bar plot."""
from pathlib import Path
import hashlib
import json
import runpy

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "paper_figure_sources/make_figures.py"
assert source.is_file(), "The handoff must contain the original plot source"
out = ROOT / "outputs/revision_2026-10-17/reused_figures"
out.mkdir(parents=True, exist_ok=True)
namespace = runpy.run_path(str(source))
names = ("figure2_uncertainty", "figure3_prior_coverage", "figure4_esm_probe")
context = namespace[names[0]].__globals__
assert all(namespace[name].__globals__ is context for name in names)
context.update(DATA=ROOT / "outputs", OUT=out,
               spearman=lambda x, y: float(spearmanr(x, y).statistic))
# The historical plotting helper used ordinal ranks; the revised definition
# requires average ranks for ties. No model, data, axes or styling is changed.
assert abs(context["spearman"]([1, 1, 2], [1, 2, 3]) - np.sqrt(3) / 2) < 1e-12
for name in names:
    namespace[name]()
for stem in names:
    for suffix in (".png", ".pdf"):
        assert (out / (stem + suffix)).is_file()
rows = json.loads((ROOT / "outputs/trajectory/full-597/traj_results.json").read_text())["rows"]
error = [r["err"] for r in rows]
rho = {key: context["spearman"]([r[key] for r in rows], error)
       for key in ("spread_cdr3", "adiv_q1", "adiv_half", "adiv_final", "div_final")}
baseline = json.loads((ROOT / "outputs/revision_2026-10-16/revision_statistics.json").read_text())
assert len(rows) == 597 and abs(rho["spread_cdr3"] - baseline["uncertainty"]["trajectory"]["spread_rho"]) < 1e-12
report = {"status": "PASS", "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
          "functions": list(names), "n_trajectory": len(rows), "spearman_average_ties": rho,
          "outputs": [name + suffix for name in names for suffix in (".png", ".pdf")],
          "superseded_benchmark_bar_not_run": True}
(out / "reproduction.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
