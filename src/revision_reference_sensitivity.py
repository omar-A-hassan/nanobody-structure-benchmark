"""Saved-score sensitivity excluding references changed by alternate-atom selection."""
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP
import argparse
import csv
import hashlib
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-17/reference_sensitivity"
MODELS = ["SimpleFold", "IgFold", "NanoBodyBuilder2"]
PATHS = ["outputs/structure/simplefold/metrics/per_target_scores.json",
         "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
         "outputs/structure/nbb2/per_target_scores_nbb2.json"]
REGIONS = ["cdr1", "cdr2", "cdr3", "framework"]
BOOTS, SEED = 4000, 20261017


def rounded(value):
    return str(Decimal(str(round(value, 10))).quantize(Decimal(".001"), rounding=ROUND_HALF_UP))


def paired_intervals(deltas, blocks):
    rng = np.random.default_rng(SEED)
    draws = np.empty((BOOTS, deltas.shape[1]))
    for b in range(BOOTS):
        selected = np.concatenate([blocks[j] for j in rng.integers(len(blocks), size=len(blocks))])
        draws[b] = np.median(deltas[selected], axis=0)
    return np.quantile(draws, [.025, .975], axis=0).T.tolist()


assert rounded(4.3645) == "4.365"
assert paired_intervals(np.ones((4, 1)), [np.array([0, 1]), np.array([2, 3])]) == [[1., 1.]]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--expanded-cohort", type=Path, help="B2 final expanded sensitivity CSV; must contain tid")
args = parser.parse_args()
source = ROOT / "outputs/revision_2026-10-17/reference_reconstruction/summary.json"
reconstruction = json.loads(source.read_text())
assert not reconstruction["provisional"] and reconstruction["archive_md5_verified"]
assert reconstruction["frozen_records_checked"] == reconstruction["frozen_exact_matches"] == 597
changed = set(reconstruction["reference_policy_impact"]["changed_tids"])
targets = {t["pdb"] + t["H"]: t for t in json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())}
scores = {model: {r["tid"]: r for r in json.loads((ROOT / path).read_text())}
          for model, path in zip(MODELS, PATHS)}
paired = sorted(set.intersection(*(set(s) for s in scores.values())))
assert len(paired) == 588
old_audit = json.loads((ROOT / "outputs/revision_2026-10-16/training_overlap_audit.json").read_text())
old_filtered = sorted(old_audit["common_lower_risk_subset"]["tids"])
assert len(old_filtered) == 398
cohorts = {"original_paired_588": paired, "previous_lower_identified_risk_398": old_filtered}
if args.expanded_cohort:
    with args.expanded_cohort.open() as stream:
        expanded = sorted(r["tid"] for r in csv.DictReader(stream))
    assert len(expanded) == len(set(expanded)) and set(expanded) <= set(paired)
    cohorts["expanded_lower_identified_risk"] = expanded
for name, tids in list(cohorts.items()):
    cohorts[name + "_without_altloc_policy_changes"] = [tid for tid in tids if tid not in changed]

report = {"scope": "Sensitivity of historical saved scores, NOT rescoring predictions against changed references",
          "bootstrap": {"replicates": BOOTS, "seed": SEED, "method": "paired PDB-block percentile; exploratory, no multiplicity correction"},
          "reference_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
          "changed_reference_targets": sorted(changed), "cohorts": {}}
summary_rows, member_rows = [], []
for name, tids in cohorts.items():
    assert tids
    pdbs = sorted({targets[tid]["pdb"] for tid in tids})
    blocks = [np.array([i for i, tid in enumerate(tids) if targets[tid]["pdb"] == pdb]) for pdb in pdbs]
    block = {"n": len(tids), "pdb_count": len(pdbs), "tids": tids,
             "changed_reference_members": sorted(set(tids) & changed), "models": {}, "paired_h3_effects": {}}
    for model in MODELS:
        block["models"][model] = {}
        for mode in ("chain", "fw"):
            for region in REGIONS:
                key = f"{mode}_{region}"
                values = np.array([scores[model][tid][key] for tid in tids])
                result = {"median": float(np.median(values)), "display": rounded(np.median(values)),
                          "mean": float(values.mean()), "sd": float(values.std(ddof=1)),
                          "above6_count": int((values > 6).sum()), "above6_pct": float(100 * (values > 6).mean())}
                block["models"][model][key] = result
                summary_rows.append({"cohort": name, "n": len(tids), "model": model, "metric": key, **result})
    keys, columns = [], []
    for mode in ("chain", "fw"):
        for comparator in MODELS[1:]:
            keys.append(f"SimpleFold minus {comparator}: {mode}_cdr3")
            columns.append([scores["SimpleFold"][tid][f"{mode}_cdr3"] - scores[comparator][tid][f"{mode}_cdr3"] for tid in tids])
    deltas = np.array(columns).T
    intervals = paired_intervals(deltas, blocks)
    for key, values, interval in zip(keys, deltas.T, intervals):
        block["paired_h3_effects"][key] = {"median": float(np.median(values)), "ci95": interval}
    report["cohorts"][name] = block
    member_rows += [{"cohort": name, "tid": tid, "reference_selection_changed": tid in changed} for tid in tids]
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "reference_sensitivity.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
for filename, rows in (("reference_sensitivity_summary.csv", summary_rows), ("cohort_members.csv", member_rows)):
    with (OUT / filename).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
print(json.dumps({name: {"n": row["n"], "chain_h3_medians": {m: row["models"][m]["chain_cdr3"]["display"] for m in MODELS}}
                  for name, row in report["cohorts"].items()}, indent=2))
