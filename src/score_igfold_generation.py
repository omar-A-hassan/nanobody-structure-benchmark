"""Score the retained IgFold generation coordinates and compare them with the primary record.

The scoring procedure is the one used for the SimpleFold and NanoBodyBuilder2 arms: matched
N/CA/C/CB atoms, IMGT regions, an unweighted proper-rotation Kabsch fit, whole-chain and
framework fits, and sequence-identity residue correspondence rather than index arithmetic.
"""
from pathlib import Path
from statistics import median
import argparse
import csv
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "outputs/structure/igfold/metrics/per_target_scores_igfold.json"
AA3to1 = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
                  "ARNDCQEGHILKMFPSTWYV"))
ATOMS = {"N": 0, "CA": 1, "C": 2, "CB": 3}
REGIONS = ["cdr1", "cdr2", "cdr3", "framework"]
FIELDS = [f"{fit}_{region}" for fit in ("chain", "fw") for region in REGIONS]


def kabsch(P, Q):
    Pm, Qm = P.mean(0), Q.mean(0)
    V, _, W = np.linalg.svd((P - Pm).T @ (Q - Qm))
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ W))])
    return V @ D @ W, Pm, Qm


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(-1).mean()))


def read_backbone(path):
    """Return ordered residue keys, one-letter sequence and the matched-atom coordinates."""
    seen, order, records = {}, [], []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in ATOMS:
            continue
        key = (line[21], line[22:26].strip(), line[26].strip())
        if key not in seen:
            seen[key] = len(order)
            order.append(AA3to1.get(line[17:20].strip(), "X"))
        records.append((seen[key], ATOMS[name],
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return "".join(order), records


def score_target(target, path):
    pseq, records = read_backbone(path)
    offset = target["hseq"].find(pseq)
    if offset < 0:
        return None, f"sequence not found: predicted {len(pseq)} vs reference {len(target['hseq'])}"
    length = len(target["hseq"])
    pred = np.full((length, 4, 3), np.nan)
    for residue, atom, xyz in records:
        pred[residue + offset, atom] = xyz
    reference = np.array([[[np.nan] * 3 if v is None else v for v in r] for r in target["hca"]], float)

    regions = {}
    for key in ("cdr1", "cdr2", "cdr3"):
        mask = np.zeros(length, bool)
        low, high = target[key]
        mask[low:high] = True
        regions[key] = mask
    regions["framework"] = ~(regions["cdr1"] | regions["cdr2"] | regions["cdr3"])

    ok = np.isfinite(reference[..., 0]) & np.isfinite(pred[..., 0])
    if ok[:, 1].sum() < 60:
        return None, "insufficient matched CA coverage"

    row = {"instance": target["instance"], "tid": target["pdb"] + target["H"], "L": length,
           "cdr3_len": int(target["cdr3"][1] - target["cdr3"][0]),
           "n_atoms_matched": int(ok.sum()), "seq_offset": int(offset), "pred_len": len(pseq),
           "method": target["method"], "holo": target["holo"], "ab_split": target["ab_split"],
           "resolution": target["resolution"]}
    for mode, fit in (("chain", np.ones(length, bool)), ("fw", regions["framework"])):
        selected = ok & fit[:, None]
        if selected.sum() < 60:
            continue
        rotation, center, reference_center = kabsch(pred[selected], reference[selected])
        aligned = np.full_like(pred, np.nan)
        aligned[ok] = (pred[ok] - center) @ rotation + reference_center
        for region in REGIONS:
            chosen = ok & regions[region][:, None]
            row[f"{mode}_{region}"] = round(rmsd(aligned[chosen], reference[chosen]), 3) if chosen.sum() >= 3 else None
    return row, None


def self_test():
    """A perfect framework with a 5 A displaced loop must read 0.000 and 5.000."""
    rng = np.random.default_rng(17)
    reference = rng.normal(size=(20, 4, 3))
    pred = reference.copy()
    loop = np.arange(20) >= 16
    pred[loop] += [5, 0, 0]
    ok = np.ones((20, 4), bool)
    fit = ~loop
    rotation, center, reference_center = kabsch(pred[ok & fit[:, None]], reference[ok & fit[:, None]])
    aligned = (pred - center) @ rotation + reference_center
    assert abs(rmsd(aligned[ok & fit[:, None]], reference[ok & fit[:, None]])) < 1e-9
    assert abs(rmsd(aligned[ok & loop[:, None]], reference[ok & loop[:, None]]) - 5) < 1e-9
    rotated = pred @ np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]]) + [50, -30, 11]
    rotation, center, reference_center = kabsch(rotated[ok & fit[:, None]], reference[ok & fit[:, None]])
    moved = (rotated - center) @ rotation + reference_center
    assert abs(rmsd(moved[ok & loop[:, None]], reference[ok & loop[:, None]]) - 5) < 1e-8


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", required=True, help="directory of per-target PDB files")
    parser.add_argument("--out", required=True, help="output directory")
    args = parser.parse_args()
    self_test()
    print("scorer self-test passed")

    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    coordinates = Path(args.coordinates)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows, problems = [], {}
    for target in targets:
        tid = target["pdb"] + target["H"]
        path = coordinates / f"{tid}.pdb"
        if not path.exists():
            problems[tid] = "missing coordinate file"
            continue
        row, problem = score_target(target, path)
        if problem:
            problems[tid] = problem
        else:
            rows.append(row)
    print(f"scored {len(rows)} / {len(targets)}")
    if problems:
        print("problems:", problems)

    violations = [r["tid"] for r in rows if r["seq_offset"] != 0 or r["pred_len"] != r["L"]]
    print(f"residue-preservation violations: {len(violations)}" + (f" {violations[:5]}" if violations else " (none)"))

    columns = (["instance", "tid", "L", "cdr3_len", "n_atoms_matched", "seq_offset", "pred_len",
                "method", "holo", "ab_split", "resolution"] + FIELDS)
    (out / "per_target_scores_igfold.json").write_text(json.dumps(rows, indent=1) + "\n")
    with (out / "per_target_scores_igfold.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)

    primary = {r["tid"]: r for r in json.loads(PRIMARY.read_text())}
    generated = {r["tid"]: r for r in rows}
    common = sorted(set(primary) & set(generated))
    comparison = {"n_primary": len(primary), "n_generated": len(generated), "n_common": len(common),
                  "scoring_problems": problems, "residue_preservation_violations": violations,
                  "fields": {}}
    with (out / "comparison_vs_primary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["tid", "metric", "primary_score", "generated_score", "generated_minus_primary"])
        for tid in common:
            for field in FIELDS:
                a, b = primary[tid][field], generated[tid][field]
                writer.writerow([tid, field, a, b, round(b - a, 3)])
    for field in FIELDS:
        a = [primary[t][field] for t in common]
        b = [generated[t][field] for t in common]
        changes = [y - x for x, y in zip(a, b)]
        comparison["fields"][field] = {
            "primary_median": round(median(a), 3), "generated_median": round(median(b), 3),
            "median_paired_change": round(median(changes), 3),
            "median_absolute_change": round(median(abs(d) for d in changes), 3),
            "maximum_absolute_change": round(max(abs(d) for d in changes), 3),
            "exact_values": sum(abs(d) < 5e-4 for d in changes), "n": len(common)}
    (out / "comparison_vs_primary.json").write_text(json.dumps(comparison, indent=2) + "\n")

    print(f"\n{'metric':<18}{'primary':>9}{'generated':>11}{'med |d|':>9}{'max |d|':>9}{'exact':>7}")
    for field in FIELDS:
        c = comparison["fields"][field]
        print(f"{field:<18}{c['primary_median']:>9.3f}{c['generated_median']:>11.3f}"
              f"{c['median_absolute_change']:>9.3f}{c['maximum_absolute_change']:>9.3f}{c['exact_values']:>7}")


if __name__ == "__main__":
    main()
