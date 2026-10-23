"""Check retained visualization coordinates with the original notebook's Kabsch fit."""
from pathlib import Path
import ast
import csv
import hashlib
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-17/reproducibility_handoff"
OUT.mkdir(parents=True, exist_ok=True)
notebook = json.loads((ROOT / "kaggle/notebooks/legacy/sf-nbench-full.ipynb").read_text())
functions = []
for cell in notebook["cells"]:
    if cell["cell_type"] != "code":
        continue
    tree = ast.parse("".join(cell["source"]))
    functions += [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "kabsch"]
assert len(functions) == 1
scope = {"np": np}
exec(compile(ast.Module(body=functions, type_ignores=[]), "original_notebook_kabsch", "exec"), scope)
kabsch = scope["kabsch"]
atoms = ["N", "CA", "C", "CB"]
aa = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
              "ARNDCQEGHILKMFPSTWYV"))


def score(pred, ref, regions):
    ok = np.isfinite(ref).all(-1) & np.isfinite(pred).all(-1)
    result = {"n_atoms_matched": int(ok.sum())}
    for mode, mask in (("chain", np.ones(len(ref), bool)), ("fw", regions["framework"])):
        selected = ok & mask[:, None]
        assert selected.sum() >= 3
        rotation, center, reference_center = kabsch(pred[selected], ref[selected])
        aligned = (pred - center) @ rotation + reference_center
        for name, region in regions.items():
            selected_region = ok & region[:, None]
            delta = aligned[selected_region] - ref[selected_region]
            result[f"{mode}_{name}"] = float(np.sqrt(np.mean(np.sum(delta ** 2, axis=-1))))
    return result


rng = np.random.default_rng(17)
reference = rng.normal(size=(20, 4, 3))
regions = {"cdr3": np.arange(20) >= 16, "framework": np.arange(20) < 16}
shifted = reference.copy()
shifted[regions["cdr3"]] += [5, 0, 0]
test = score(shifted, reference, regions)
assert abs(test["fw_framework"]) < 1e-12 and abs(test["fw_cdr3"] - 5) < 1e-12
rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
test_rotated = score(shifted @ rotation + [50, -30, 11], reference, regions)
assert all(abs(test[key] - test_rotated[key]) < 1e-10 for key in test)
targets = {t["pdb"] + t["H"]: t for t in json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())}
rows = []
for model, folder, suffix, report_name, score_key in (
    ("SimpleFold", "simplefold", "_pred.pdb", "viz_report.json", "rescored_fw_aligned"),
    ("IgFold", "igfold", "_igfold.pdb", "igfold_viz_report.json", "igfold_fw_aligned"),
    ("NanoBodyBuilder2", "nbb2", "_nbb2.pdb", "nbb2_viz_report.json", "nbb2_fw_aligned"),
):
    directory = ROOT / "outputs/visualization" / folder
    old_reports = {r["tid"]: r for r in json.loads((directory / report_name).read_text())}
    for path in sorted(directory.glob("*" + suffix)):
        tid = path.name.removesuffix(suffix)
        target = targets[tid]
        residues, coordinates = {}, {}
        for line in path.read_text().splitlines():
            if not line.startswith("ATOM") or line[12:16].strip() not in atoms:
                continue
            assert line[16] in (" ", "A"), (path.name, "unexpected alternate location")
            key = line[21:27]
            residues.setdefault(key, aa[line[17:20]])
            slot = (key, atoms.index(line[12:16].strip()))
            assert slot not in coordinates, (path.name, "duplicate atom", slot)
            coordinates[slot] = [float(line[i:i + 8]) for i in (30, 38, 46)]
        sequence = "".join(residues.values())
        offset = target["hseq"].find(sequence)
        assert offset >= 0 and target["hseq"].count(sequence) == 1
        index = {key: i + offset for i, key in enumerate(residues)}
        pred = np.full((len(target["hseq"]), 4, 3), np.nan)
        for (key, atom), xyz in coordinates.items():
            pred[index[key], atom] = xyz
        ref = np.array([[[np.nan] * 3 if a is None else a for a in residue] for residue in target["hca"]])
        masks = {}
        for region in ("cdr1", "cdr2", "cdr3"):
            lo, hi = target[region]
            masks[region] = (np.arange(len(ref)) >= lo) & (np.arange(len(ref)) < hi)
        masks["framework"] = ~(masks["cdr1"] | masks["cdr2"] | masks["cdr3"])
        scores = score(pred, ref, masks)
        old = old_reports[tid][score_key]
        max_difference = max(abs(scores["fw_" + k] - v) for k, v in old.items())
        assert max_difference < .011, (model, tid, max_difference)
        rows.append({"model": model, "tid": tid, "source": str(path.relative_to(ROOT)),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "maximum_difference_from_two_decimal_viz_report": max_difference, **scores})
assert len(rows) == 18
with (OUT / "retained_visualization_rescoring.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
report = {"structures": len(rows), "source": "six separately generated visualization examples per model",
          "original_notebook_kabsch_function_reused": True,
          "maximum_absolute_difference_from_two_decimal_reports_A": max(r["maximum_difference_from_two_decimal_viz_report"] for r in rows),
          "synthetic_region_mask_and_rigid_rotation_tests": "PASS",
          "meaning": "saved-coordinate to regional-score reproduction only; not original full-benchmark prediction recovery",
          "exact_prepared_7R73_inputs_available": False}
(OUT / "retained_visualization_rescoring.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
