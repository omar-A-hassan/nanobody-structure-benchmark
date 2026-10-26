"""Score the four supplied prepared 7R73 VHH inputs against the frozen reference."""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import io
import json

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from Bio.SVDSuperimposer import SVDSuperimposer

ROOT = Path(__file__).resolve().parents[1]
ATOMS = ("N", "CA", "C", "CB")
INPUTS = {
    "NanoBodyBuilder2": "7R73A_NBB2_DOCK.pdb",
    "IgFold": "7R73A_IgFold_DOCK.pdb",
    "SimpleFold": "7R73A_SimpleFold_DOCK.pdb",
    "Experimental": "7R73_experimental_DOCK.pdb",
}
NOTEBOOK = ROOT / "kaggle/notebooks/legacy/sf-nbench-full.ipynb"
TARGETS = ROOT / "data/inputs/targets_vhh.json"
OUTPUT = ROOT / "outputs/docking/7R73_prepared_regional_rmsd"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def original_kabsch():
    functions = []
    for cell in json.loads(NOTEBOOK.read_text())["cells"]:
        if cell["cell_type"] == "code":
            tree = ast.parse("".join(cell["source"]))
            functions.extend(node for node in tree.body
                             if isinstance(node, ast.FunctionDef) and node.name == "kabsch")
    assert len(functions) == 1, "Expected one original Kabsch function"
    scope = {"np": np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), "original_kabsch", "exec"), scope)
    return scope["kabsch"]


def read_vhh(path, expected_sequence):
    residues, coordinates = {}, {}
    heavy_atoms = 0
    seen_atoms = set()
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        assert line[21] == "H" and line[16] == " ", "Unexpected chain or alternate location"
        key = line[21:27]
        aa = protein_letters_3to1[line[17:20]]
        assert residues.setdefault(key, aa) == aa, "Inconsistent residue identity"
        atom = line[12:16].strip()
        assert (key, atom) not in seen_atoms, "Duplicate atom"
        seen_atoms.add((key, atom))
        assert line[76:78].strip() not in {"H", "D"}, "Expected heavy atoms only"
        heavy_atoms += 1
        if atom in ATOMS:
            xyz = [float(line[i:i + 8]) for i in (30, 38, 46)]
            assert np.isfinite(xyz).all(), "Non-finite coordinate"
            coordinates[key, ATOMS.index(atom)] = xyz
    assert "".join(residues.values()) == expected_sequence, "VHH sequence/order mismatch"
    index = {key: i for i, key in enumerate(residues)}
    result = np.full((len(expected_sequence), len(ATOMS), 3), np.nan)
    for (key, slot), xyz in coordinates.items():
        result[index[key], slot] = xyz
    return result, heavy_atoms


def region_masks(target):
    index = np.arange(len(target["hseq"]))
    masks = {name: (index >= target[name][0]) & (index < target[name][1])
             for name in ("cdr1", "cdr2", "cdr3")}
    masks["framework"] = ~(masks["cdr1"] | masks["cdr2"] | masks["cdr3"])
    return masks


def framework_scores(prediction, reference, masks, kabsch):
    matched = np.isfinite(reference).all(-1) & np.isfinite(prediction).all(-1)
    fit = matched & masks["framework"][:, None]
    assert fit.sum() >= 3 and np.linalg.matrix_rank(reference[fit] - reference[fit].mean(0)) >= 2
    rotation, center, reference_center = kabsch(prediction[fit], reference[fit])
    aligned = (prediction - center) @ rotation + reference_center
    check = SVDSuperimposer()
    check.set(reference[fit], prediction[fit])
    check.run()
    other_rotation, translation = check.get_rotran()
    independently_aligned = prediction @ other_rotation + translation
    maximum_difference = float(np.max(abs(aligned[matched] - independently_aligned[matched])))
    assert maximum_difference < 1e-9, "Independent alignment check failed"
    scores = {}
    for name, mask in masks.items():
        selected = matched & mask[:, None]
        assert selected.sum() >= 3, f"Insufficient matched atoms: {name}"
        scores[name] = float(np.sqrt(np.mean(np.sum((aligned[selected] - reference[selected]) ** 2, axis=-1))))
    return scores, int(matched.sum()), int(fit.sum()), maximum_difference


def calculate():
    target = next(t for t in json.loads(TARGETS.read_text()) if t["pdb"] + t["H"] == "7R73A")
    reference = np.array([[[np.nan] * 3 if atom is None else atom for atom in residue]
                          for residue in target["hca"]], dtype=float)
    masks = region_masks(target)
    kabsch = original_kabsch()
    rows = []
    for model, filename in INPUTS.items():
        path = ROOT / "data/docking" / filename
        prediction, heavy_atoms = read_vhh(path, target["hseq"])
        scores, matched, fit_atoms, difference = framework_scores(prediction, reference, masks, kabsch)
        rows.append({"model": model, "file": str(path.relative_to(ROOT)), "sha256": sha(path),
                     "residues": len(prediction), "heavy_atoms": heavy_atoms,
                     "matched_atoms": matched, "framework_fit_atoms": fit_atoms,
                     **scores, "independent_alignment_max_difference_A": difference})
    return {"target": "7R73A", "alignment": "framework", "atom_names": list(ATOMS),
            "cdr_spans_zero_based_half_open": {name: target[name] for name in ("cdr1", "cdr2", "cdr3")},
            "reference": str(TARGETS.relative_to(ROOT)), "reference_sha256": sha(TARGETS),
            "kabsch_source": str(NOTEBOOK.relative_to(ROOT)), "kabsch_source_sha256": sha(NOTEBOOK),
            "scoring_source_sha256": sha(Path(__file__)),
            "independent_check": "BioPython SVDSuperimposer; coordinate tolerance 1e-9 angstrom",
            "scope": "Prepared monomer regional RMSD only; no DockQ or epitope rescore and no HDOCK job-file identity check",
            "rows": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare with saved scores without replacing them")
    args = parser.parse_args()
    report = calculate()
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(report["rows"][0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(report["rows"])
    if args.check:
        saved = json.loads(OUTPUT.with_suffix(".json").read_text())
        assert {k: v for k, v in saved.items() if k != "rows"} == {k: v for k, v in report.items() if k != "rows"}
        assert len(saved["rows"]) == len(report["rows"])
        numeric = {"cdr1", "cdr2", "cdr3", "framework", "independent_alignment_max_difference_A"}
        for old, new in zip(saved["rows"], report["rows"]):
            assert old.keys() == new.keys()
            for key in old:
                assert (np.isclose(old[key], new[key], atol=1e-10, rtol=0) if key in numeric else old[key] == new[key]), key
        saved_csv = io.StringIO(newline="")
        csv_writer = csv.DictWriter(saved_csv, fieldnames=list(saved["rows"][0]), lineterminator="\n")
        csv_writer.writeheader()
        csv_writer.writerows(saved["rows"])
        assert OUTPUT.with_suffix(".csv").read_bytes() == saved_csv.getvalue().encode(), "Saved CSV/JSON disagree"
    else:
        OUTPUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
        OUTPUT.with_suffix(".csv").write_bytes(stream.getvalue().encode())
    print(json.dumps({"prepared_inputs": len(report["rows"]), "saved_score_check": args.check,
                      "rows": report["rows"]}, indent=2))


if __name__ == "__main__":
    main()
