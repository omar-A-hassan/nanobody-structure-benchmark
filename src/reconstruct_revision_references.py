"""B3 reconstruction from SAbDab2 coordinates; never edits frozen/shared inputs.

Run with isolated uv: uv run --no-project --python 3.12 --with biopython==1.85 --with numpy==2.2.6 python
src/reconstruct_revision_references.py --archive PATH
Use --partial while the archive is downloading; its outputs are provisional.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import io
import gzip
import json
import platform
import shutil
import sys
import tarfile
import traceback
from urllib.request import Request, urlopen

import Bio
import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-17/reference_reconstruction"
CACHE = ROOT / "outputs/revision_2026-10-16/sources"
FROZEN = ROOT / "data/inputs/targets_vhh.json"
ATOMS = ["N", "CA", "C", "CB"]
REGIONS = [(27, 38), (56, 65), (105, 117)]
ARCHIVE_URL = "https://zenodo.org/records/20083995/files/splits.tar.gz?download=1"
ARCHIVE_MD5 = "0dbb4cc499e9eb77f14008b232f2c38c"
ARCHIVE_BYTES = 876381859


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def checksums(path):
    hashes = {k: hashlib.new(k) for k in ("md5", "sha256")}
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            for h in hashes.values():
                h.update(block)
    return {k: h.hexdigest() for k, h in hashes.items()}


def norm(value):
    return "" if value in (".", "?", " ") else value.strip()


def parse_atoms(payload, chain, occupancy_aware=True):
    """First deposited model; explicit replay or highest-occupancy atom selection."""
    d = MMCIF2Dict(io.StringIO(payload.decode()))
    cols = {k: d["_atom_site." + k] for k in (
        "group_PDB", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
        "label_comp_id", "label_atom_id", "label_alt_id", "occupancy",
        "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")}
    model = cols["pdbx_PDB_model_num"][0]
    label_seq = d.get("_atom_site.label_seq_id", ["."]*len(cols["group_PDB"]))
    residues, alternates = {}, []
    for i, c in enumerate(cols["auth_asym_id"]):
        polymer = cols["group_PDB"][i] == "ATOM" or (cols["group_PDB"][i] == "HETATM" and label_seq[i] not in (".", "?"))
        if c != chain or not polymer or cols["pdbx_PDB_model_num"][i] != model:
            continue
        key = (int(cols["auth_seq_id"][i]), norm(cols["pdbx_PDB_ins_code"][i]))
        aa = protein_letters_3to1.get(cols["label_comp_id"][i], "X")
        res = residues.setdefault(key, {"aa": aa, "resname": cols["label_comp_id"][i], "atoms": {}})
        if res["aa"] != aa:
            raise ValueError(f"disordered residue identity at {key}: {res['aa']} / {aa}")
        name = cols["label_atom_id"][i]
        if name not in ATOMS:
            continue
        candidate = {"occupancy": float(cols["occupancy"][i]),
                     "altloc": norm(cols["label_alt_id"][i]),
                     "coord": [float(cols[k][i]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]}
        previous = res["atoms"].get(name)
        if not occupancy_aware or previous is None or candidate["occupancy"] > previous["occupancy"]:
            res["atoms"][name] = candidate
        if candidate["altloc"]:
            alternates.append({"imgt": list(key), "atom": name, **candidate})
    return residues, alternates, model


def reconstruct(row, payload, occupancy_aware=True):
    numbering = json.loads(row["VH_numbering_list"])
    seq = row["VH_numerable_seq"]
    assert "".join(aa for aa, _ in numbering) == seq
    residues, alternates, model = parse_atoms(payload, row["Hchain"], occupancy_aware)
    target = {"instance": row["INSTANCE"], "pdb": row["PDB_ID"][-4:].upper(),
              "H": row["Hchain"], "sabdab_id": row["SABDAB_ID"],
              "cdrh3_cluster": row["cdrh3_cluster"], "hseq": seq,
              "numbering": "IMGT", "atoms": ATOMS.copy()}
    issues = []
    for j, (low, high) in enumerate(REGIONS, 1):
        positions = [i for i, (_, (n, _)) in enumerate(numbering) if low <= n <= high]
        span = [positions[0], positions[-1] + 1] if positions else None
        target[f"cdr{j}"] = span
        target[f"CDRH{j}"] = row[f"CDRH{j}"]
        if not positions:
            issues.append({"condition": "empty_cdr", "region": j})
        elif positions != list(range(*span)):
            issues.append({"condition": "noncontiguous_cdr", "region": j, "positions": positions})
        if span is not None and seq[slice(*span)] != row[f"CDRH{j}"]:
            issues.append({"condition": "cdr_string_mismatch", "region": j})
    # This is a new explicit domain-coverage guard, not a recovered historical threshold.
    framework_coverage = {f"fr{j}": sum(low <= n <= high for _, (n, _) in numbering)
                          for j, (low, high) in enumerate(((1, 26), (39, 55), (66, 104), (118, 128)), 1)}
    if any(n == 0 for n in framework_coverage.values()):
        issues.append({"condition": "missing_framework_segment", "coverage": framework_coverage})
    target.update({k: row[k] for k in ("resolution", "method", "holo", "ab_split", "PDBdepo")})
    coords, selected = [], []
    for i, (aa, (n, insertion)) in enumerate(numbering):
        key = (n, norm(insertion))
        res = residues.get(key)
        if res is None:
            issues.append({"condition": "missing_numbered_residue", "index": i, "imgt": list(key)})
        elif res["aa"] != aa:
            issues.append({"condition": "residue_identity_mismatch", "index": i,
                           "imgt": list(key), "expected": aa, "observed": res["aa"], "source_resname": res["resname"]})
        quad = []
        for atom in ATOMS:
            a = res["atoms"].get(atom) if res else None
            quad.append([round(v, 3) for v in a["coord"]] if a else None)
            if a and a["altloc"]:
                selected.append({"index": i, "imgt": list(key), "atom": atom, **a})
        coords.append(quad)
        if any(a is None for a in quad[:3]):
            issues.append({"condition": "incomplete_backbone", "index": i, "imgt": list(key)})
        if aa == "G" and quad[3] is not None:
            issues.append({"condition": "glycine_has_cb", "index": i})
    target["n_resolved"] = sum(q[1] is not None for q in coords)
    target["hca"] = coords
    evidence = {"issues": issues, "source_model": model, "source_residues": len(residues),
                "numbered_residues": len(numbering), "alternate_atom_candidates": alternates,
                "framework_numbered_residue_counts": framework_coverage,
                "selected_alternate_atoms": selected,
                "source_residue_keys": [list(k) for k in residues]}
    return target, evidence


def compare(target, frozen):
    differences = [k for k in frozen if k != "hca" and target.get(k) != frozen[k]]
    old, new = frozen["hca"], target["hca"]
    null_mismatches, unequal_components, components, max_delta = 0, 0, 0, 0.
    mismatches = []
    assert len(old) == len(new)
    for i, (a, b) in enumerate(zip(old, new)):
        for j, (x, y) in enumerate(zip(a, b)):
            if (x is None) != (y is None):
                null_mismatches += 1
                mismatches.append({"residue_index": i, "atom": ATOMS[j], "frozen": x, "reconstructed": y})
            elif x is not None:
                components += 3
                unequal_components += sum(v != w for v, w in zip(x, y))
                max_delta = max(max_delta, max(abs(v-w) for v, w in zip(x, y)))
                if x != y:
                    mismatches.append({"residue_index": i, "atom": ATOMS[j], "frozen": x, "reconstructed": y})
    return {"metadata_differences": differences, "null_mask_mismatches": null_mismatches,
            "sequence_match": target["hseq"] == frozen["hseq"],
            "cdr_spans_match": all(target[f"cdr{j}"] == frozen[f"cdr{j}"] for j in (1, 2, 3)),
            "coordinate_components_compared": components, "unequal_coordinate_components": unequal_components,
            "max_abs_coordinate_difference_angstrom": max_delta, "coordinate_mismatches": mismatches,
            "exact_match": not differences and not null_mismatches and not unequal_components}


def fit_coordinates(p, q, fit_mask):
    a, b = p[fit_mask], q[fit_mask]
    pc, qc = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a-pc).T @ (b-qc))
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    return (p-pc) @ rotation+qc


def reference_displacement(old, new):
    def array(t):
        return np.array([[[np.nan]*3 if a is None else a for a in q] for q in t["hca"]])

    p, q = array(new), array(old)
    matched = np.isfinite(p).all(-1) & np.isfinite(q).all(-1)
    changed = matched & np.any(p != q, axis=-1)
    regions = {}
    for j in (1, 2, 3):
        mask = np.zeros(len(p), dtype=bool)
        mask[slice(*old[f"cdr{j}"])] = True
        regions[f"cdr{j}"] = mask
    regions["framework"] = ~np.logical_or.reduce(list(regions.values()))
    regions["full"] = np.ones(len(p), dtype=bool)
    result = {"tid": old["pdb"]+old["H"], "instance": old["instance"],
              "changed_atom_slots": int(changed.sum()), "changed_residues": int(changed.any(-1).sum()),
              "raw_max_atom_displacement": float(np.linalg.norm(p[changed]-q[changed], axis=-1).max()),
              "raw_full_ca_rmsd": float(np.sqrt(np.mean(np.sum((p[:, 1]-q[:, 1])**2, axis=-1))))}
    for region, mask in regions.items():
        result[f"changed_atom_slots_{region}"] = int((changed & mask[:, None]).sum())
    for mode in ("chain", "fw"):
        aligned = fit_coordinates(p, q, matched if mode == "chain" else matched & regions["framework"][:, None])
        squared = np.sum((aligned-q)**2, axis=-1)
        for region, mask in regions.items():
            use = matched & mask[:, None]
            result[f"{mode}_{region}_nca_c_cb_rmsd"] = float(np.sqrt(squared[use].mean()))
            ca = use[:, 1]
            result[f"{mode}_{region}_ca_rmsd"] = float(np.sqrt(squared[ca, 1].mean()))
        for j, atom in enumerate(ATOMS):
            result[f"{mode}_full_{atom}_rmsd"] = float(np.sqrt(squared[matched[:, j], j].mean()))
        result[f"{mode}_max_atom_displacement"] = float(np.sqrt(squared[matched].max()))
    return result


def csvout(path, rows):
    if rows:
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def impact_analysis(records, corrected, frozen, prefix):
    models = {name: {r["tid"]: r for r in json.loads((ROOT / path).read_text())}
              for name, path in {
                  "SimpleFold": "outputs/structure/simplefold/metrics/per_target_scores.json",
                  "IgFold": "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
                  "NanoBodyBuilder2": "outputs/structure/nbb2/per_target_scores_nbb2.json"}.items()}
    cohorts = {name: {r["tid"] for r in csv.DictReader((CACHE.parent / file).open())}
               for name, file in {
                   "lower_risk_398": "common_lower_risk_subset.csv",
                   "uncertainty_36": "uncertainty_groups_36.csv",
                   "coverage_40": "coverage_groups_40.csv"}.items()}
    paired = set.intersection(*(set(m) for m in models.values()))
    cohorts["paired_588"] = paired
    cohorts["flagged_partial_or_obsolete"] = {"7XQVB", "9RITB"}
    cohorts["framework_extremes"] = {"9G13A", "8V9XA"}
    proper_by_id = {t["instance"]: t for t in corrected}
    rows, score_rows, atom_rows = [], [], []
    for record in records:
        instance = record["target"]["instance"]
        if instance not in frozen:
            continue
        proper = proper_by_id[instance]
        if proper["hca"] == frozen[instance]["hca"]:
            continue
        row = reference_displacement(frozen[instance], proper)
        old_atoms = {(a["index"], a["atom"]): a for a in record["evidence"]["selected_alternate_atoms"]}
        new_atoms = {(a["index"], a["atom"]): a for a in record["evidence"]["occupancy_selected_alternate_atoms"]}
        for a in compare(proper, frozen[instance])["coordinate_mismatches"]:
            key = (a["residue_index"], a["atom"])
            before, after = old_atoms.get(key), new_atoms.get(key)
            atom_rows.append({"tid": row["tid"], "residue_index": key[0], "atom": key[1],
                              "frozen_altloc": before["altloc"] if before else "",
                              "frozen_occupancy": before["occupancy"] if before else None,
                              "occupancy_selected_altloc": after["altloc"] if after else "",
                              "selected_occupancy": after["occupancy"] if after else None,
                              "frozen_coord": json.dumps(a["frozen"]),
                              "occupancy_selected_coord": json.dumps(a["reconstructed"]),
                              "diagnosis": "alternate atom selection" if before and after else "requires inspection"})
        row.update({"in_"+name: row["tid"] in ids for name, ids in cohorts.items()})
        row.update({"saved_score_"+name: row["tid"] in scores for name, scores in models.items()})
        rows.append(row)
        for model, scores in models.items():
            if row["tid"] in scores:
                score_rows.append({"tid": row["tid"], "model": model,
                                   "status": "potentially reference-policy affected; NOT rescored",
                                   "saved_chain_cdr3": scores[row["tid"]]["chain_cdr3"],
                                   "saved_fw_cdr3": scores[row["tid"]]["fw_cdr3"],
                                   "saved_chain_framework": scores[row["tid"]]["chain_framework"]})
    csvout(OUT / (prefix+"reference_policy_displacements.csv"), rows)
    csvout(OUT / (prefix+"changed_atom_diagnosis.csv"), atom_rows)
    csvout(OUT / (prefix+"potentially_affected_model_scores.csv"), score_rows)
    changed = {r["tid"] for r in rows}
    impact = {"changed_frozen_targets": len(changed), "changed_tids": sorted(changed),
              "changed_atom_slots": len(atom_rows),
              "changed_atoms_with_lower_occupancy_in_frozen": sum(r["frozen_occupancy"] < r["selected_occupancy"]
                  for r in atom_rows if r["frozen_occupancy"] is not None and r["selected_occupancy"] is not None),
              "changed_atoms_requiring_additional_diagnosis": sum(r["diagnosis"] == "requires inspection" for r in atom_rows),
              "saved_score_records_potentially_affected": len(score_rows),
              "by_model": {name: len(changed & set(scores)) for name, scores in models.items()},
              "by_cohort": {name: {"changed": len(changed & ids), "total": len(ids),
                                    "changed_tids": sorted(changed & ids)} for name, ids in cohorts.items()},
              "model_predictions_rescored": 0,
              "scope": "Reference-to-reference displacement only. Saved-score/cohort intersections identify potential impact; new model RMSDs, ranking changes and threshold crossings are not measured."}
    impact["directly_changed_regions"] = {region: sum(r[f"changed_atom_slots_{region}"] > 0 for r in rows)
                                         for region in ("cdr1", "cdr2", "cdr3", "framework")}
    impact["reference_rmsd_maxima_angstrom"] = {key: max((r[key] for r in rows), default=0.)
        for key in ("chain_full_nca_c_cb_rmsd", "fw_full_nca_c_cb_rmsd", "chain_full_ca_rmsd", "fw_full_ca_rmsd",
                    "chain_cdr1_nca_c_cb_rmsd", "fw_cdr1_nca_c_cb_rmsd", "chain_cdr2_nca_c_cb_rmsd",
                    "fw_cdr2_nca_c_cb_rmsd", "chain_cdr3_nca_c_cb_rmsd", "fw_cdr3_nca_c_cb_rmsd")}
    dump(OUT / (prefix+"reference_policy_impact.json"), impact)
    return impact


def write_report(summary, excluded):
    s, impact = summary, summary["reference_policy_impact"]
    exclusions = "\n".join(f"| {e['tid']} | {', '.join(e['current_conditions'])} | unknown |" for e in excluded)
    cohorts = "\n".join(f"| {name} | {v['changed']} / {v['total']} |"
                         for name, v in impact["by_cohort"].items())
    model_counts = ", ".join(f"{k}: {v}" for k, v in impact["by_model"].items())
    maxima = "\n".join(f"| {name} | {v:.9f} |" for name, v in impact["reference_rmsd_maxima_angstrom"].items())
    report = f"""# Reference reconstruction and alternate-atom policy

## Result and scope

This is newly reconstructed code, not recovered contemporaneous extraction code.
All {s['reconstructed_records']} source representatives were rebuilt. The saved 597-target
replay was reloaded and equals the frozen targets as JSON values, including every
metadata field, sequence, span, coordinate and null. The frozen benchmark
references and saved model scores are unchanged.

Source: [SAbDab2 Zenodo 20083995]({ARCHIVE_URL}); {ARCHIVE_BYTES:,} bytes.
Full archive MD5 `{s['archive_hashes']['md5']}` and SHA256
`{s['archive_hashes']['sha256']}` were computed. The MD5 matches Zenodo.
All three retained metadata members equal the 16 October cached bytes.

Selection: 3,414 single-domain rows → 3,350 SD-H → 2,404 deposited strictly after
2020-05-01 → 604 CDR-H3 cluster representatives. Each representative is the
lowest numerical resolution, retaining original CSV order on ties. There is no
additional SAbDab train/test filter. Frozen order is preserved in the replay.

## Exact frozen validation

- Sequences and complete metadata: {s['frozen_exact_matches']} / 597 exact.
- CDR spans: {s['cdr_spans_checked']} / 1,791, derived from IMGT 27–38, 56–65,
  105–117 including insertion codes; all CDR strings agree.
- Residues: {s['residues_compared']:,}; numeric coordinate components:
  {s['coordinate_components_compared']:,}; unequal components: {s['unequal_coordinate_components']}.
- Null atom slots: {s['null_atom_slots_compared']:,}; mask differences: {s['null_mask_mismatches']}.
- Maximum absolute coordinate component difference: {s['max_abs_coordinate_difference_angstrom']} Å.
- Raw author chain IDs and author residue numbers/insertion codes identify IMGT
  residues. The first deposited model, ATOM records and polymer HETATM records
  are read. Noncanonical residue names are preserved as X for strict identity
  checking; MSE is not silently converted to methionine.
- Coordinates are rounded to three decimal places and compared exactly, without
  superposition or tolerance. The exact comparison tests the stored precision;
  it does not claim sub-millångström agreement with an unavailable original program.

## New reference-policy finding

Frozen replay requires the last source atom row to win when alternate locations
duplicate an atom. This policy is inferred from raw-versus-frozen agreement;
it is not a recovered statement about the historical program. Highest occupancy
per atom, with first source row on ties and shared blank-altloc atoms retained,
changes {impact['changed_frozen_targets']} / 597 references ({impact['changed_atom_slots']}
atom slots, {s['occupancy_aware_unequal_coordinate_components']} scalar components).
Of these, {impact['changed_atoms_with_lower_occupancy_in_frozen']} frozen atom selections have
lower occupancy. Maximum scalar change is {s['occupancy_aware_max_abs_difference_angstrom']:.9f} Å.
Rounding is applied identically to both policies; the detailed atom diagnosis
distinguishes alternate selection from any remaining unexplained difference.
Remaining atom differences requiring diagnosis: {impact['changed_atoms_requiring_additional_diagnosis']}.

The occupancy-aware 604-record file is a separate comparison artifact; it is not
a replacement benchmark. This finding prevents describing the frozen references
as having been built with highest-occupancy altloc selection. Model performance
after changing reference policy remains unmeasured.

Reference-to-reference RMSDs use matched N/CA/C/CB atoms with proper-rotation
Kabsch fits on either the complete chain or framework (excluding all three CDRs).
Regional and full-chain Cα RMSDs use those same transforms, not a separate Cα fit.
Every changed target has full and regional RMSDs, per-atom RMSDs and maximum
displacements in `reference_policy_displacements.csv`. Values below are maxima
across changed targets, in Å:

| Metric | Maximum |
|---|---:|
{maxima}

## Saved scores and sensitivity cases potentially affected

There are {impact['saved_score_records_potentially_affected']} saved target/model score records
on changed references ({model_counts}). These are reference-policy exposure
counts, not measured score differences. No prediction was rescored. Original
benchmark-wide prediction CIFs are not available in the inspected structure
output directories; a small visualization PDB collection does not establish
complete prediction provenance. No new model ranking, threshold crossing,
bootstrap interval or robustness result is claimed.

| Existing cohort | Changed references / cohort size |
|---|---:|
{cohorts}

Exact target lists are in `reference_policy_impact.json`; original saved scores
on affected targets are in `potentially_affected_model_scores.csv`. Directly
changed-region counts (which can overlap): {json.dumps(impact['directly_changed_regions'])}.

## Seven exclusions: current evidence versus history

Current rejection requires nonempty contiguous CDR spans, agreement with source
CDR strings, presence of every numbered residue, matching residue identities,
complete N/CA/C backbone coordinates, no glycine Cβ and at least one numbered
residue in each of IMGT framework segments 1–26, 39–55, 66–104 and 118–128.
The framework-segment guard is newly specified: elementary residue/span/atom
checks alone accepted 8XTYA, whose numbering starts at 30 and lacks framework 1.
This guard must not be represented as a recovered historical coverage threshold.
Missing non-glycine Cβ
alone is retained and represented by null, consistent with frozen data.
These are explicit reconstructed checks, not recovered historical thresholds.

| Representative | Current reproducible conditions | Historical per-ID reason |
|---|---|---|
{exclusions}

Full per-residue failures and numbering keys are in `exclusions.json`.
The retained aggregate record lists three non-contiguous CDR-H1 spans,
three coordinate/sequence identity mismatches and one coverage failure.
No verified historical per-ID mapping or original extraction program was found.
Empty and non-contiguous CDRs remain distinct current conditions; aggregate
historical category names are not assigned to individual IDs.

## Reproduction and files

From the benchmark root, with the retained archive, run:

```sh
uv run --no-project --python 3.12 --with-requirements environments/reference.lock python src/reconstruct_revision_references.py --archive outputs/revision_2026-10-17/reference_reconstruction/sources/splits.tar.gz
```

Omit `--offline` once if dependencies are absent. Python {s['python']}, Biopython
{s['biopython']} and NumPy {s['numpy']} were used. `--self-check` tests selection,
insertion, model/chain filtering, masks, spans and rigid fitting. `--download`
fetches this exact archive using eight validated HTTP ranges; full MD5 is required.
`--partial` outputs are explicitly provisional and must not be used as completion evidence.

- `targets_reconstructed_597.json`: exact frozen replay, original target order.
- `representatives_604.json`: all reconstructed targets with source keys, issues,
  selected alternates, source row/ties and coordinate SHA256 values.
- `occupancy_aware_604.json`: separate occupancy-aware coordinates.
- `frozen_validation.json`, `occupancy_aware_frozen_comparison.json`: per-target comparisons.
- `changed_atom_diagnosis.csv`: per-atom alternate labels, occupancies and both coordinates.
- `reference_policy_displacements.csv`, `reference_policy_impact.json`: impact evidence.
- `exclusions.json`: reconstructed failures; historical per-ID reasons remain null.
- `source_manifest.json`, `output_manifest.json`, `sources/download_receipt.json`:
  provenance and checksums; `run.log`, `progress.log`, `failures.log` (if present): run history.

This reconstruction makes no claim to recover its historical predecessor and
recalculates no model scores. The upstream release and extracted member hashes
define the reconstruction sources.
"""
    (OUT / "REPORT.md").write_text(report)


def self_check():
    header = "data_test\nloop_\n" + "\n".join("_atom_site."+k for k in (
        "group_PDB", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_comp_id",
        "label_atom_id", "label_alt_id", "occupancy", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")) + "\n"
    lines = ["ATOM H 111 A ALA CA A .4 1 2 3 1", "ATOM H 111 A ALA CA B .6 4 5 6 1",
             "ATOM H 111 A ALA CA C .6 7 8 9 1", "ATOM H 111 A ALA N . 1 0 0 0 1",
             "ATOM H 111 A ALA CA D 1 9 9 9 2", "ATOM X 111 A ALA CA . 1 8 8 8 1"]
    r, a, model = parse_atoms((header+"\n".join(lines)+"\n").encode(), "H")
    assert r[(111, "A")]["atoms"]["CA"]["coord"] == [4., 5., 6.]
    assert r[(111, "A")]["atoms"]["N"]["altloc"] == ""
    assert len(a) == 3 and model == "1"
    replay, _, _ = parse_atoms((header+"\n".join(lines)+"\n").encode(), "H", False)
    assert replay[(111, "A")]["atoms"]["CA"]["coord"] == [7., 8., 9.]
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    moved = points @ np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])+7.
    assert np.allclose(fit_coordinates(moved, points, np.ones(4, bool)), points, atol=1e-12)
    # A null CB must stay null through both policies; every numbered insertion is a separate residue.
    numbering = [["G", [n, " "]] for n in (27, 56, 105)]
    fake = {"INSTANCE": "pdb_0000test_H_+", "PDB_ID": "pdb_0000test", "Hchain": "H",
            "SABDAB_ID": "test", "cdrh3_cluster": "test", "VH_numbering_list": json.dumps(numbering),
            "VH_numerable_seq": "GGG", "CDRH1": "G", "CDRH2": "G", "CDRH3": "G",
            **{k: "test" for k in ("resolution", "method", "holo", "ab_split", "PDBdepo")}}
    gly = [f"ATOM H {n} ? GLY {atom} . 1 0 0 0 1" for n in (27, 56, 105) for atom in ATOMS[:3]]
    target, evidence = reconstruct(fake, (header+"\n".join(gly)+"\n").encode())
    assert [x["condition"] for x in evidence["issues"]] == ["missing_framework_segment"]
    assert all(q[3] is None for q in target["hca"])
    assert target["cdr1"] == [0, 1] and target["cdr3"] == [2, 3]
    changed = json.loads(json.dumps(target))
    changed["hca"][0][3] = [0., 0., 0.]
    assert compare(changed, target)["null_mask_mismatches"] == 1


def download_archive():
    """Range requests retain transport evidence and verify the complete archive."""
    source = OUT / "sources"
    source.mkdir(parents=True, exist_ok=True)
    chunks = [(i, start, min(start+110_000_000, ARCHIVE_BYTES)-1)
              for i, start in enumerate(range(0, ARCHIVE_BYTES, 110_000_000))]

    def fetch(chunk):
        i, start, end = chunk
        path = source / f"archive_range_{i:02d}.bin"
        req = Request(ARCHIVE_URL, headers={"Range": f"bytes={start}-{end}"})
        with urlopen(req, timeout=120) as response:
            assert response.status == 206
            assert response.headers["Content-Range"] == f"bytes {start}-{end}/{ARCHIVE_BYTES}"
            with path.open("wb") as f:
                shutil.copyfileobj(response, f, length=1024*1024)
            assert path.stat().st_size == end-start+1
            receipt = {"range": [start, end], "status": response.status,
                       "content_range": response.headers["Content-Range"], **checksums(path)}
        dump(source / f"archive_range_{i:02d}.json", receipt)
        print(f"Downloaded range {i}: {path.stat().st_size} bytes", flush=True)
        return path

    with ThreadPoolExecutor(max_workers=8) as pool:
        paths = list(pool.map(fetch, chunks))
    archive = source / "splits.tar.gz"
    with archive.open("wb") as out:
        for path in paths:
            with path.open("rb") as f:
                shutil.copyfileobj(f, out, length=1024*1024)
    hashes = checksums(archive)
    assert archive.stat().st_size == ARCHIVE_BYTES and hashes["md5"] == ARCHIVE_MD5
    dump(source / "download_receipt.json", {"url": ARCHIVE_URL, "bytes": ARCHIVE_BYTES, **hashes})
    print("Full archive checksum verified", hashes, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--partial", action="store_true", help="Provisional extraction from an incomplete download")
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--download", action="store_true", help="Fetch and checksum this exact source release")
    args = parser.parse_args()
    self_check()
    if args.download:
        download_archive()
        return
    if args.self_check:
        print("Self-check passed: occupancy, ties, replay, blank altloc, insertion, chain/model, null mask, sequence/CDR spans and rigid fit")
        return
    if not args.partial and args.archive is None:
        args.archive = OUT / "sources/splits.tar.gz"
    OUT.mkdir(parents=True, exist_ok=True)
    logfile = (OUT / ("progress.log" if args.partial else "run.log")).open("a", buffering=1)

    def log(message):
        line = datetime.now(timezone.utc).isoformat() + " " + message
        print(line, flush=True)
        print(line, file=logfile, flush=True)

    log("START " + " ".join(sys.argv))
    csv.field_size_limit(10_000_000)
    source = list(csv.DictReader((CACHE / "ab_split_sd.csv").open()))
    heavy = [r for r in source if r["type"] == "SD-H"]
    eligible = [r for r in heavy if r["PDBdepo"] > "2020-05-01"]
    groups = {}
    for r in eligible:
        groups.setdefault(r["cdrh3_cluster"], []).append(r)
    representatives = [min(g, key=lambda r: float(r["resolution"])) for g in groups.values()]
    assert (len(source), len(heavy), len(eligible), len(representatives)) == (3414, 3350, 2404, 604)
    frozen_list = json.loads(FROZEN.read_text())
    frozen_hash_at_start = checksums(FROZEN)["sha256"]
    frozen = {t["instance"]: t for t in frozen_list}
    assert len(frozen) == len(frozen_list) == 597
    wanted = {r["INSTANCE"] for r in representatives}
    assert set(frozen) <= wanted
    coords_dir = OUT / "sources/coordinates"
    coords_dir.mkdir(parents=True, exist_ok=True)
    archive_hashes = None
    members = {}
    if args.archive:
        if not args.partial:
            archive_hashes = checksums(args.archive)
            assert archive_hashes["md5"] == ARCHIVE_MD5, archive_hashes
            log("Archive full MD5 verified " + ARCHIVE_MD5)
        try:
            with tarfile.open(args.archive, "r|gz") as archive:
                for member in archive:
                    leaf = Path(member.name).name
                    if member.isfile() and ((leaf.endswith(".cif") and Path(leaf).stem in wanted)
                                             or leaf in ("README.md", "ab_split_sd.csv", "ab_split.csv")):
                        payload = archive.extractfile(member).read()
                        assert len(payload) == member.size
                        path = coords_dir / leaf if leaf.endswith(".cif") else OUT / "sources" / leaf
                        path.write_bytes(payload)
                        members[leaf] = {"archive_member": member.name, "bytes": len(payload)}
                        if leaf.endswith(".csv") or leaf == "README.md":
                            assert payload == (gzip.open(CACHE / "ab_split.csv.gz", "rb").read() if leaf == "ab_split.csv" else (CACHE / leaf).read_bytes()), f"cached source differs: {leaf}"
        except (EOFError, tarfile.ReadError) as exc:
            if not args.partial:
                raise
            log("Provisional archive boundary: " + str(exc))
    records, validations, excluded, occupancy_validations, corrected = [], [], [], [], []
    conditions = Counter()
    for r in representatives:
        instance = r["INSTANCE"]
        path = coords_dir / (instance + ".cif")
        if not path.exists():
            continue
        payload = path.read_bytes()
        target, evidence = reconstruct(r, payload, occupancy_aware=False)
        proper, proper_evidence = reconstruct(r, payload, occupancy_aware=True)
        corrected.append(proper)
        evidence["altloc_policy"] = "last source atom row, exact frozen replay; NOT recommended occupancy selection"
        evidence["occupancy_selected_alternate_atoms"] = proper_evidence["selected_alternate_atoms"]
        evidence.update({"instance": instance, "tid": target["pdb"]+target["H"],
                         "retained_in_frozen": instance in frozen,
                         "historical_individual_rejection_reason": None,
                         "current_conditions": sorted({x["condition"] for x in evidence["issues"]}),
                         "current_primary_condition": evidence["issues"][0]["condition"] if evidence["issues"] else None,
                         "source_coordinate_sha256": hashlib.sha256(payload).hexdigest(),
                         "source_csv_row_1based": source.index(r)+2,
                         "minimum_resolution_tie_count": sum(float(x["resolution"]) == float(r["resolution"])
                                                             for x in groups[r["cdrh3_cluster"]])})
        conditions.update(evidence["current_conditions"])
        if instance in frozen:
            evidence["frozen_validation"] = compare(target, frozen[instance])
            validations.append({"instance": instance, **evidence["frozen_validation"]})
            occupancy_validations.append({"instance": instance, **compare(proper, frozen[instance])})
        else:
            excluded.append(evidence)
        records.append({"target": target, "evidence": evidence})
        if len(records) % 50 == 0:
            log(f"Reconstructed {len(records)}/604; frozen checked {len(validations)}/597")
    prefix = "partial_" if args.partial else ""
    dump(OUT / (prefix + "representatives_604.json"), records)
    dump(OUT / (prefix + "exclusions.json"), excluded)
    by_id = {r["target"]["instance"]: r["target"] for r in records}
    dump(OUT / (prefix + "targets_reconstructed_597.json"), [by_id[t["instance"]] for t in frozen_list if t["instance"] in by_id])
    dump(OUT / (prefix + "frozen_validation.json"), validations)
    dump(OUT / (prefix + "occupancy_aware_604.json"), corrected)
    dump(OUT / (prefix + "occupancy_aware_frozen_comparison.json"), occupancy_validations)
    impact = impact_analysis(records, corrected, frozen, prefix)
    summary = {"label": "RECONSTRUCTION, not recovered historical extraction code",
               "validation_status": "PROVISIONAL" if args.partial else "PENDING",
               "provisional": args.partial, "source_url": ARCHIVE_URL,
               "archive_hashes": archive_hashes, "archive_expected_md5": ARCHIVE_MD5,
               "archive_md5_verified": archive_hashes is not None,
               "selection_counts": {"single_domain": len(source), "sd_h": len(heavy),
                                    "post_deposition_filter": len(eligible), "representatives": len(representatives)},
               "representative_rule": "lowest numerical resolution; original CSV order on ties",
               "reconstructed_records": len(records), "frozen_records_checked": len(validations),
               "frozen_exact_matches": sum(v["exact_match"] for v in validations),
               "sequences_exact": sum(v["sequence_match"] for v in validations),
               "span_sets_exact": sum(v["cdr_spans_match"] for v in validations),
               "cdr_spans_checked": len(validations)*3,
               "coordinate_components_compared": sum(v["coordinate_components_compared"] for v in validations),
               "residues_compared": sum(len(frozen[v["instance"]]["hseq"]) for v in validations),
               "null_atom_slots_compared": sum(a is None for v in validations for q in frozen[v["instance"]]["hca"] for a in q),
               "unequal_coordinate_components": sum(v["unequal_coordinate_components"] for v in validations),
               "null_mask_mismatches": sum(v["null_mask_mismatches"] for v in validations),
               "max_abs_coordinate_difference_angstrom": max((v["max_abs_coordinate_difference_angstrom"] for v in validations), default=0),
               "current_condition_counts_records": dict(conditions),
               "current_primary_condition_counts": dict(Counter(r["evidence"]["current_primary_condition"]
                   for r in records if r["evidence"]["current_primary_condition"] is not None)),
               "new_domain_coverage_guard": "at least one numbered residue in each framework segment; explicitly reconstructed, not a recovered historical threshold",
               "current_accepted": sum(not r["evidence"]["issues"] for r in records),
               "current_rejected": sum(bool(r["evidence"]["issues"]) for r in records),
               "excluded_records_checked": len(excluded),
               "frozen_replay_altloc_policy": "last source atom row wins; inferred from raw-versus-frozen comparison",
               "occupancy_aware_altloc_policy": "highest occupancy per atom; first source row on ties; blank shared atoms retained",
               "occupancy_aware_frozen_exact_matches": sum(v["exact_match"] for v in occupancy_validations),
               "occupancy_aware_unequal_coordinate_components": sum(v["unequal_coordinate_components"] for v in occupancy_validations),
               "occupancy_aware_max_abs_difference_angstrom": max((v["max_abs_coordinate_difference_angstrom"] for v in occupancy_validations), default=0),
               "historical_aggregate_only": {"noncontiguous_cdr1": 3, "residue_identity_mismatch": 3, "coverage_failure": 1},
               "historical_per_id_mapping_available": False,
               "reference_policy_impact": impact,
               "python": platform.python_version(), "biopython": Bio.__version__, "numpy": np.__version__}
    dump(OUT / (prefix + "summary.json"), summary)
    if not args.partial:
        assert len(records) == 604 and len(validations) == 597 and len(excluded) == 7
        assert all(v["exact_match"] for v in validations), "Frozen differences; see frozen_validation.json"
        assert all(not r["evidence"]["issues"] for r in records if r["evidence"]["retained_in_frozen"])
        assert all(r["issues"] for r in excluded), "An old exclusion is accepted by reconstructed conditions"
        assert json.loads((OUT / "targets_reconstructed_597.json").read_text()) == frozen_list
        assert len(json.loads((OUT / "representatives_604.json").read_text())) == 604
        assert len(json.loads((OUT / "occupancy_aware_604.json").read_text())) == 604
        assert checksums(FROZEN)["sha256"] == frozen_hash_at_start
        summary["validation_status"] = "PASS_EXACT_FROZEN_REPLAY; OCCUPANCY_POLICY_DIFFERENCES_DISCLOSED"
        summary["saved_outputs_reloaded_and_validated"] = True
        dump(OUT / "summary.json", summary)
        log("PASS saved outputs reloaded: all 597 frozen references exactly match; all 604 representatives reconstructed")
    else:
        log(f"PROVISIONAL {summary['frozen_exact_matches']}/{len(validations)} frozen references exactly match")
    evidence_inputs = [FROZEN, CACHE / "README.md", CACHE / "ab_split_sd.csv", CACHE / "ab_split.csv.gz", Path(__file__),
                       ROOT / "src/revision_source_metadata.py",
                       ROOT / "outputs/structure/simplefold/metrics/per_target_scores.json",
                       ROOT / "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
                       ROOT / "outputs/structure/nbb2/per_target_scores_nbb2.json"]
    evidence_inputs.extend(CACHE.parent / name for name in ("common_lower_risk_subset.csv", "uncertainty_groups_36.csv", "coverage_groups_40.csv"))
    dump(OUT / (prefix + "source_manifest.json"), {"inputs": {
        str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, **checksums(p)}
        for p in evidence_inputs},
        "archive": {"url": ARCHIVE_URL, "hashes": archive_hashes, "expected_md5": ARCHIVE_MD5,
                    "expected_bytes": ARCHIVE_BYTES}, "extracted_members": members})
    if not args.partial:
        write_report(summary, excluded)
        dump(OUT / "output_manifest.json", {p.name: {"bytes": p.stat().st_size, **checksums(p)}
             for p in sorted(OUT.iterdir()) if p.suffix in (".json", ".csv", ".md")
             and not p.name.startswith("partial_") and p.name != "output_manifest.json"})
    log("END " + json.dumps(summary))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        if OUT.exists():
            with (OUT / "failures.log").open("a") as f:
                f.write(datetime.now(timezone.utc).isoformat() + "\n" + traceback.format_exc() + "\n")
        raise
