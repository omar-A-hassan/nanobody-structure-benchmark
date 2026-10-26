"""Screen available source complexes for antigen and heavy-domain interface overlap.

This is an entry-linked source audit, not proof that a checkpoint saw antigens.
Frozen benchmark scores and the published VHH/H3 screen are never overwritten.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import time

from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-26/complex_overlap"
SOURCES = ROOT / "outputs/revision_2026-10-16/sources"
POLYMERS = {"PROTEIN", "PEPTIDE", "DNA", "RNA"}
MODELS = ("SimpleFold_pre_cutoff_proxy", "IgFold", "NanoBodyBuilder2")
THRESHOLD = 0.9
CUTOFF = 5.0


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def table(path, rows):
    if rows:
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def aligner():
    return PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-5, extend_gap_score=-0.5)


def identity(a, b, engine):
    if a == b:
        return 1.0
    counts = engine.align(a, b)[0].counts()
    return counts.identities / (counts.identities + counts.mismatches + counts.gaps)


def residue_map(query, source, engine):
    """Map source sequence indices onto homologous query indices, including substitutions."""
    if query == source:
        return dict(enumerate(range(len(query))))
    alignment = engine.align(query, source)[0]
    return {j: i for q, s in zip(*alignment.aligned)
            for i, j in zip(range(*q), range(*s))}


def load_sources():
    with gzip.open(SOURCES / "ab_split.csv.gz", "rt", newline="") as f:
        rows = list(csv.DictReader(f))
    nbb = {r["PDB"].upper() for r in csv.DictReader((SOURCES / "nbb2_nanobody_data.csv").open())
           if r["set"] == "train/val"}
    igfold = set(json.loads((SOURCES / "igfold_experimental_sequence_pool.json").read_text())["pdb_ids"])
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    target_ids = {t["instance"]: t["pdb"] + t["H"] for t in targets}
    entries = []
    coverage = []
    for row in rows:
        pdb = row["PDB_ID"][-4:].upper()
        models = []
        if row["PDBdepo"] <= "2020-05-01":
            models.append(MODELS[0])
        if pdb in igfold:
            models.append(MODELS[1])
        if pdb in nbb:
            models.append(MODELS[2])
        tid = target_ids.get(row["INSTANCE"], "")
        if not models and not tid:
            continue
        if not row["VH_numerable_seq"]:
            continue
        parts = [row[k].split("/") for k in ("agchains", "agtypes", "agresolvedseqs", "ag_index_list")]
        assert len({len(p) for p in parts}) == 1, row["INSTANCE"]
        polymer_count = 0
        for chain, kind, sequence, indices in zip(*parts):
            if kind not in POLYMERS or not sequence:
                continue
            numbering = json.loads(indices)
            assert sequence == "".join(a for a, _ in numbering), (row["INSTANCE"], chain)
            assert len(sequence) >= 5
            polymer_count += 1
            entries.append({"key": row["INSTANCE"] + ":" + chain, "instance": row["INSTANCE"],
                            "pdb": pdb, "tid": tid, "chain": chain, "kind": kind,
                            "sequence": sequence, "antigen_numbering": numbering,
                            "hchain": row["Hchain"], "vhh": row["VH_numerable_seq"],
                            "heavy_numbering": json.loads(row["VH_numbering_list"]), "models": models})
        if tid:
            coverage.append({"tid": tid, "source_holo": row["holo"],
                             "retained_antigen_types": row["agtypes"],
                             "polymer_antigen_chains": polymer_count})
    return entries, coverage


def screen():
    OUT.mkdir(parents=True, exist_ok=True)
    entries, coverage = load_sources()
    queries = [e for e in entries if e["tid"]]
    training = [e for e in entries if e["models"]]
    groups = defaultdict(list)
    for e in training:
        groups[(e["kind"], e["sequence"])].append(e["key"])
    query_groups = defaultdict(list)
    for e in queries:
        query_groups[(e["kind"], e["sequence"])].append(e["key"])
    engine = aligner()
    matches = []
    started = time.monotonic()
    for i, ((kind, query), keys) in enumerate(query_groups.items(), 1):
        for (other_kind, other), source_keys in groups.items():
            # Protein/peptide labels describe the same alphabet; nucleic acids stay separate.
            if (kind in {"PROTEIN", "PEPTIDE"}) != (other_kind in {"PROTEIN", "PEPTIDE"}):
                continue
            if kind not in {"PROTEIN", "PEPTIDE"} and kind != other_kind:
                continue
            if min(len(query), len(other)) / max(len(query), len(other)) < THRESHOLD:
                continue
            score = identity(query, other, engine)
            if score >= THRESHOLD:
                matches.extend({"query_key": q, "source_key": s, "antigen_identity": score}
                               for q in keys for s in source_keys)
        if i % 20 == 0 or i == len(query_groups):
            print(f"Antigen sequence screen {i}/{len(query_groups)}; matches={len(matches)}; seconds={time.monotonic()-started:.1f}", flush=True)
    needed = {m[k] for m in matches for k in ("query_key", "source_key")}
    needed.update(e["key"] for e in queries)
    save(OUT / "entries.json", [e for e in entries if e["key"] in needed])
    save(OUT / "antigen_matches.json", matches)
    table(OUT / "benchmark_complex_coverage.csv", coverage)
    table(OUT / "antigen_matches.csv", matches)
    save(OUT / "screen_inputs.json", {"threshold": THRESHOLD, "distance_cutoff_A": CUTOFF,
         "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
         [SOURCES / "ab_split.csv.gz", SOURCES / "nbb2_nanobody_data.csv",
          SOURCES / "igfold_experimental_sequence_pool.json", ROOT / "data/inputs/targets_vhh.json"]},
         "query_antigen_chains": len(queries), "query_unique_typed_sequences": len(query_groups),
         "source_antigen_chains": len(training), "source_unique_typed_sequences": len(groups),
         "matches": len(matches), "elapsed_seconds": time.monotonic()-started})


def normalize(value):
    return "" if value.strip() in {".", "?", ""} else value.strip()


def coordinate_path(instance):
    # Source chain IDs are case-sensitive; the default macOS filesystem is not.
    return OUT / "coordinates" / (hashlib.sha256(instance.encode()).hexdigest()[:20] + ".cif")


def extract_contacts(payload, entry):
    if entry["hchain"] == entry["chain"]:
        raise ValueError(f"Antigen and heavy chain share an identifier: {entry['key']}")
    d = MMCIF2Dict(io.StringIO(payload.decode()))
    c = {k: d["_atom_site." + k] for k in
         ("auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_atom_id", "type_symbol",
          "label_comp_id", "occupancy", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")}
    actual_chains = set(c["auth_asym_id"])
    chain_map = {chain: chain for chain in (entry["hchain"], entry["chain"])}
    if not set(chain_map) <= actual_chains:
        raise ValueError(f"Missing exact case-sensitive source chain: {entry['key']}")
    reverse_chain_map = {v: k for k, v in chain_map.items()}
    assert len(reverse_chain_map) == 2
    maps = {}
    for chain, numbering in [(entry["hchain"], entry["heavy_numbering"]),
                             (entry["chain"], entry["antigen_numbering"])]:
        maps[chain] = {(int(n), normalize(ins)): j for j, (_, (n, ins)) in enumerate(numbering)}
    heavy_expected = {j: aa for j, (aa, _) in enumerate(entry["heavy_numbering"])}
    atoms = {}
    model = c["pdbx_PDB_model_num"][0]
    for i, actual_chain in enumerate(c["auth_asym_id"]):
        chain = reverse_chain_map.get(actual_chain)
        if chain not in maps or c["pdbx_PDB_model_num"][i] != model or c["type_symbol"][i] in {"H", "D"}:
            continue
        residue = maps[chain].get((int(c["auth_seq_id"][i]), normalize(c["pdbx_PDB_ins_code"][i])))
        if residue is None:
            continue
        if chain == entry["hchain"]:
            assert protein_letters_3to1.get(c["label_comp_id"][i]) == heavy_expected[residue], f"Heavy residue mismatch: {entry['key']} position {residue}"
        key = chain, residue, c["label_atom_id"][i]
        occupancy = float(c["occupancy"][i])
        if not np.isfinite(occupancy) or occupancy <= 0:
            continue
        if key not in atoms or occupancy > atoms[key][0]:
            atoms[key] = (occupancy, [float(c[k][i]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")])
    arrays = {}
    for chain in maps:
        selected = [(key[1], value[1]) for key, value in atoms.items() if key[0] == chain]
        if not selected:
            raise ValueError(f"No positive-occupancy matching heavy atoms: {entry['key']} chain {chain}")
        arrays[chain] = ([r for r, _ in selected], np.array([v for _, v in selected]))
    h_ids, h_xyz = arrays[entry["hchain"]]
    a_ids, a_xyz = arrays[entry["chain"]]
    neighbours = cKDTree(h_xyz).query_ball_tree(cKDTree(a_xyz), CUTOFF)
    contacts = {(h_ids[i], a_ids[j]) for i, js in enumerate(neighbours) for j in js}
    return {"contacts": sorted(contacts), "epitope": sorted({j for _, j in contacts}),
            "paratope": sorted({i for i, _ in contacts}), "chain_mapping": chain_map,
            "coordinate_sha256": hashlib.sha256(payload).hexdigest()}


def overlap(query, source, mapping):
    """Unaligned source positions remain distinct and count against similarity."""
    mapped = {mapping[x] if x in mapping else ("unmapped", x) for x in source}
    common = len(set(query) & mapped)
    union = len(set(query) | mapped)
    return {"intersection": common, "query_size": len(query), "source_size": len(source),
            "jaccard": common / union if union else None,
            "query_recall": common / len(query) if query else None}


def contact_similarity(query, source, hmap, amap):
    mapped = {(hmap[i], amap[j]) if i in hmap and j in amap else ("unmapped", i, j)
              for i, j in source}
    shared = len(set(map(tuple, query)) & mapped)
    union = len(set(map(tuple, query)) | mapped)
    return {"intersection": shared, "query_size": len(query), "source_size": len(source),
            "jaccard": shared / union if union else None,
            "query_recall": shared / len(query) if query else None}


def contacts(archive):
    entries = {e["key"]: e for e in json.loads((OUT / "entries.json").read_text())}
    matches = json.loads((OUT / "antigen_matches.json").read_text())
    needed = {m[k] for m in matches for k in ("query_key", "source_key")}
    needed.update(k for k, e in entries.items() if e["tid"])
    entries = {k: e for k, e in entries.items() if k in needed}
    save(OUT / "entries.json", list(entries.values()))
    by_instance = defaultdict(list)
    for key in sorted(needed):
        by_instance[entries[key]["instance"]].append(key)
    results = {}
    errors = []

    def process(instance, payload):
        for key in by_instance.pop(instance):
            try:
                results[key] = extract_contacts(payload, entries[key])
            except (ValueError, KeyError, AssertionError) as ex:
                errors.append({"key": key, "error": str(ex)})
        if len(results) % 50 < 3:
            print(f"Contact extraction {len(results)} interfaces; remaining files={len(by_instance)}", flush=True)

    for instance in list(by_instance):
        path = coordinate_path(instance)
        if path.exists():
            process(instance, path.read_bytes())
    if by_instance:
        assert archive and archive.is_file(), "Supply the local verified SAbDab2 archive with --archive"
        assert archive.stat().st_size == 876381859, "Unexpected source archive size"
        assert sha(archive) == "54e3b9cdae5f5ff57a4210866f87e86686203556eea6a033200ff55c0f81e2b7"
        with tarfile.open(archive, "r|gz") as tar:
            for member in tar:
                instance = Path(member.name).stem
                if member.isfile() and member.name.endswith(".cif") and instance in by_instance:
                    payload = tar.extractfile(member).read()
                    # Retain only the source files actually used by this new analysis.
                    path = coordinate_path(instance)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(payload)
                    process(instance, payload)
                if not by_instance:
                    break
    assert not by_instance, f"Source archive missing {len(by_instance)} required files"
    save(OUT / "native_contacts.json", results)
    save(OUT / "contact_extraction_errors.json", errors)
    engine = aligner()
    comparisons = []
    for i, m in enumerate(matches):
        q, s = entries[m["query_key"]], entries[m["source_key"]]
        if q["key"] not in results or s["key"] not in results:
            continue
        qc, sc = results[q["key"]], results[s["key"]]
        if not qc["contacts"] and not sc["contacts"]:
            continue
        amap = residue_map(q["sequence"], s["sequence"], engine)
        hmap = residue_map(q["vhh"], s["vhh"], engine)
        epitope = overlap(qc["epitope"], sc["epitope"], amap)
        interface = contact_similarity(qc["contacts"], sc["contacts"], hmap, amap)
        for model in s["models"]:
            comparisons.append({"tid": q["tid"], "model": model, **m,
                                "same_pdb": q["pdb"] == s["pdb"],
                                "epitope_jaccard": epitope["jaccard"],
                                "epitope_query_recall": epitope["query_recall"],
                                "interface_jaccard": interface["jaccard"],
                                "interface_query_recall": interface["query_recall"],
                                "query_contact_pairs": interface["query_size"],
                                "source_contact_pairs": interface["source_size"]})
        if (i+1) % 500 == 0:
            print(f"Mapped interface comparisons {i+1}/{len(matches)}", flush=True)
    table(OUT / "interface_comparisons.csv", comparisons)
    summarize(entries, matches, results, errors, comparisons)


def summarize(entries, matches, results, errors, comparisons):
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    common = {r["tid"] for r in json.loads((ROOT / "outputs/structure/nbb2/per_target_scores_nbb2.json").read_text())}
    screened = set(json.loads((ROOT / "outputs/revision_2026-10-17/overlap_extension/overlap_extension_summary.json").read_text())["cohorts"]["extended_lower_risk"]["tids"])
    flags = []
    for t in targets:
        tid = t["pdb"] + t["H"]
        qkeys = [k for k, e in entries.items() if e["tid"] == tid]
        for model in MODELS:
            selected = [m for m in matches if m["query_key"] in qkeys and model in entries[m["source_key"]]["models"]]
            measured = [c for c in comparisons if c["tid"] == tid and c["model"] == model]
            distinct = [c for c in measured if not c["same_pdb"]]
            flags.append({"tid": tid, "model": model, "common_588": tid in common,
                          "screened_393": tid in screened, "polymer_chains": len(qkeys),
                          "measured_interfaces": sum(k in results for k in qkeys),
                          "interfaces_with_contacts": sum(k in results and bool(results[k]["contacts"]) for k in qkeys),
                          "antigen_ge90": bool(selected),
                          "exact_antigen": any(m["antigen_identity"] == 1 for m in selected),
                          "distinct_pdb_antigen_ge90": any(entries[m["source_key"]]["pdb"] != t["pdb"] for m in selected),
                          "max_epitope_jaccard": max((c["epitope_jaccard"] for c in measured), default=None),
                          "max_interface_jaccard": max((c["interface_jaccard"] for c in measured), default=None),
                          "distinct_pdb_max_epitope_jaccard": max((c["epitope_jaccard"] for c in distinct), default=None),
                          "distinct_pdb_max_interface_jaccard": max((c["interface_jaccard"] for c in distinct), default=None)})
    table(OUT / "per_target_model_complex_overlap.csv", flags)
    counts = {}
    for name, ids in [("full_597", {t["pdb"]+t["H"] for t in targets}), ("common_588", common), ("screened_393", screened)]:
        counts[name] = {}
        for model in MODELS:
            rows = [r for r in flags if r["tid"] in ids and r["model"] == model]
            counts[name][model] = {"targets": len(rows), "polymer_antigen_evaluable": sum(r["polymer_chains"]>0 for r in rows),
                "native_heavy_interface_evaluable": sum(r["measured_interfaces"]>0 for r in rows),
                **{k: sum(r[k] for r in rows) for k in ["antigen_ge90", "exact_antigen", "distinct_pdb_antigen_ge90"]},
                "epitope_jaccard_ge50": sum(r["max_epitope_jaccard"] is not None and r["max_epitope_jaccard"]>=.5 for r in rows),
                "interface_jaccard_ge50": sum(r["max_interface_jaccard"] is not None and r["max_interface_jaccard"]>=.5 for r in rows)}
    all_entries, coverage = load_sources()
    source_coverage = {}
    for model in MODELS:
        pool = [e for e in all_entries if model in e["models"]]
        source_coverage[model] = {"polymer_antigen_chain_records": len(pool),
            "pdb_entries": len({e["pdb"] for e in pool}),
            "unique_typed_sequences": len({(e["kind"], e["sequence"]) for e in pool})}
    bad_instances = {e["key"].split(":")[0] for e in errors if "positive-occupancy" in e["error"]}
    excluded = sorted(t["pdb"]+t["H"] for t in targets if t["instance"] in bad_instances)
    quality_ids = common - set(excluded)
    quality = {"excluded_targets": excluded, "common_targets_retained": len(quality_ids), "medians": {}}
    for model, path in [("SimpleFold", "outputs/structure/simplefold/metrics/per_target_scores.json"),
                        ("IgFold", "outputs/structure/igfold/metrics/per_target_scores_igfold.json"),
                        ("NanoBodyBuilder2", "outputs/structure/nbb2/per_target_scores_nbb2.json")]:
        rows = [r for r in json.loads((ROOT/path).read_text()) if r["tid"] in quality_ids]
        quality["medians"][model] = {key: float(np.median([r[key] for r in rows]))
                                     for key in ("chain_cdr3", "fw_cdr3")}
    save(OUT / "summary.json", {"status": "COMPLETED_AVAILABLE_SOURCE_COMPLEX_AUDIT", "counts": counts,
        "source_pool_coverage": source_coverage,
        "coordinate_quality_sensitivity": quality,
        "source_interface_records": len(results), "extraction_errors": errors, "interface_comparisons": len(comparisons),
        "methods": {"antigen_sequence_identity": "global affine-gap optimal-score alignment; identical positions/all columns including gaps; match2/mismatch-1/open-5/extend-0.5",
                    "antigen_threshold": THRESHOLD, "native_contact_cutoff_A": CUTOFF,
                    "epitope": "antigen residues within cutoff of heavy-variable-domain heavy atoms",
                    "interface": "heavy-domain/antigen residue contact pairs; sequence-align both partners; unmapped residues retained in union",
                    "altloc": "positive occupancy only; highest atom occupancy; first source row on ties; first deposited model",
                    "source_file_identity": "exact case-sensitive archive member; case-sensitive instance hashed for portable cache filenames"},
        "limits": ["Entry-linked experimental-source complexes; not verified checkpoint training-chain membership or antigen input exposure.",
                   "SimpleFold source pool is a pre-cutoff antibody-complex proxy, not its full experimental/distilled training set.",
                   "Only source-retained polymer antigens can be sequence-screened; nonpolymeric or removed partners are not classified as nonoverlapping.",
                   "Only antigen pairs passing global90 identity receive epitope/contact-pattern comparisons; lower-identity structural analogues are not screened.",
                   "Heavy-domain contact-pattern similarity is not a DockQ score or a rigid-body docking prediction.",
                   "Protein/peptide types share an alphabet; DNA/RNA are screened separately."]})
    print(json.dumps(counts, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["screen", "contacts"])
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    if args.phase == "screen":
        screen()
    else:
        contacts(args.archive)


if __name__ == "__main__":
    main()
