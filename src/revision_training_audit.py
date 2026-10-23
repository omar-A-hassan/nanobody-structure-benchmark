"""Conservative, incomplete training-overlap audit using public structural sources."""
from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import zipfile
import gzip

import numpy as np
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-16"
SRC = OUT / "sources"
csv.field_size_limit(10_000_000)
targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
with gzip.open(SRC / "ab_split.csv.gz", "rt", newline="") as f:
    source = list(csv.DictReader(f))
source_map = {r["INSTANCE"]: r for r in source}
assert len(source_map) == len(source) == 15641
for t in targets:
    r = source_map[t["instance"]]
    assert t["hseq"] == r["VH_numerable_seq"] and t["CDRH3"] == r["CDRH3"]
sd = [r for r in source if r["type"] == "SD-H"]
eligible = [r for r in sd if r["PDBdepo"] > "2020-05-01"]
assert len(sd) == 3350 and len(eligible) == 2404
assert len({r["cdrh3_cluster"] for r in eligible}) == 604
assert len({t["cdrh3_cluster"] for t in targets}) == 597
source_groups = {}
for r in eligible:
    source_groups.setdefault(r["cdrh3_cluster"], []).append(r)
for t in targets:
    assert t["instance"] == min(source_groups[t["cdrh3_cluster"]], key=lambda r: float(r["resolution"]))["INSTANCE"]
retained_instances = {t["instance"] for t in targets}
rejected = [min(g, key=lambda r: float(r["resolution"]))["INSTANCE"] for g in source_groups.values()
            if not {r["INSTANCE"] for r in g} & retained_instances]
assert len(rejected) == 7
links = [{"tid": t["pdb"] + t["H"], **{k: source_map[t["instance"]][k] for k in
          ("Hseq_expected", "VH_numbering_list", "agchains", "agtypes", "agresolvedseqs", "agexpectedseqs")}}
         for t in targets]
with (OUT / "source_complex_links.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(links[0])); w.writeheader(); w.writerows(links)

with (SRC / "nbb2_nanobody_data.csv").open() as f:
    nbb = list(csv.DictReader(f))
nbb_pdbs = {r["PDB"].upper() for r in nbb if r["set"] == "train/val"}
nbb_rows = [r for r in sd if r["PDB_ID"][-4:].upper() in nbb_pdbs]
nbb_sequences = sorted({(r["VH_numerable_seq"], r["CDRH3"]) for r in nbb_rows})
sf_proxy = [r for r in sd if r["PDBdepo"] <= "2020-05-01"]
sf_sequences = sorted({(r["VH_numerable_seq"], r["CDRH3"]) for r in sf_proxy})

pool = json.loads((SRC / "igfold_experimental_sequence_pool.json").read_text())
igseq = {tuple(row) for row in pool["heavy_h3_pairs"]}
igpdbs = set(pool["pdb_ids"])
assert len(igseq) == 2649 and len(igpdbs) == 4275

aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-5, extend_gap_score=-.5)
def ident(a, b):
    if not a or not b:
        return 0.
    c = aligner.align(a, b)[0].counts()
    return c.identities / (c.identities + c.mismatches + c.gaps)
assert ident("AAAA", "AAA") == .75

report = {"source_funnel_verified": {"all": len(source), "sd_h": len(sd), "post_deposition_filter": len(eligible),
           "source_cdr3_clusters": 604, "retained": 597, "rejected_representatives": rejected,
           "representative_rule": "lowest resolution; original source order for tied values; verified for every retained target"},
          "source_pdb_labels_not_verified_per_checkpoint": True,
          "nbb_source": {"train_val_pdbs": len(nbb_pdbs), "mapped_sd_h_pdbs": len({r["PDB_ID"] for r in nbb_rows}),
          "mapped_sequences": len(nbb_sequences), "unmapped_pdbs": sorted(nbb_pdbs - {r["PDB_ID"][-4:].upper() for r in nbb_rows})},
          "igfold_source": {"experimental_pdbs": len(igpdbs), "unique_heavy_h3_pairs": len(igseq)},
          "simplefold_proxy_source": {"n_sd_h_structures": len(sf_proxy), "n_sequence_pairs": len(sf_sequences),
             "description": "SAbDab2 SD-H entries deposited no later than 2020-05-01; conservative candidate pool, not confirmed SimpleFold training membership"},
          "limitations": ["IgFold distilled OAS structures not sequence-audited",
             "SimpleFold SwissProt/AFESM distilled structures not sequence-audited",
             "SimpleFold pre-cutoff candidate pool covers only SD-H entries in this SAbDab2 ML release, not its full experimental training set",
             "NBB2 sequences recovered only where public manifest PDBs map to SAbDab2 SD-H entries",
             "Antigen, epitope and interface overlap not evaluated",
             "Public manifests may differ from exact historical checkpoint membership"]}
records = []
for label, seqs, pdbs in (("nbb2", nbb_sequences, nbb_pdbs), ("igfold", sorted(igseq), igpdbs),
                        ("simplefold_pre_cutoff_proxy", sf_sequences, {r["PDB_ID"][-4:].upper() for r in sf_proxy})):
    cache = OUT / f"{label}_nearest_training.csv"
    if cache.exists():
        with cache.open() as f:
            block = list(csv.DictReader(f))
        assert {r["tid"] for r in block} == {t["pdb"] + t["H"] for t in targets}
        if label == "igfold" and any(r.get("cdr3_definition") != "Chothia93-102_IMGT-equivalent" for r in block):
            targetmap = {t["pdb"] + t["H"]: t for t in targets}
            for r in block:
                r["max_cdr3_identity"] = max(ident(targetmap[r["tid"]]["CDRH3"], b) for a, b in seqs)
                r["cdr3_definition"] = "Chothia93-102_IMGT-equivalent"
    else:
        block = []
        for i, t in enumerate(targets):
            vh_best = max(ident(t["hseq"], a) for a, b in seqs)
            h3_best = max(ident(t["CDRH3"], b) for a, b in seqs)
            block.append({"tid": t["pdb"] + t["H"], "pdb_overlap": t["pdb"].upper() in pdbs,
                          "max_vhh_identity": vh_best, "max_cdr3_identity": h3_best})
            if (i + 1) % 100 == 0:
                print(label, i + 1, "/597", flush=True)
    if label == "igfold":
        for r in block:
            r["cdr3_definition"] = "Chothia93-102_IMGT-equivalent"
    with cache.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(block[0])); w.writeheader(); w.writerows(block)
    report[label] = {"pdb_overlap_all": sum(str(r["pdb_overlap"]) == "True" for r in block),
                     "vhh_ge90_all": sum(float(r["max_vhh_identity"]) >= .9 for r in block),
                     "cdr3_ge80_all": sum(float(r["max_cdr3_identity"]) >= .8 for r in block)}
    records.append({r["tid"]: r for r in block})

paired_scores = [json.loads((ROOT / p).read_text()) for p in
    ("outputs/structure/simplefold/metrics/per_target_scores.json", "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
     "outputs/structure/nbb2/per_target_scores_nbb2.json")]
scores = [{r["tid"]: r for r in rows} for rows in paired_scores]
paired = sorted(set.intersection(*(set(s) for s in scores)))
meta = {t["pdb"] + t["H"]: t for t in targets}
releases_path = SRC / "rcsb_release_dates.json"
releases = json.loads(releases_path.read_text()) if releases_path.exists() else {}
filtered = []
subset_rows = []
for tid in paired:
    no_pdb = all(str(d[tid]["pdb_overlap"]) == "False" for d in records)
    low_vh = all(float(d[tid]["max_vhh_identity"]) < .9 for d in records)
    low_h3 = all(float(d[tid]["max_cdr3_identity"]) < .8 for d in records)
    release = releases.get(meta[tid]["pdb"].upper(), {}).get("initial_release_date", "")[:10]
    common_atoms = len({s[tid]["n_atoms_matched"] for s in scores}) == 1
    subset_rows.append({"tid": tid, "no_known_training_pdb": no_pdb, "vhh_below90": low_vh,
                       "cdr3_below80": low_h3, "initial_release_date": release,
                       "pdb_obsolete": meta[tid]["pdb"] == "9RIT", "same_matched_atom_count": common_atoms})
    if no_pdb and low_vh and low_h3 and meta[tid]["pdb"] != "9RIT" and common_atoms:
        filtered.append(tid)
with (OUT / "training_filter_flags.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(subset_rows[0])); w.writeheader(); w.writerows(subset_rows)

cache = np.load(OUT / "pairwise_identity.npz")
positions = {t: i for i, t in enumerate(cache["tids"].tolist())}
vhh_matrix = cache["vh"]
parent = {t: t for t in filtered}
def find(t):
    while parent[t] != t:
        parent[t] = parent[parent[t]]; t = parent[t]
    return t
for i, tid in enumerate(filtered):
    for other in filtered[:i]:
        if vhh_matrix[positions[tid], positions[other]] >= .9:
            parent[find(tid)] = find(other)
groups = {}
for t in filtered:
    groups.setdefault(find(t), []).append(t)
representatives = sorted(min(g, key=lambda t: (float(meta[t]["resolution"]), t)) for g in groups.values())
report["common_lower_risk_subset"] = {"n": len(representatives), "tids": representatives,
        "description": "No known IgFold experimental/NBB2 train-val PDB; <90% VHH and <80% CDR3 identity to audited sequences including potential pre-cutoff SimpleFold VHH candidates; 90% VHH connected-component representatives; obsolete 9RIT and unequal matched-atom counts excluded. NOT proven training-independent."}
for label, ids in (("all_known_pdb_excluded", [r["tid"] for r in subset_rows if r["no_known_training_pdb"]]),
                   ("lower_risk_nonredundant", representatives),
                   ("excluding_obsolete_and_partial", [t for t in paired if t not in ("9RITB", "7XQVB")])):
    rng = np.random.default_rng(20261016)
    ix = rng.integers(len(ids), size=(4000, len(ids)))
    result = {"n": len(ids), "medians": {}, "paired_effects": {}}
    pdb_blocks = {}
    for i, tid in enumerate(ids):
        pdb_blocks.setdefault(meta[tid]["pdb"], []).append(i)
    blocks = list(pdb_blocks.values())
    block_rng = np.random.default_rng(20261016)
    block_draws = [np.concatenate([blocks[j] for j in block_rng.integers(len(blocks), size=len(blocks))]) for _ in range(4000)]
    result["n_pdbs"] = len(blocks)
    for model, s in zip(("SimpleFold", "IgFold", "NanoBodyBuilder2"), scores):
        result["medians"][model] = {key: float(np.median([s[t][key] for t in ids]))
             for key in ("chain_cdr1", "chain_cdr2", "chain_cdr3", "chain_framework", "fw_cdr1", "fw_cdr2", "fw_cdr3", "fw_framework")}
    for j, model in ((1, "IgFold"), (2, "NanoBodyBuilder2")):
        for mode in ("chain", "fw"):
            v = np.array([scores[0][t][f"{mode}_cdr3"] - scores[j][t][f"{mode}_cdr3"] for t in ids])
            result["paired_effects"][f"SF minus {model}: {mode}"] = {"median": float(np.median(v)),
                "ci95": np.quantile(np.median(v[ix], axis=1), [.025, .975]).tolist(),
                "pdb_block_ci95": np.quantile([np.median(v[k]) for k in block_draws], [.025, .975]).tolist()}
    report[label] = result
with gzip.open(SRC / "ab_split.csv.gz", "rb") as stream:
    ab_hash = hashlib.file_digest(stream, "sha256").hexdigest()
report["sources_sha256"] = {"ab_split.csv": ab_hash, "nbb2_nanobody_data.csv": hashlib.sha256((SRC / "nbb2_nanobody_data.csv").read_bytes()).hexdigest(), "igfold_sabdab.zip": pool["source_archive"]["sha256"]}
with (OUT / "common_lower_risk_subset.csv").open("w", newline="") as f:
    fields = ["tid", "pdb", "public_release_date", "model", "chain_cdr1", "chain_cdr2", "chain_cdr3", "chain_framework", "fw_cdr1", "fw_cdr2", "fw_cdr3", "fw_framework"]
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for tid in representatives:
        for model, s in zip(("SimpleFold", "IgFold", "NanoBodyBuilder2"), scores):
            w.writerow({"tid": tid, "pdb": meta[tid]["pdb"], "public_release_date": releases[meta[tid]["pdb"]]["initial_release_date"][:10],
                        "model": model, **{k: s[tid][k] for k in fields[4:]}})
(OUT / "training_overlap_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({k: v for k, v in report.items() if k in ("nbb2", "igfold", "all_known_pdb_excluded", "lower_risk_nonredundant")}, indent=2))
