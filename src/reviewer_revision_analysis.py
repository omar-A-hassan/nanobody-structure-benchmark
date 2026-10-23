"""Offline reviewer analyses. Run with uv; original scores and figures are not changed."""
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from Bio.Align import PairwiseAligner
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-16"
OUT.mkdir(parents=True, exist_ok=True)
SEED, BOOTS = 20261016, 4000
MODELS = ["SimpleFold", "IgFold", "NanoBodyBuilder2"]
PATHS = ["outputs/structure/simplefold/metrics/per_target_scores.json",
         "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
         "outputs/structure/nbb2/per_target_scores_nbb2.json"]
REGIONS = ["cdr1", "cdr2", "cdr3", "framework"]


def load(path):
    return json.loads((ROOT / path).read_text())


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def csvout(name, rows):
    assert rows
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def rounded(x):
    return str(Decimal(str(round(x, 10))).quantize(Decimal(".001"), rounding=ROUND_HALF_UP))


def ci(values):
    return np.quantile(values, [.025, .975]).tolist()


def rho(a, b):
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def partial(a, b, length):
    x, y, z = map(rankdata, (a, b, length))
    design = np.column_stack([np.ones(len(z)), z])
    x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(np.corrcoef(x, y)[0, 1])


def auc(pos, neg):
    ranks = rankdata(np.r_[pos, neg])
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2)
                 / (len(pos) * len(neg)))


aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-5, extend_gap_score=-.5)


def identity(a, b):
    """Global identity = identical columns / all columns, including terminal gaps."""
    counts = aligner.align(a, b)[0].counts()
    return counts.identities / (counts.identities + counts.mismatches + counts.gaps)


assert identity("AAAA", "AAAA") == 1
assert identity("AAAA", "AAAT") == .75
assert identity("AAAA", "AAA") == .75
assert rho([1, 1, 2], [2, 2, 3]) == 1
assert auc([1, 1], [1, 1]) == .5
assert rounded(3.2664999999999997) == "3.267"

targets = load("data/inputs/targets_vhh.json")
rows = {m: {r["tid"]: r for r in load(p)} for m, p in zip(MODELS, PATHS)}
by_instance = {r["instance"]: r["tid"] for r in rows["SimpleFold"].values()}
meta = {by_instance[t["instance"]]: t for t in targets}
allids = sorted(meta)
paired = sorted(set.intersection(*(set(r) for r in rows.values())))
assert len(allids) == 597 and len(paired) == 588
for t in targets:
    assert len(t["hseq"]) == len(t["hca"])
    for j in (1, 2, 3):
        lo, hi = t[f"cdr{j}"]
        assert t["hseq"][lo:hi] == t[f"CDRH{j}"]
    for aa, atoms in zip(t["hseq"], t["hca"]):
        assert all(a is not None for a in atoms[:3])
        if aa == "G":
            assert atoms[3] is None


def summary(ids):
    rng = np.random.default_rng(SEED)
    idx = rng.integers(len(ids), size=(BOOTS, len(ids)))
    result = {"n": len(ids), "models": {}, "paired_differences": {}}
    for model in MODELS:
        result["models"][model] = {}
        for mode in ("chain", "fw"):
            for region in REGIONS:
                key = f"{mode}_{region}"
                v = np.array([rows[model][tid][key] for tid in ids])
                result["models"][model][key] = {
                    "median": float(np.median(v)), "display": rounded(np.median(v)),
                    "median_ci95": ci(np.median(v[idx], axis=1)),
                    "mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                    "q25_q75": np.quantile(v, [.25, .75]).tolist(),
                    "above6_count": int((v > 6).sum()),
                    "above6_pct": float(100 * (v > 6).mean())}
        for other in MODELS[:MODELS.index(model)]:
            for mode in ("chain", "fw"):
                key = f"{mode}_cdr3"
                delta = np.array([rows[other][t][key] - rows[model][t][key] for t in ids])
                result["paired_differences"][f"{other} minus {model}: {key}"] = {
                    "median": float(np.median(delta)),
                    "ci95": ci(np.median(delta[idx], axis=1))}
    return result


results = {"bootstrap": {"replicates": BOOTS, "seed": SEED,
                         "method": "paired target bootstrap; percentile interval"},
           "full_paired": summary(paired),
           "excluding_partial_7XQVB": summary([t for t in paired if t != "7XQVB"]),
           "excluding_obsolete_9RITB": summary([t for t in paired if t != "9RITB"]),
           "excluding_partial_and_obsolete": summary([t for t in paired if t not in ("7XQVB", "9RITB")])}
results["missing_reference_cb"] = dict(Counter(aa for t in targets
    for aa, atoms in zip(t["hseq"], t["hca"]) if atoms[3] is None))
results["composition"] = {}
for label, ids in (("full", allids), ("paired", paired)):
    results["composition"][label] = {
        "n": len(ids), "holo_source_annotation": dict(Counter(meta[t]["holo"] for t in ids)),
        "method": dict(Counter(meta[t]["method"] for t in ids)),
        "unique_pdb": len({meta[t]["pdb"] for t in ids}),
        "unique_vhh": len({meta[t]["hseq"] for t in ids}),
        "unique_cdr3": len({meta[t]["CDRH3"] for t in ids}),
        "resolution_range": [min(float(meta[t]["resolution"]) for t in ids),
                             max(float(meta[t]["resolution"]) for t in ids)]}
groups = {}
for i, tid in enumerate(paired):
    groups.setdefault(meta[tid]["pdb"], []).append(i)
blocks = list(groups.values())
rng = np.random.default_rng(SEED)
draws = [np.concatenate([blocks[j] for j in rng.integers(len(blocks), size=len(blocks))])
         for _ in range(BOOTS)]
results["paired_pdb_block_bootstrap"] = {"n_pdbs": len(blocks), "effects": {}}
for name, d in results["full_paired"]["paired_differences"].items():
    pair, key = name.split(": ")
    a, b = pair.split(" minus ")
    v = np.array([rows[a][t][key] - rows[b][t][key] for t in paired])
    results["paired_pdb_block_bootstrap"]["effects"][name] = {
        "median": d["median"], "ci95": ci([np.median(v[ix]) for ix in draws])}
manifest = [{"tid": tid, **{k: meta[tid][k] for k in
             ("instance", "pdb", "H", "sabdab_id", "cdrh3_cluster", "hseq", "CDRH1",
              "CDRH2", "CDRH3", "resolution", "method", "holo", "ab_split", "PDBdepo")},
             "paired588": tid in paired} for tid in allids]
csvout("benchmark_manifest.csv", manifest)
(OUT / "benchmark_sequences.fasta").write_text("".join(
    f">{tid}\n{meta[tid]['hseq']}\n" for tid in allids))
csvout("region_scores_588.csv", [{"tid": t, "model": m, **{f"{mode}_{r}": rows[m][t][f"{mode}_{r}"]
       for mode in ("chain", "fw") for r in REGIONS}} for t in paired for m in MODELS])

print("Scoring and manifest complete; checking sequence redundancy", flush=True)
cache = OUT / "pairwise_identity.npz"
if cache.exists():
    cached = np.load(cache)
    assert cached["tids"].tolist() == allids
    vh, h3 = cached["vh"], cached["h3"]
else:
    vh, h3 = np.eye(len(allids)), np.eye(len(allids))
    for i, t in enumerate(allids):
        for j in range(i):
            u = allids[j]
            vh[i, j] = vh[j, i] = identity(meta[t]["hseq"], meta[u]["hseq"])
            h3[i, j] = h3[j, i] = identity(meta[t]["CDRH3"], meta[u]["CDRH3"])
    np.savez_compressed(cache, tids=np.array(allids), vh=vh, h3=h3)


def components(matrix, threshold, ids):
    # Connected components are conservative: no cross-component pair exceeds the threshold.
    parent = {i: i for i in ids}
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    positions = {t: i for i, t in enumerate(allids)}
    for i, t in enumerate(ids):
        for u in ids[:i]:
            if matrix[positions[t], positions[u]] >= threshold:
                parent[find(t)] = find(u)
    groups = {}
    for t in ids:
        groups.setdefault(find(t), []).append(t)
    groups = sorted(groups.values(), key=lambda x: min(x))
    representatives = [min(g, key=lambda t: (float(meta[t]["resolution"]), t)) for g in groups]
    return groups, representatives


results["redundancy"] = {}
clusterrows = []
for mode, matrix, threshold in (("vh90", vh, .9), ("vh95", vh, .95), ("cdr3_80", h3, .8)):
    groups, representatives = components(matrix, threshold, paired)
    results["redundancy"][mode] = {"n_groups": len(groups), "largest": max(map(len, groups)),
        "threshold": threshold, "representative_rule": "lowest numerical resolution, then target ID",
        "representatives": representatives, "scores": summary(representatives)}
    for g, rep in zip(groups, representatives):
        clusterrows.extend({"definition": mode, "cluster": rep, "tid": t,
                            "representative": t == rep} for t in g)
csvout("sequence_clusters.csv", clusterrows)
pos = {t: i for i, t in enumerate(allids)}
csvout("nearest_benchmark_neighbors.csv", [{"tid": t,
       "nearest_vhh": allids[int(np.argmax(np.where(np.arange(597) != pos[t], vh[pos[t]], -1)))],
       "max_vhh_identity": float(np.max(np.delete(vh[pos[t]], pos[t]))),
       "max_cdr3_identity": float(np.max(np.delete(h3[pos[t]], pos[t])))} for t in allids])

print("Redundancy complete; uncertainty controls", flush=True)
results["uncertainty"] = {}
for name, path in (("trajectory", "outputs/trajectory/full-597/traj_results.json"),
                   ("seed777_replication", "outputs/tnp/sf/tnp_sf.json")):
    data = load(path)["rows"]
    x = np.array([r["spread_cdr3"] for r in data])
    y = np.array([r["err"] for r in data])
    z = np.array([r["cdr3_len"] for r in data])
    rng = np.random.default_rng(SEED)
    estimates = []
    for _ in range(BOOTS):
        ix = rng.integers(len(data), size=len(data))
        estimates.append([rho(x[ix], y[ix]), partial(x[ix], y[ix], z[ix])])
    estimates = np.array(estimates)
    results["uncertainty"][name] = {"n": len(data), "spread_rho": rho(x, y),
        "spread_rho_ci95": ci(estimates[:, 0]), "length_rho": rho(z, y),
        "spread_length_rho": rho(x, z), "partial_rank_rho": partial(x, y, z),
        "partial_rank_ci95": ci(estimates[:, 1])}
hard = load("outputs/hard/run/hard_results.json")["rows"]
hp = np.array([[r["spread_cdr3"], r["cdr3_len"]] for r in hard if r["group"] == "hard"])
hn = np.array([[r["spread_cdr3"], r["cdr3_len"]] for r in hard if r["group"] == "easy"])
rng = np.random.default_rng(SEED)
estimates = []
for _ in range(BOOTS):
    p, n = hp[rng.integers(len(hp), size=len(hp))], hn[rng.integers(len(hn), size=len(hn))]
    a, b = auc(p[:, 0], n[:, 0]), auc(p[:, 1], n[:, 1])
    estimates.append([a, b, a - b])
estimates = np.array(estimates)
results["hard_easy_auc"] = {"hard_n": len(hp), "easy_n": len(hn),
    "spread_auc": auc(hp[:, 0], hn[:, 0]), "length_auc": auc(hp[:, 1], hn[:, 1]),
    "spread_ci95": ci(estimates[:, 0]), "length_ci95": ci(estimates[:, 1]),
    "auc_difference_ci95": ci(estimates[:, 2])}
csvout("uncertainty_groups_36.csv", [{**r, "selection_chain_cdr3": rows["SimpleFold"][r["tid"]]["chain_cdr3"]} for r in hard])
hit = load("outputs/hitrate/run/hitrate_results.json")
arm = next(iter(hit["arms"].values()))["rows"]
csvout("coverage_groups_40.csv", [{"tid": r["tid"], "group": r["group"],
       "selection_chain_cdr3": rows["SimpleFold"][r["tid"]]["chain_cdr3"],
       "selection_chain_framework": rows["SimpleFold"][r["tid"]]["chain_framework"]} for r in arm])
labels = load("data/inputs/tnp_crystals.json")["rows"]
csvout("probe_label_inclusion.csv", [{"tid": r["tid"], "label": r["kink"],
       "included": r["kink"] in ("kinked", "extended")} for r in labels])
results["probe_label_counts"] = dict(Counter(r["kink"] for r in labels))

nbsource = OUT / "sources/nbb2_nanobody_data.csv"
if nbsource.exists():
    with nbsource.open() as f:
        source = list(csv.DictReader(f))
    training = {r["PDB"].upper() for r in source if r["set"] == "train/val"}
    test = {r["PDB"].upper() for r in source if r["set"] == "test"}
    csvout("nbb2_pdb_overlap.csv", [{"tid": t, "pdb": meta[t]["pdb"],
           "in_author_train_val": meta[t]["pdb"].upper() in training,
           "in_author_test": meta[t]["pdb"].upper() in test} for t in allids])
    removed = [t for t in paired if meta[t]["pdb"].upper() in training]
    retained = [t for t in paired if t not in removed]
    results["nbb2_pdb_audit"] = {"source_train_val_pdbs": len(training),
        "source_test_pdbs": len(test), "overlap_all597": sum(meta[t]["pdb"].upper() in training for t in allids),
        "overlap_paired588": len(removed), "removed": removed, "retained_scores": summary(retained)}
    _, reps = components(vh, .9, retained)
    results["nbb2_pdb_audit"]["vh90_nonredundant_scores"] = summary(reps)
    results["nbb2_pdb_audit"]["vh90_representatives"] = reps

results["7R73A_benchmark_only_not_verified_docking_inputs"] = {m: rows[m]["7R73A"] for m in MODELS}
results["source_sha256"] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                           for p in PATHS + ["data/inputs/targets_vhh.json"]}
dump("revision_statistics.json", results)
csvout("benchmark_summary_588.csv", [{"model": m, "region_and_alignment": key, **v,
        "median_ci95": json.dumps(v["median_ci95"]), "q25_q75": json.dumps(v["q25_q75"])}
       for m, block in results["full_paired"]["models"].items() for key, v in block.items()])

print("Drawing benchmark distributions", flush=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
colors = ["#0077BB", "#009988", "#CC6677"]
for mode, title in (("chain", "Whole-chain alignment"), ("fw", "Framework alignment")):
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.8), layout="constrained")
    for ax, region, label in zip(axes.flat, REGIONS, ("CDR-H1", "CDR-H2", "CDR-H3", "Framework")):
        for model, color in zip(MODELS, colors):
            v = np.sort([rows[model][t][f"{mode}_{region}"] for t in paired])
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=color, label=model)
        ax.set(title=label, xlabel="RMSD (Å)", ylabel="Cumulative fraction", ylim=(0, 1.03))
        ax.grid(alpha=.2)
    ax = axes[1, 1]
    diffs = results["paired_pdb_block_bootstrap"]["effects"]
    names = [f"SimpleFold minus {m}: {mode}_cdr3" for m in MODELS[1:]]
    for j, name in enumerate(names):
        d = diffs[name]
        ax.plot(d["ci95"], [j, j], color=colors[j + 1], lw=2)
        ax.plot(d["median"], j, "o", color=colors[j + 1])
    ax.axvline(0, color=".5", lw=1, ls="--")
    ax.set(yticks=[0, 1], yticklabels=["SF − IgFold", "SF − NBB2"],
           xlabel="Paired median CDR-H3 difference (Å)", title="Effect estimates and 95% CIs", ylim=(-.6, 1.6))
    axes[1, 2].axis("off")
    handles, legends = axes[0, 0].get_legend_handles_labels()
    axes[1, 2].legend(handles, legends, loc="center", frameon=False)
    axes[1, 2].text(.5, .12, "588 targets scored by all three models\nN, Cα, C and Cβ; IMGT regions\nPaired PDB-block bootstrap; 4,000 draws", ha="center", transform=axes[1, 2].transAxes, fontsize=8)
    for letter, ax in zip("abcde", axes.flat):
        ax.text(-.12, 1.08, f"({letter})", transform=ax.transAxes, fontweight="bold")
    fig.suptitle(title, fontsize=12)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"benchmark_distributions_{mode}.{ext}", dpi=200, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
assert results["probe_label_counts"] == {"other": 79, "extended": 188, "kinked": 330}
assert results["full_paired"]["models"]["SimpleFold"]["chain_cdr3"]["display"] == "2.944"
print(json.dumps({"paired": len(paired), "redundancy_groups": {k: v["n_groups"] for k, v in results["redundancy"].items()},
                  "uncertainty": results["uncertainty"], "auc": results["hard_easy_auc"],
                  "nbb2_audit": {k: v for k, v in results.get("nbb2_pdb_audit", {}).items() if k in ("overlap_all597", "overlap_paired588")}}, indent=2))
