"""Frozen-embedding linear probe: does ESM2-3B linearly encode nanobody CDR-H3
kinked/extended conformation? Decides controller-vs-side-stream for SimpleFold.

Reads outputs/esm-probe/{esm_probe_embeddings.npz,esm_probe_meta.json} (from the
potatoguy44/esm-probe Kaggle kernel) and data/inputs/targets_vhh.json (for framework
sequences, used only to build a leakage-safe split).
"""
import json

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

RNG = np.random.default_rng(0)

# ---- load ---------------------------------------------------------------
meta = json.load(open("outputs/esm-probe/esm_probe_meta.json"))
npz = np.load("outputs/esm-probe/esm_probe_embeddings.npz")
cdr3_mean, anchor = npz["cdr3_mean"].astype(np.float32), npz["anchor"].astype(np.float32)
y = np.array([1 if r["kink"] == "kinked" else 0 for r in meta])
tids = [r["tid"] for r in meta]
print(f"n={len(y)}  kinked={y.sum()}  extended={(1-y).sum()}")

# ---- framework-cluster grouping (germline proxy), for leakage-safe splits ----
# ponytail: k-mer Jaccard clustering on the framework sequence, not a real IGHV
# germline call (would need IgBLAST/ANARCI-germline). Cheap, dependency-free,
# and it does the one job the split needs: same-germline-like nanobodies must not
# straddle a train/test boundary, since attention lets framework context leak
# into CDR3-pooled embeddings.
targets = {f"{t['pdb']}{t['H']}": t for t in json.load(open("data/inputs/targets_vhh.json"))}


def framework_seq(t):
    seq = t["hseq"]
    mask = [True] * len(seq)
    for lo, hi in (t["cdr1"], t["cdr2"], t["cdr3"]):
        for i in range(lo, hi):
            mask[i] = False
    return "".join(c for c, m in zip(seq, mask) if m)


K = 5
kmers = {}
for tid in tids:
    fw = framework_seq(targets[tid])
    kmers[tid] = {fw[i : i + K] for i in range(len(fw) - K + 1)}

parent = {tid: tid for tid in tids}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


THRESH = 0.5
n = len(tids)
for i in range(n):
    si = kmers[tids[i]]
    for j in range(i + 1, n):
        sj = kmers[tids[j]]
        inter = len(si & sj)
        if inter == 0:
            continue
        jac = inter / len(si | sj)
        if jac >= THRESH:
            union(tids[i], tids[j])

cluster_of = {tid: find(tid) for tid in tids}
from collections import Counter

sizes = Counter(cluster_of.values())
print(f"framework clusters (Jaccard>={THRESH}): {len(sizes)}, "
      f"largest={max(sizes.values())}, singletons={sum(1 for v in sizes.values() if v==1)}")
groups = np.array([list(sizes.keys()).index(cluster_of[tid]) for tid in tids])

N_SPLITS = 5
assert len(sizes) >= N_SPLITS, f"only {len(sizes)} clusters, need >= {N_SPLITS}"

# ---- probe fitting --------------------------------------------------------
gkf = GroupKFold(n_splits=N_SPLITS)


def oof_auc(X, y, groups, model_fn):
    proba = cross_val_predict(model_fn(), X, y, groups=groups, cv=gkf, method="predict_proba")
    return roc_auc_score(y, proba[:, 1])


def shuffled_auc(X, y, groups, model_fn, n_shuffles):
    aucs = []
    for _ in range(n_shuffles):
        ys = RNG.permutation(y)
        proba = cross_val_predict(model_fn(), X, ys, groups=groups, cv=gkf, method="predict_proba")
        aucs.append(roc_auc_score(ys, proba[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def lr():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))


def mlp():
    return make_pipeline(
        StandardScaler(),
        MLPClassifier(hidden_layer_sizes=(64,), alpha=1e-2, max_iter=2000, random_state=0),
    )


NLAYER = cdr3_mean.shape[1]  # 37
results = {"per_layer_linear_cdr3mean": []}
for layer in range(NLAYER):
    auc = oof_auc(cdr3_mean[:, layer, :], y, groups, lr)
    results["per_layer_linear_cdr3mean"].append({"layer": layer, "auc": round(float(auc), 4)})

best = max(results["per_layer_linear_cdr3mean"], key=lambda r: r["auc"])
print("\nper-layer linear probe (CDR3-mean), best:", best)
for r in results["per_layer_linear_cdr3mean"]:
    print(r)

# mean-across-layers variant (uniform combine, proxy for "typical" ESM readout)
mean_layers = cdr3_mean.mean(axis=1)
auc_mean_layers = oof_auc(mean_layers, y, groups, lr)
print(f"\nmean-across-layers linear AUC: {auc_mean_layers:.4f}")

# anchor-residue variant (IMGT 105 + 117), best layer only, both linear and mlp
anchor_best = anchor[:, best["layer"], :, :].reshape(len(y), -1)
auc_anchor_lin = oof_auc(anchor_best, y, groups, lr)
print(f"anchor-residue linear AUC (layer {best['layer']}): {auc_anchor_lin:.4f}")

# majority-class baseline (the only well-defined baseline under GroupKFold --
# a per-cluster "germline lookup" is undefined here because GroupKFold holds
# each cluster out entirely, so a held-out cluster has zero training exposure
# by construction)
maj_auc = 0.5  # constant predictor, AUC is undefined/chance by definition
maj_acc = max(y.mean(), 1 - y.mean())
print(f"majority-class baseline: accuracy={maj_acc:.4f}  AUC=0.5 (by definition)")

# ---- control task (Hewitt & Liang): shuffle labels, refit, same splits ----
print("\ncontrol task (shuffled labels), best linear variant, 50 shuffles:")
ctrl_mean, ctrl_std = shuffled_auc(cdr3_mean[:, best["layer"], :], y, groups, lr, 50)
print(f"  control AUC: {ctrl_mean:.4f} +/- {ctrl_std:.4f}")
gap = best["auc"] - ctrl_mean
print(f"  selectivity gap (real - control): {gap:.4f}")

# ---- nonlinear probe (only meaningful to inspect if linear underperforms) ----
print("\nnonlinear (MLP) probe, best linear layer, 5-fold GroupKFold:")
auc_mlp = oof_auc(cdr3_mean[:, best["layer"], :], y, groups, mlp)
print(f"  MLP AUC: {auc_mlp:.4f}")
print("control task, MLP, 20 shuffles:")
ctrl_mlp_mean, ctrl_mlp_std = shuffled_auc(cdr3_mean[:, best["layer"], :], y, groups, mlp, 20)
print(f"  MLP control AUC: {ctrl_mlp_mean:.4f} +/- {ctrl_mlp_std:.4f}")
gap_mlp = auc_mlp - ctrl_mlp_mean

# ---- pre-registered decision ----------------------------------------------
print("\n=== DECISION (pre-registered thresholds: AUC>=0.70 AND gap>=0.15) ===")
linear_pass = best["auc"] >= 0.70 and gap >= 0.15
mlp_pass = auc_mlp >= 0.70 and gap_mlp >= 0.15
print(f"linear:    AUC={best['auc']:.4f} (layer {best['layer']}), gap={gap:.4f} -> {'PASS' if linear_pass else 'FAIL'}")
print(f"nonlinear: AUC={auc_mlp:.4f}, gap={gap_mlp:.4f} -> {'PASS' if mlp_pass else 'FAIL'}")
if linear_pass:
    print("TIER 1: signal linearly present -> controller")
elif mlp_pass:
    print("TIER 2: signal present, nonlinear only -> controller (trunk is nonlinear downstream)")
else:
    print("TIER 3: no signal in frozen ESM2 embedding -> side-stream needed")

out = {
    "n": int(len(y)), "n_kinked": int(y.sum()), "n_extended": int((1 - y).sum()),
    "n_clusters": len(sizes), "largest_cluster": max(sizes.values()),
    "per_layer_linear_cdr3mean": results["per_layer_linear_cdr3mean"],
    "best_layer": best["layer"], "best_layer_auc": round(float(best["auc"]), 4),
    "mean_across_layers_auc": round(float(auc_mean_layers), 4),
    "anchor_best_layer_auc": round(float(auc_anchor_lin), 4),
    "majority_baseline_acc": round(float(maj_acc), 4),
    "control_mean": round(ctrl_mean, 4), "control_std": round(ctrl_std, 4),
    "selectivity_gap": round(float(gap), 4),
    "mlp_auc": round(float(auc_mlp), 4),
    "mlp_control_mean": round(ctrl_mlp_mean, 4), "mlp_control_std": round(ctrl_mlp_std, 4),
    "mlp_gap": round(float(gap_mlp), 4),
    "linear_pass": bool(linear_pass), "mlp_pass": bool(mlp_pass),
}
json.dump(out, open("outputs/esm-probe/esm_probe_results.json", "w"), indent=2)
print("\nwrote outputs/esm-probe/esm_probe_results.json")
