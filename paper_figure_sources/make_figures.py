"""Create the four evidence figures used in main.tex."""

from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PAPER = Path(__file__).parent
DATA = PAPER.parent / "kaggle-nanobody-benchmark" / "outputs"
OUT = PAPER / "figures"
OUT.mkdir(exist_ok=True)

COLORS = {"SimpleFold": "#0077BB", "IgFold": "#009988", "NanoBodyBuilder2": "#EE7733"}
MODELS = list(COLORS)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.dpi": 300})


def load(path):
    return json.loads(path.read_text())


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def rank(values):
    order = np.argsort(values, kind="mergesort")
    result = np.empty(len(values), dtype=float)
    result[order] = np.arange(len(values), dtype=float)
    return result


def spearman(x, y):
    return float(np.corrcoef(rank(np.asarray(x)), rank(np.asarray(y)))[0, 1])


def figure1_benchmark():
    paths = {
        "SimpleFold": DATA / "structure/simplefold/metrics/per_target_scores.json",
        "IgFold": DATA / "structure/igfold/metrics/per_target_scores_igfold.json",
        "NanoBodyBuilder2": DATA / "structure/nbb2/per_target_scores_nbb2.json",
    }
    rows = {name: {r["tid"]: r for r in load(path)} for name, path in paths.items()}
    tids = set.intersection(*(set(value) for value in rows.values()))
    medians = {
        "CDR-H3": [float(np.median([rows[name][tid]["chain_cdr3"] for tid in tids])) for name in MODELS],
        "Framework": [float(np.median([rows[name][tid]["chain_framework"] for tid in tids])) for name in MODELS],
    }
    tail = [100 * sum(rows[name][tid]["chain_cdr3"] > 6 for tid in tids) / len(tids)
            for name in MODELS]
    assert len(tids) == 588
    assert np.allclose(medians["CDR-H3"], [2.944, 2.8845, 2.4835], atol=0.001)
    assert np.allclose(medians["Framework"], [0.8205, 0.8265, 0.793], atol=0.001)
    assert np.allclose(tail, [14.5, 8.0, 7.8], atol=0.2)

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.25), gridspec_kw={"wspace": 0.46})
    x = np.arange(len(MODELS))
    for ax, region, ylabel in [(axes[0], "CDR-H3", "Median RMSD (Å)"),
                               (axes[1], "Framework", "Median RMSD (Å)")]:
        bars = ax.bar(x, medians[region], color=[COLORS[m] for m in MODELS], width=0.68)
        ax.set_xticks(x, ["SimpleFold", "IgFold", "NBB2"], rotation=25, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, max(medians[region]) * 1.20)
        ax.grid(axis="y", alpha=0.20, linewidth=0.6)
        ax.set_axisbelow(True)
        for bar, value in zip(bars, medians[region]):
            ax.text(bar.get_x() + bar.get_width() / 2, value + max(medians[region]) * 0.035, f"{value:.2f}",
                    ha="center", va="bottom", fontsize=7)
    bars = axes[2].bar(x, tail, color=[COLORS[m] for m in MODELS], width=0.68)
    axes[2].set_xticks(x, ["SimpleFold", "IgFold", "NBB2"], rotation=25, ha="right")
    axes[2].set_ylabel("Targets > 6 Å (%)")
    axes[2].set_ylim(0, 17.5)
    axes[2].grid(axis="y", alpha=0.20, linewidth=0.6)
    axes[2].set_axisbelow(True)
    for bar, value in zip(bars, tail):
        axes[2].text(bar.get_x() + bar.get_width() / 2, value + 0.45, f"{value:.1f}",
                     ha="center", va="bottom", fontsize=7)
    for label, ax in zip(("(a)", "(b)", "(c)"), axes):
        ax.text(0.01, 0.99, label, transform=ax.transAxes, fontweight="bold", va="top")
    fig.subplots_adjust(bottom=0.28, top=0.98)
    save(fig, "figure1_benchmark")


def figure2_uncertainty():
    rows = load(DATA / "trajectory/full-597/traj_results.json")["rows"]
    spread = [r["spread_cdr3"] for r in rows]
    error = [r["err"] for r in rows]
    rho_spread = spearman(spread, error)
    assert abs(rho_spread - 0.413) < 0.01
    labels = ["Q1", "Half", "End", "Spread", "No align."]
    values = [spearman([r["adiv_q1"] for r in rows], error),
              spearman([r["adiv_half"] for r in rows], error),
              spearman([r["adiv_final"] for r in rows], error),
              rho_spread,
              spearman([r["div_final"] for r in rows], error)]
    assert np.allclose(values, [0.443, 0.411, 0.397, 0.413, -0.025], atol=0.01)

    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.25), gridspec_kw={"wspace": 0.45},
                             constrained_layout=True)
    axes[0].scatter(spread, error, s=9, alpha=0.42, color="#0077BB", edgecolors="none")
    axes[0].set_xlabel("CDR-H3 sample spread (Å)")
    axes[0].set_ylabel("CDR-H3 RMSD (Å)")
    axes[0].set_ylim(0, max(error) * 1.06)
    axes[0].text(0.96, 0.94, f"Spearman ρ = {rho_spread:.3f}", transform=axes[0].transAxes,
                 ha="right", va="top", bbox={"facecolor": "white", "edgecolor": "0.8", "pad": 2})
    axes[0].grid(alpha=0.20, linewidth=0.6)
    axes[0].set_axisbelow(True)

    colors = ["#0077BB", "#0077BB", "#0077BB", "#009988", "#666666"]
    axes[1].bar(np.arange(len(values)), values, color=colors, width=0.68)
    axes[1].axhline(0, color="black", linewidth=0.7)
    axes[1].set_xticks(np.arange(len(values)), labels, rotation=22, ha="right")
    axes[1].tick_params(axis="x", labelsize=7.5)
    axes[1].set_ylabel("Spearman ρ with CDR-H3 RMSD")
    axes[1].set_ylim(-0.08, 0.52)
    axes[1].grid(axis="y", alpha=0.20, linewidth=0.6)
    axes[1].set_axisbelow(True)
    axes[0].text(0.01, 0.99, "(a)", transform=axes[0].transAxes, fontweight="bold", va="top")
    axes[1].text(0.01, 0.99, "(b)", transform=axes[1].transAxes, fontweight="bold", va="top")
    save(fig, "figure2_uncertainty")


def figure3_prior_coverage():
    hit = load(DATA / "hitrate/run/hitrate_results.json")
    taus = sorted(hit["arms"], key=float)
    hard_hit, easy_hit, hard_best = [], [], []
    for tau in taus:
        rows = hit["arms"][tau]["rows"]
        hard = [r for r in rows if r["group"] == "hard"]
        easy = [r for r in rows if r["group"] == "easy"]
        hard_hit.append(100 * sum(r["fw_best"] < 2 for r in hard) / len(hard))
        easy_hit.append(100 * sum(r["fw_best"] < 2 for r in easy) / len(easy))
        hard_best.append(float(np.median([r["fw_best"] for r in hard])))
    assert np.allclose(hard_hit, [0.0, 3.333, 3.333, 6.667], atol=0.01)
    assert np.allclose(easy_hit, [100.0, 100.0, 100.0, 90.0], atol=0.01)
    assert np.allclose(hard_best, [7.599, 6.6435, 4.951, 4.3645], atol=0.001)

    x = np.arange(len(taus))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.25), gridspec_kw={"wspace": 0.45},
                             constrained_layout=True)
    axes[0].bar(x - width / 2, easy_hit, width, label="Easy controls", color="#33BBEE")
    axes[0].bar(x + width / 2, hard_hit, width, label="Difficult targets", color="#EE7733")
    axes[0].set_xticks(x, taus)
    axes[0].set_xlabel("Noise setting τ")
    axes[0].set_ylabel("Hit rate below 2 Å (%)")
    axes[0].set_ylim(0, 120)
    axes[0].legend(frameon=False, loc="upper left", ncol=2, handlelength=1.2, columnspacing=0.8)
    axes[0].grid(axis="y", alpha=0.20, linewidth=0.6)
    axes[0].set_axisbelow(True)
    axes[1].plot(x, hard_best, marker="o", color="#EE7733", linewidth=1.8)
    axes[1].axhline(2, color="0.35", linestyle="--", linewidth=0.9, label="2 Å target")
    axes[1].set_xticks(x, taus)
    axes[1].set_xlabel("Noise setting τ")
    axes[1].set_ylabel("Median best CDR-H3 RMSD (Å)")
    axes[1].set_ylim(0, 8.4)
    axes[1].legend(frameon=False, loc="upper right")
    axes[1].grid(alpha=0.20, linewidth=0.6)
    axes[1].set_axisbelow(True)
    axes[0].text(0.01, 0.99, "(a)", transform=axes[0].transAxes, fontweight="bold", va="top")
    axes[1].text(0.16, 0.99, "(b)", transform=axes[1].transAxes, fontweight="bold", va="top")
    save(fig, "figure3_prior_coverage")


def figure4_esm_probe():
    result = load(DATA / "esm-probe/esm_probe_results.json")
    values = np.asarray([row["auc"] for row in result["per_layer_linear_cdr3mean"]], dtype=float)
    best = int(result["best_layer"])
    assert len(values) == 37
    assert abs(values[best] - result["best_layer_auc"]) < 1e-8
    assert abs(result["control_mean"] - 0.5) < 0.01

    fig, ax = plt.subplots(figsize=(3.35, 2.25))
    layers = np.arange(len(values))
    ax.plot(layers, values, color="#0077BB", marker="o", markersize=2.8, linewidth=1.5,
            label="Linear probe")
    ax.axhline(0.5, color="#666666", linestyle="--", linewidth=0.9, label="Shuffled-label control")
    ax.scatter([best], [values[best]], s=38, color="#EE7733", zorder=3)
    ax.set_xlabel("ESM2 layer")
    ax.set_ylabel("Five-fold ROC AUC")
    ax.set_ylim(0.45, 0.95)
    ax.set_xlim(-1, 37)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=2, handlelength=1.5, columnspacing=0.8)
    ax.grid(alpha=0.20, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.text(0.03, 0.98, "(a)", transform=ax.transAxes, fontweight="bold", va="top")
    ax.annotate(f"layer {best}", xy=(best, values[best]), xytext=(best - 7, values[best] - 0.045),
                arrowprops={"arrowstyle": "-", "color": "#555555", "linewidth": 0.7}, fontsize=7)
    fig.subplots_adjust(bottom=0.25, top=0.84)
    save(fig, "figure4_esm_probe")


if __name__ == "__main__":
    figure1_benchmark()
    figure2_uncertainty()
    figure3_prior_coverage()
    figure4_esm_probe()
    print(f"Wrote figures to {OUT}")
