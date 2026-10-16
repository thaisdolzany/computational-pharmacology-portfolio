"""Figures for the solubility models, plus a train/test shift summary table.

Design rules: Okabe-Ito colours that keep 3:1 contrast on white, colour always
paired with a second cue (marker shape or direct label), every interval named
(95% bootstrap CI or SD) and n stated in each panel.
"""
import ast
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

BLUE, VERMILLION, GREEN = "#0072B2", "#D55E00", "#009E73"
SPLIT_STYLE = {"random": dict(color=BLUE, marker="o", label="Random split"),
               "scaffold": dict(color=VERMILLION, marker="s", label="Scaffold split")}
MODEL_LABEL = {"mean": "Training mean", "esol": "ESOL (refit)", "rf": "Random forest", "xgb": "XGBoost"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
    "axes.labelsize": 9, "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
    "savefig.dpi": 300, "savefig.facecolor": "white", "pdf.fonttype": 42,
})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}")
    plt.close(fig)


def fig_model_comparison(metrics, noise_rms):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), layout="constrained", sharey=True)
    order = ["mean", "esol", "rf", "xgb"]
    ypos = {m: i for i, m in enumerate(order)}
    for ax, metric, title in [(axes[0], "rmse", "a  Test RMSE (log units, lower is better)"),
                              (axes[1], "r2", "b  Test R²")]:
        for split, style in SPLIT_STYLE.items():
            sub = metrics[metrics.split == split].set_index("model").loc[order]
            offset = -0.13 if split == "random" else 0.13
            y = np.array([ypos[m] for m in order]) + offset
            lo = sub[metric] - sub[f"{metric}_ci"].map(lambda c: c[0])
            hi = sub[f"{metric}_ci"].map(lambda c: c[1]) - sub[metric]
            ax.errorbar(sub[metric], y, xerr=[lo, hi], fmt=style["marker"], color=style["color"],
                        ms=5, capsize=2.5, lw=1, label=style["label"])
        ax.set_title(title)
        ax.grid(axis="x", lw=0.4, color="#dddddd")
    axes[0].axvline(noise_rms, ls="--", lw=0.9, color="#555555")
    axes[0].text(noise_rms, 3.55, " between-source\n spread (RMS SD)", fontsize=7, color="#555555", va="top")
    axes[0].set_xlim(left=0)
    axes[1].axvline(0, lw=0.6, color="#999999")
    axes[0].set_yticks(range(len(order)), [MODEL_LABEL[m] for m in order])
    axes[0].set_ylim(-0.6, 3.6)
    axes[1].legend(loc="lower right")
    n = metrics.groupby("split").n_test.first()
    fig.supxlabel(f"Points: estimate; bars: 95% bootstrap CI (2,000 resamples). "
                  f"Test n = {n['random']:,} (random), {n['scaffold']:,} (scaffold).", fontsize=7.5)
    save(fig, "fig1_model_comparison")


def fig_pred_vs_obs(preds, metrics):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), layout="constrained", sharex=True, sharey=True)
    lims = (-13.5, 2.5)
    for ax, split, tag in [(axes[0], "random", "a"), (axes[1], "scaffold", "b")]:
        d = preds[(preds.split == split) & (preds.model == "xgb")]
        m = metrics[(metrics.split == split) & (metrics.model == "xgb")].iloc[0]
        ax.fill_between(lims, np.array(lims) - 1, np.array(lims) + 1, color="#eeeeee", lw=0, label="±1 log unit")
        ax.plot(lims, lims, color="#333333", lw=0.8, label="y = x")
        ax.scatter(d.observed, d.predicted, s=6, alpha=0.45, color=SPLIT_STYLE[split]["color"],
                   marker=SPLIT_STYLE[split]["marker"], lw=0)
        within = np.mean(np.abs(d.predicted - d.observed) <= 1)
        ax.text(0.04, 0.96, f"n = {len(d):,}\nRMSE = {m.rmse:.2f}\nR² = {m.r2:.2f}\n"
                            f"{within:.0%} within ±1 log unit",
                transform=ax.transAxes, va="top", fontsize=8)
        ax.set_title(f"{tag}  XGBoost, {split} split")
        ax.set(xlim=lims, ylim=lims, aspect="equal")
    axes[0].set_ylabel("Predicted logS (log10 mol/L)")
    fig.supxlabel("Observed logS (log10 mol/L)", fontsize=9)
    axes[1].legend(loc="lower right")
    save(fig, "fig2_predicted_vs_observed")


def fig_error_vs_similarity(ad):
    fig, ax = plt.subplots(figsize=(4.2, 3.0), layout="constrained")
    x = np.arange(len(ad))
    lo = ad.rmse - ad.rmse_ci.map(lambda c: c[0])
    hi = ad.rmse_ci.map(lambda c: c[1]) - ad.rmse
    ax.errorbar(x, ad.rmse, yerr=[lo, hi], fmt="s-", color=VERMILLION, ms=5, capsize=2.5, lw=1)
    for xi, (r, n) in enumerate(zip(ad.rmse, ad.n)):
        ax.annotate(f"n={n}", (xi, r), textcoords="offset points", xytext=(6, 4), fontsize=7, color="#555555")
    ax.set_xticks(x, ad.similarity_bin, rotation=0)
    ax.set_xlabel("Tanimoto similarity to nearest training molecule")
    ax.set_ylabel("Test RMSE (log units)")
    ax.set_ylim(bottom=0)
    ax.set_title("XGBoost, scaffold split (95% bootstrap CI)", fontsize=8.5)
    ax.grid(axis="y", lw=0.4, color="#dddddd")
    save(fig, "fig3_error_vs_similarity")


def fig_importance(imp, top=15):
    d = imp.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    ax.barh(d.feature, d.rmse_increase_mean, xerr=d.rmse_increase_sd, color=BLUE,
            error_kw=dict(lw=0.8, capsize=2))
    ax.set_xlabel("Increase in test RMSE when permuted (log units)")
    ax.set_title("XGBoost, scaffold split (mean ± SD, 5 repeats)", fontsize=8.5)
    ax.grid(axis="x", lw=0.4, color="#dddddd")
    save(fig, "fig4_permutation_importance")


def split_shift_table(preds):
    """Median and max of logS, MW and logP in train vs test for each split."""
    meta = pd.read_csv(ROOT / "data" / "processed" / "molecules.csv", keep_default_na=False)
    rows = []
    for split in ("random", "scaffold"):
        test_rows = preds[preds.split == split].row.unique()
        is_test = meta.index.isin(test_rows)
        for part, sub in (("train", meta[~is_test]), ("test", meta[is_test])):
            row = {"split": split, "set": part, "n": len(sub)}
            for col, label in (("logS", "logS"), ("esol_mw", "MW"), ("esol_logp", "cLogP")):
                v = sub[col].astype(float)
                row[f"{label}_median"], row[f"{label}_max"] = round(v.median(), 2), round(v.max(), 1)
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "split_shift.csv", index=False)
    return table


def main():
    FIGURES.mkdir(exist_ok=True)
    metrics = pd.read_csv(RESULTS / "metrics.csv")
    for c in ("rmse_ci", "r2_ci", "mae_ci"):
        metrics[c] = metrics[c].map(ast.literal_eval)
    preds = pd.read_csv(RESULTS / "test_predictions.csv")
    ad = pd.read_csv(RESULTS / "error_by_similarity_scaffold_xgb.csv")
    ad["rmse_ci"] = ad.rmse_ci.map(ast.literal_eval)
    imp = pd.read_csv(RESULTS / "permutation_importance_scaffold_xgb.csv")
    noise = json.loads((RESULTS / "curation_log.json").read_text())["noise_floor"]["rms_between_source_sd"]

    fig_model_comparison(metrics, noise)
    fig_pred_vs_obs(preds, metrics)
    fig_error_vs_similarity(ad)
    fig_importance(imp)
    print(split_shift_table(preds).to_string(index=False))
    print("figures written to", FIGURES)


if __name__ == "__main__":
    main()
