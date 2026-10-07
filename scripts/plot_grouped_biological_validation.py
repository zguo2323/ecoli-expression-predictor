"""Plot grouped-CV/bootstrap metrics for the biological-context report."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


LABELS = {
    "B4_sequence_only": "Sequence only",
    "B3_historical_MFE": "Historical MFE",
    "B6_TSS_accessibility": "TSS accessibility",
    "Ridge_B6_TSS_accessibility": "Ridge · TSS accessibility",
    "B7_OSTIR": "OSTIR",
}
COLORS = {
    "B4_sequence_only": "#6B7280",
    "B3_historical_MFE": "#D97706",
    "B6_TSS_accessibility": "#007F73",
    "Ridge_B6_TSS_accessibility": "#2563EB",
    "B7_OSTIR": "#7C3AED",
}
SCENARIOS = ["random", "promoter_holdout", "rbs_holdout", "double_unseen"]
SCENARIO_LABELS = ["Random", "Promoter", "RBS", "Double unseen"]


def plot_grouped_summary(ci_path: Path, output_path: Path) -> None:
    ci = pd.read_csv(ci_path)
    ci = ci.loc[~ci.baseline_or_delta.str.startswith("DELTA:")].copy()
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.5), sharex="col")
    settings = [
        ("log_translation_proxy", "spearman", "Translation proxy · Spearman ρ", 0.0, 1.0),
        ("log_translation_proxy", "mae_log", "Translation proxy · log-scale MAE", None, None),
        ("log_prot", "spearman", "Total protein · Spearman ρ", 0.0, 1.0),
        ("log_prot", "mae_log", "Total protein · log-scale MAE", None, None),
    ]
    x_base = range(len(SCENARIOS))
    offsets = {"B4_sequence_only": -0.24, "B3_historical_MFE": -0.12,
               "B6_TSS_accessibility": 0.0,
               "Ridge_B6_TSS_accessibility": 0.12, "B7_OSTIR": 0.24}
    handles = {}
    for ax, (target, metric, title, ymin, ymax) in zip(axes.flat, settings):
        sub = ci.loc[(ci.target == target) & (ci.metric == metric)]
        available = list(LABELS) if target == "log_translation_proxy" else list(LABELS)[:-1]
        for baseline in available:
            method = sub.loc[sub.baseline_or_delta == baseline].set_index("scenario")
            method = method.reindex(SCENARIOS)
            x = [v + offsets[baseline] for v in x_base]
            y = method.bootstrap_median.to_numpy(float)
            lo = y - method.ci_2_5.to_numpy(float)
            hi = method.ci_97_5.to_numpy(float) - y
            points = ax.errorbar(x, y, yerr=[lo, hi], fmt="o", markersize=5,
                                 capsize=3, linewidth=1.25, color=COLORS[baseline],
                                 label=LABELS[baseline])
            handles[baseline] = points
        ax.set_title(title, fontsize=11, loc="left")
        ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        if ymin is not None:
            ax.set_ylim(ymin, ymax)
        if metric == "mae_log":
            ax.set_ylabel("Lower is better")
        else:
            ax.set_ylabel("Higher is better")
    for ax in axes[-1]:
        ax.set_xticks(list(x_base), SCENARIO_LABELS)
    fig.legend([handles[key] for key in LABELS], [LABELS[key] for key in LABELS],
               loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.99))
    fig.suptitle("Biological-context models across held-out groups", y=1.035,
                 fontsize=15, fontweight="semibold")
    fig.text(0.01, 0.01,
             "Points and bars show bootstrap medians and 95% percentile intervals. "
             "Group resampling is used for biological holdouts; random split is diagnostic. "
             "OSTIR is evaluated on translation proxy only.", fontsize=8.7, color="#4B5563")
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ci", default="reports/biological_context/group_bootstrap_ci.csv")
    parser.add_argument("--output", default="reports/biological_context/grouped_validation_summary.png")
    args = parser.parse_args()
    plot_grouped_summary(Path(args.ci), Path(args.output))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
