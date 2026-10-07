"""Plot held-out B6 feature-group permutation importance."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

SCENARIOS = ["random", "promoter_holdout", "rbs_holdout", "double_unseen"]
SCENARIO_LABELS = ["Random", "Promoter", "RBS", "Double unseen"]
GROUP_LABELS = {
    "promoter_recognition": "Promoter recognition",
    "rbs_recognition": "RBS / SD recognition",
    "transcript_accessibility": "Transcript accessibility",
}
COLORS = {
    "promoter_recognition": "#D97706",
    "rbs_recognition": "#007F73",
    "transcript_accessibility": "#2563EB",
}


def plot_summary(summary_path: Path, output_path: Path) -> None:
    data = pd.read_csv(summary_path)
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.5), sharex="col")
    settings = [
        ("log_translation_proxy", "mean_delta_spearman",
         "Translation proxy · Spearman drop", "Higher means greater model reliance"),
        ("log_translation_proxy", "mean_delta_mae_log",
         "Translation proxy · log-MAE increase", "Higher means greater model reliance"),
        ("log_prot", "mean_delta_spearman",
         "Total protein · Spearman drop", "Higher means greater model reliance"),
        ("log_prot", "mean_delta_mae_log",
         "Total protein · log-MAE increase", "Higher means greater model reliance"),
    ]
    x = range(len(SCENARIOS))
    offsets = {
        "promoter_recognition": -0.16,
        "rbs_recognition": 0.0,
        "transcript_accessibility": 0.16,
    }
    handles = {}
    for ax, (target, metric, title, ylabel) in zip(axes.flat, settings):
        subset = data.loc[data.target == target]
        for group, label in GROUP_LABELS.items():
            rows = (subset.loc[subset.feature_group == group]
                    .set_index("scenario").reindex(SCENARIOS))
            y = rows[metric].to_numpy(float)
            sd_col = ("sd_across_folds_spearman" if metric.endswith("spearman")
                      else "sd_across_folds_mae_log")
            err = rows[sd_col].to_numpy(float)
            points = ax.errorbar(
                [index + offsets[group] for index in x], y, yerr=err,
                fmt="o", markersize=5, capsize=3, linewidth=1.25,
                color=COLORS[group], label=label,
            )
            handles[group] = points
        ax.axhline(0, color="#4B5563", linewidth=0.8, alpha=0.7)
        ax.set_title(title, fontsize=11, loc="left")
        ax.set_ylabel(ylabel, fontsize=8.5)
        ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[-1]:
        ax.set_xticks(list(x), SCENARIO_LABELS)
    fig.legend([handles[key] for key in GROUP_LABELS],
               [GROUP_LABELS[key] for key in GROUP_LABELS],
               loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.99))
    fig.suptitle("B6 held-out feature-group permutation importance", y=1.035,
                 fontsize=15, fontweight="semibold")
    fig.text(
        0.01, 0.01,
        "Points are means across five outer folds; bars are ±1 SD across folds. "
        "Each feature block was jointly row-permuted 30 times within each test fold. "
        "Importance measures model reliance, not causality.",
        fontsize=8.3, color="#4B5563",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        default="reports/biological_context/grouped_feature_group_permutation_summary.csv",
    )
    parser.add_argument(
        "--output",
        default="reports/biological_context/grouped_feature_group_permutation.png",
    )
    args = parser.parse_args()
    plot_summary(Path(args.summary), Path(args.output))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
