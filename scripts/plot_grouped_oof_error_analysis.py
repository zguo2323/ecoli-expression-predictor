"""Plot the RBS-holdout OOF error diagnostics."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

MODEL_LABELS = {
    "B6_TSS_accessibility": "XGBoost B6",
    "Ridge_B6_TSS_accessibility": "Ridge B6",
    "B7_OSTIR": "OSTIR",
}
MODEL_COLORS = {
    "B6_TSS_accessibility": "#007F73",
    "Ridge_B6_TSS_accessibility": "#2563EB",
    "B7_OSTIR": "#7C3AED",
}


def plot_error_diagnostics(
    paired_path: Path, bins_path: Path, output_path: Path
) -> None:
    paired = pd.read_csv(paired_path)
    bins = pd.read_csv(bins_path)
    comparison = paired.loc[
        (paired.scenario == "rbs_holdout")
        & (paired.target == "log_translation_proxy")
        & (paired.left_model == "B6_TSS_accessibility")
        & (paired.right_model == "Ridge_B6_TSS_accessibility")
    ]
    bin_data = bins.loc[
        (bins.scenario == "rbs_holdout")
        & (bins.target == "log_translation_proxy")
        & bins.baseline.isin(MODEL_LABELS)
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.0))
    for ax, component, title in (
        (axes[0, 0], "rbs", "Per unseen RBS"),
        (axes[0, 1], "promoter", "Per known promoter background"),
    ):
        subset = comparison.loc[comparison.component_type == component]
        x = subset.left_mae_log.to_numpy(float)
        y = subset.right_mae_log.to_numpy(float)
        high = max(float(x.max()), float(y.max())) * 1.04
        ax.scatter(x, y, s=23, alpha=0.72, color="#007F73", edgecolors="none")
        ax.plot([0, high], [0, high], linestyle="--", linewidth=1,
                color="#6B7280", alpha=0.8)
        ax.set_xlim(0, high)
        ax.set_ylim(0, high)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("XGBoost B6 group MAE")
        ax.set_ylabel("Ridge B6 group MAE")
        ax.set_title(title, loc="left", fontsize=11)
        ax.grid(color="#D1D5DB", linewidth=0.6, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(
            0.03, 0.96,
            f"{(y > x).sum()}/{len(subset)} groups above diagonal favor B6",
            transform=ax.transAxes, va="top", fontsize=8.5, color="#374151",
        )
        if component == "rbs":
            row = subset.nsmallest(1, "right_minus_left_mae_log").iloc[0]
            ax.annotate(
                str(row.component_id), (row.left_mae_log, row.right_mae_log),
                xytext=(12, -18), textcoords="offset points", fontsize=8,
                arrowprops={"arrowstyle": "-", "color": "#4B5563", "lw": 0.7},
            )

    for ax, metric, title, ylabel in (
        (axes[1, 0], "mean_signed_residual",
         "Signed residual by measured-expression quintile", "Actual − predicted (log units)"),
        (axes[1, 1], "mae_log",
         "Absolute error by measured-expression quintile", "Log-scale MAE"),
    ):
        for model, label in MODEL_LABELS.items():
            subset = (bin_data.loc[bin_data.baseline == model]
                      .sort_values("actual_bin"))
            ax.plot(
                subset.actual_bin, subset[metric], marker="o", linewidth=1.8,
                markersize=4.5, label=label, color=MODEL_COLORS[model],
            )
        ax.axhline(0, color="#6B7280", linewidth=0.8, alpha=0.7)
        ax.set_xticks([1, 2, 3, 4, 5], ["Q1 low", "Q2", "Q3", "Q4", "Q5 high"])
        ax.set_xlabel("Measured translation-proxy quintile")
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontsize=11)
        ax.grid(axis="y", color="#D1D5DB", linewidth=0.6, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1, 1].legend(frameon=False, fontsize=8, loc="best")

    fig.suptitle(
        "RBS-holdout OOF diagnostics: component effects and calibration bias",
        y=1.02, fontsize=14, fontweight="semibold",
    )
    fig.text(
        0.01, 0.01,
        "Top: each point is one held-out component; points above the diagonal have lower MAE for B6. "
        "Bottom: quintiles are defined within the RBS-holdout test set. "
        "Negative residual means overprediction; positive means underprediction.",
        fontsize=8.2, color="#4B5563",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--paired",
        default="reports/biological_context/grouped_oof_error_paired_component_comparison.csv",
    )
    parser.add_argument(
        "--actual-bins",
        default="reports/biological_context/grouped_oof_error_by_actual_bin.csv",
    )
    parser.add_argument(
        "--output",
        default="reports/biological_context/grouped_oof_error_analysis.png",
    )
    args = parser.parse_args()
    plot_error_diagnostics(Path(args.paired), Path(args.actual_bins), Path(args.output))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
