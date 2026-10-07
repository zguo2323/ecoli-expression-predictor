"""Grouped OOF residual diagnostics for the biological-context models."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

PREDICTION_COLUMNS = {
    "scenario", "fold", "target", "baseline", "benchmark_row_id",
    "promoter_id", "rbs_id", "actual", "predicted", "residual",
}
MAIN_MODELS = [
    "B6_TSS_accessibility", "Ridge_B6_TSS_accessibility", "B7_OSTIR",
]
PAIRS = [
    ("B6_TSS_accessibility", "Ridge_B6_TSS_accessibility"),
    ("B6_TSS_accessibility", "B7_OSTIR"),
]


def add_actual_quantile_bins(predictions: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    """Add within-scenario/target actual-value bins shared across models."""
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")
    missing = PREDICTION_COLUMNS - set(predictions.columns)
    if missing:
        raise ValueError(f"OOF predictions missing columns: {sorted(missing)}")
    actuals = predictions[
        ["scenario", "target", "benchmark_row_id", "actual"]
    ].drop_duplicates(["scenario", "target", "benchmark_row_id"]).copy()
    actuals["actual_bin"] = 1
    for _, group in actuals.groupby(["scenario", "target"], sort=False):
        bins = pd.qcut(group.actual, q=n_bins, labels=False, duplicates="drop")
        actuals.loc[group.index, "actual_bin"] = bins.fillna(0).astype(int).to_numpy() + 1
    return predictions.merge(
        actuals[["scenario", "target", "benchmark_row_id", "actual_bin"]],
        on=["scenario", "target", "benchmark_row_id"], how="left", validate="many_to_one",
    )


def _residual_metrics(frame: pd.DataFrame) -> dict:
    residual = frame.actual.to_numpy(float) - frame.predicted.to_numpy(float)
    absolute = np.abs(residual)
    return {
        "n": len(frame),
        "mae_log": float(absolute.mean()),
        "median_absolute_error_log": float(np.median(absolute)),
        "rmse_log": float(np.sqrt(np.mean(residual ** 2))),
        "mean_signed_residual": float(residual.mean()),
        "median_signed_residual": float(np.median(residual)),
        "underprediction_fraction": float(np.mean(residual > 0)),
        "overprediction_fraction": float(np.mean(residual < 0)),
    }


def summarize_actual_bins(predictions: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    binned = add_actual_quantile_bins(predictions, n_bins=n_bins)
    rows = []
    for keys, group in binned.groupby(
        ["scenario", "target", "actual_bin", "baseline"], sort=False
    ):
        rows.append(dict(zip(
            ["scenario", "target", "actual_bin", "baseline"], keys
        ), **_residual_metrics(group)))
    return pd.DataFrame(rows)


def component_error_table(predictions: pd.DataFrame) -> pd.DataFrame:
    """Errors per promoter/RBS, plus per pair for double-unseen samples."""
    frame = predictions.copy()
    rows = []
    for component in ("promoter_id", "rbs_id"):
        for keys, group in frame.groupby(
            ["scenario", "target", "baseline", component], sort=False
        ):
            scenario, target, baseline, component_id = keys
            rows.append({
                "scenario": scenario, "target": target, "baseline": baseline,
                "component_type": component.removesuffix("_id"),
                "component_id": component_id,
                **_residual_metrics(group),
            })
    double = frame.loc[frame.scenario == "double_unseen"].copy()
    double["component_id"] = double.promoter_id.astype(str) + "::" + double.rbs_id.astype(str)
    for keys, group in double.groupby(
        ["scenario", "target", "baseline", "component_id"], sort=False
    ):
        scenario, target, baseline, component_id = keys
        rows.append({
            "scenario": scenario, "target": target, "baseline": baseline,
            "component_type": "promoter_rbs_pair", "component_id": component_id,
            **_residual_metrics(group),
        })
    return pd.DataFrame(rows)


def summarize_component_errors(component_errors: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["scenario", "target", "baseline", "component_type"]
    for group_keys, group in component_errors.groupby(keys, sort=False):
        mae = group.mae_log
        rows.append(dict(zip(keys, group_keys),
                         n_components=len(group),
                         n_rows=int(group.n.sum()),
                         median_component_mae=float(mae.median()),
                         q25_component_mae=float(mae.quantile(.25)),
                         q75_component_mae=float(mae.quantile(.75)),
                         max_component_mae=float(mae.max()),
                         median_component_bias=float(group.mean_signed_residual.median()),
                         fraction_components_underpredicting=float(
                             (group.mean_signed_residual > 0).mean()
                         )))
    return pd.DataFrame(rows)


def component_outcome_variance(predictions: pd.DataFrame) -> pd.DataFrame:
    """Marginal fraction of held-out target variance associated with each ID."""
    source = predictions.loc[
        predictions.baseline == "B6_TSS_accessibility",
        ["scenario", "target", "benchmark_row_id", "actual", "promoter_id", "rbs_id"],
    ].drop_duplicates(["scenario", "target", "benchmark_row_id"])
    rows = []
    for (scenario, target), frame in source.groupby(["scenario", "target"], sort=False):
        grand_mean = float(frame.actual.mean())
        total_ss = float(((frame.actual - grand_mean) ** 2).sum())
        for component in ("promoter_id", "rbs_id"):
            means = frame.groupby(component).actual.agg(["mean", "count"])
            between_ss = float((means["count"] * (means["mean"] - grand_mean) ** 2).sum())
            rows.append({
                "scenario": scenario, "target": target,
                "component_type": component.removesuffix("_id"),
                "n_rows": len(frame), "n_components": len(means),
                "marginal_eta_squared": between_ss / total_ss if total_ss else float("nan"),
                "sd_component_mean": float(means["mean"].std()),
                "median_component_n": float(means["count"].median()),
                "q25_component_n": float(means["count"].quantile(.25)),
                "q75_component_n": float(means["count"].quantile(.75)),
            })
    return pd.DataFrame(rows)


def paired_component_comparison(predictions: pd.DataFrame) -> pd.DataFrame:
    """Compare absolute errors on identical OOF rows (positive favors B6)."""
    rows = []
    identity = ["scenario", "fold", "target", "benchmark_row_id",
                "promoter_id", "rbs_id", "actual"]
    for left, right in PAIRS:
        available = predictions.loc[predictions.baseline.isin([left, right])]
        if available.empty:
            continue
        wide = available.pivot(index=identity, columns="baseline", values="predicted").reset_index()
        if left not in wide or right not in wide:
            continue
        wide = wide.dropna(subset=[left, right])
        wide["left_abs_error"] = (wide.actual - wide[left]).abs()
        wide["right_abs_error"] = (wide.actual - wide[right]).abs()
        wide["right_minus_left_abs_error"] = wide.right_abs_error - wide.left_abs_error
        wide["left_better"] = wide.left_abs_error < wide.right_abs_error
        components = [("promoter", "promoter_id"), ("rbs", "rbs_id")]
        pair_rows = wide.loc[wide.scenario == "double_unseen"].copy()
        pair_rows["pair_id"] = pair_rows.promoter_id.astype(str) + "::" + pair_rows.rbs_id.astype(str)
        if not pair_rows.empty:
            components.append(("promoter_rbs_pair", "pair_id"))
        for component_type, component_column in components:
            source = pair_rows if component_type == "promoter_rbs_pair" else wide
            if source.empty:
                continue
            for keys, group in source.groupby(
                ["scenario", "target", component_column], sort=False
            ):
                scenario, target, component_id = keys
                rows.append({
                    "scenario": scenario, "target": target,
                    "left_model": left, "right_model": right,
                    "component_type": component_type, "component_id": component_id,
                    "n": len(group),
                    "left_mae_log": float(group.left_abs_error.mean()),
                    "right_mae_log": float(group.right_abs_error.mean()),
                    "right_minus_left_mae_log": float(group.right_minus_left_abs_error.mean()),
                    "fraction_left_better": float(group.left_better.mean()),
                })
    return pd.DataFrame(rows)


def summarize_paired_components(paired: pd.DataFrame) -> pd.DataFrame:
    """Summarize paired MAE deltas across biological component IDs."""
    keys = ["scenario", "target", "left_model", "right_model", "component_type"]
    rows = []
    for group_keys, group in paired.groupby(keys, sort=False):
        delta = group.right_minus_left_mae_log
        rows.append(dict(zip(keys, group_keys),
                         n_components=len(group),
                         median_component_delta_mae=float(delta.median()),
                         q25_component_delta_mae=float(delta.quantile(.25)),
                         q75_component_delta_mae=float(delta.quantile(.75)),
                         fraction_components_left_better=float((delta > 0).mean()),
                         median_fraction_rows_left_better=float(group.fraction_left_better.median())))
    return pd.DataFrame(rows)


def select_extreme_cases(predictions: pd.DataFrame, n_cases: int = 20) -> pd.DataFrame:
    """Return the largest absolute residuals per scenario/target/model."""
    if n_cases < 1:
        raise ValueError("n_cases must be positive")
    frame = predictions.loc[predictions.baseline.isin(MAIN_MODELS)].copy()
    frame["absolute_residual"] = frame.residual.abs()
    frame["error_direction"] = np.where(
        frame.residual > 0, "underprediction",
        np.where(frame.residual < 0, "overprediction", "exact"),
    )
    return (frame.sort_values("absolute_residual", ascending=False)
            .groupby(["scenario", "target", "baseline"], sort=False)
            .head(n_cases)
            .reset_index(drop=True))


def run_analysis(predictions: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build reproducible residual, component, paired, and extreme-case tables."""
    missing = PREDICTION_COLUMNS - set(predictions.columns)
    if missing:
        raise ValueError(f"OOF predictions missing columns: {sorted(missing)}")
    errors = predictions.loc[predictions.baseline.isin(MAIN_MODELS)].copy()
    bins = summarize_actual_bins(errors)
    components = component_error_table(errors)
    paired = paired_component_comparison(errors)
    extremes = select_extreme_cases(errors)
    return {
        "actual_bins": bins,
        "component_errors": components,
        "component_summary": summarize_component_errors(components),
        "target_component_variance": component_outcome_variance(errors),
        "paired_component_comparison": paired,
        "paired_component_summary": summarize_paired_components(paired),
        "extreme_cases": extremes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictions", default="reports/biological_context/grouped_oof_predictions.csv"
    )
    parser.add_argument("--output-dir", default="reports/biological_context")
    parser.add_argument("--extreme-cases", type=int, default=20)
    args = parser.parse_args()
    predictions = pd.read_csv(args.predictions)
    reports = run_analysis(predictions)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    names = {
        "actual_bins": "grouped_oof_error_by_actual_bin.csv",
        "component_errors": "grouped_oof_error_by_component.csv",
        "component_summary": "grouped_oof_error_component_summary.csv",
        "target_component_variance": "grouped_oof_target_variance_by_component.csv",
        "paired_component_comparison": "grouped_oof_error_paired_component_comparison.csv",
        "paired_component_summary": "grouped_oof_error_paired_component_summary.csv",
        "extreme_cases": "grouped_oof_error_extreme_cases.csv",
    }
    for key, frame in reports.items():
        if key == "extreme_cases":
            frame = select_extreme_cases(predictions, args.extreme_cases)
        frame.to_csv(output / names[key], index=False)
    print("Wrote grouped OOF residual diagnostics:")
    for key, frame in reports.items():
        print(f"  {names[key]}: {len(frame)} rows")
    print("\nRBS-holdout translation-proxy component comparison (B6 vs Ridge):")
    paired_summary = reports["paired_component_summary"]
    print(paired_summary.loc[
        (paired_summary.scenario == "rbs_holdout")
        & (paired_summary.target == "log_translation_proxy")
        & (paired_summary.right_model == "Ridge_B6_TSS_accessibility")
    ][["component_type", "n_components", "median_component_delta_mae",
       "fraction_components_left_better", "median_fraction_rows_left_better"]]
      .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
