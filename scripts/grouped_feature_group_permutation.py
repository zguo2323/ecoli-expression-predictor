"""Held-out, grouped-feature permutation importance for B6 across outer folds.

Each biological feature block is permuted jointly on a fold's held-out rows.
This preserves dependencies within a block while breaking its row-level
association with the target. It measures model reliance, not a causal effect.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.biological_benchmark import ACCESS_FEATURES, BASIC_FEATURES, _fit_xgb, _metrics

FEATURE_SET = BASIC_FEATURES + ACCESS_FEATURES
FEATURE_GROUPS = {
    "promoter_recognition": [
        "gc_promoter", "score_minus10", "score_minus35",
        "spacer_length", "spacer_optimal",
    ],
    "rbs_recognition": [
        "gc_rbs", "score_sd", "sd_spacing", "sd_spacing_optimal",
    ],
    "transcript_accessibility": [
        "sd_unpaired_probability", "start_unpaired_probability",
        "sd_start_opening_energy",
    ],
}
PERMUTATION_UNITS = {
    "promoter_recognition": ["promoter_id"],
    "rbs_recognition": ["rbs_id"],
    "transcript_accessibility": ["promoter_id", "rbs_id"],
}


def validate_feature_groups(groups=FEATURE_GROUPS, features=FEATURE_SET) -> None:
    flattened = [feature for group in groups.values() for feature in group]
    if len(flattened) != len(set(flattened)):
        raise ValueError("A feature may belong to only one permutation group")
    if set(flattened) != set(features):
        missing = sorted(set(features) - set(flattened))
        extra = sorted(set(flattened) - set(features))
        raise ValueError(f"Feature groups must partition B6 features; missing={missing}, extra={extra}")


def permute_feature_group(frame: pd.DataFrame, features: list[str],
                          rng: np.random.Generator,
                          unit_columns: list[str] | None = None) -> pd.DataFrame:
    """Jointly permute a block across biological units, preserving repeats.

    Feature values must be constant within each unit. The same source unit is
    assigned to every row belonging to a destination unit.
    """
    if unit_columns is None:
        unit_columns = []
    permuted = frame.copy()
    if not unit_columns:
        order = rng.permutation(len(frame))
        permuted.loc[:, features] = frame.iloc[order][features].to_numpy()
        return permuted

    missing = set(unit_columns + features) - set(frame.columns)
    if missing:
        raise ValueError(f"Permutation frame missing columns: {sorted(missing)}")
    distinct_counts = frame.groupby(unit_columns, dropna=False)[features].nunique(dropna=False)
    if (distinct_counts > 1).any().any():
        raise ValueError("Feature group values must be constant within permutation units")
    units = frame[unit_columns].drop_duplicates().reset_index(drop=True)
    values = (frame.groupby(unit_columns, sort=False, dropna=False)[features]
              .first().reset_index())
    unit_index = pd.MultiIndex.from_frame(units[unit_columns])
    row_index = pd.MultiIndex.from_frame(frame[unit_columns])
    row_codes = unit_index.get_indexer(row_index)
    order = rng.permutation(len(units))
    source_values = values[features].to_numpy()
    permuted.loc[:, features] = source_values[order[row_codes]]
    return permuted


def permutation_metric_deltas(actual, baseline_prediction, permuted_prediction):
    """Positive deltas mean the original feature group helped prediction."""
    baseline = _metrics(np.asarray(actual, dtype=float),
                        np.asarray(baseline_prediction, dtype=float))
    permuted = _metrics(np.asarray(actual, dtype=float),
                        np.asarray(permuted_prediction, dtype=float))
    return {
        "delta_spearman": baseline["spearman"] - permuted["spearman"],
        "delta_mae_log": permuted["mae_log"] - baseline["mae_log"],
        "baseline_spearman": baseline["spearman"],
        "permuted_spearman": permuted["spearman"],
        "baseline_mae_log": baseline["mae_log"],
        "permuted_mae_log": permuted["mae_log"],
    }


def run_permutation_importance(
    data: pd.DataFrame,
    assignments: pd.DataFrame,
    n_permutations: int = 30,
    random_state: int = 42,
):
    """Fit B6 per saved outer fold and permute feature groups on its test rows."""
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive")
    validate_feature_groups()
    required = set(FEATURE_SET + ["benchmark_row_id", "log_prot",
                                  "log_translation_proxy"])
    if not required.issubset(data.columns):
        raise ValueError(f"Validation data missing columns: {sorted(required - set(data.columns))}")
    needed_assignment_columns = {"scenario", "fold", "benchmark_row_id", "role"}
    if not needed_assignment_columns.issubset(assignments.columns):
        raise ValueError(
            f"Fold assignments missing columns: "
            f"{sorted(needed_assignment_columns - set(assignments.columns))}"
        )
    if data.benchmark_row_id.duplicated().any():
        raise ValueError("benchmark_row_id must be unique in validation data")

    indexed = data.set_index("benchmark_row_id", drop=False)
    rows = []
    scenario_order = list(assignments.scenario.drop_duplicates())
    for (scenario, fold), split in assignments.groupby(["scenario", "fold"], sort=False):
        train_ids = split.loc[split.role == "train", "benchmark_row_id"].astype(int).tolist()
        test_ids = split.loc[split.role == "test", "benchmark_row_id"].astype(int).tolist()
        if not train_ids or not test_ids:
            continue
        absent = (set(train_ids) | set(test_ids)) - set(indexed.index)
        if absent:
            raise ValueError(f"Fold assignments refer to missing benchmark IDs: {sorted(absent)[:5]}")
        train, test = indexed.loc[train_ids], indexed.loc[test_ids]
        for target_index, target in enumerate(("log_translation_proxy", "log_prot")):
            actual = test[target].to_numpy(float)
            baseline_prediction, model = _fit_xgb(train, test, FEATURE_SET, target)
            for group_index, (group_name, features) in enumerate(FEATURE_GROUPS.items()):
                for repeat in range(n_permutations):
                    seed = (random_state + 100_000 * group_index
                            + 10_000 * target_index + 1_000 * int(fold) + repeat
                            + 1_000_000 * scenario_order.index(scenario))
                    permuted = permute_feature_group(
                        test, features, np.random.default_rng(seed),
                        PERMUTATION_UNITS[group_name],
                    )
                    permuted_prediction = model.predict(permuted[FEATURE_SET].astype(float))
                    delta = permutation_metric_deltas(
                        actual, baseline_prediction, permuted_prediction
                    )
                    rows.append({
                        "scenario": scenario,
                        "fold": int(fold),
                        "target": target,
                        "feature_group": group_name,
                        "n_test": len(test),
                        "permutation": repeat,
                        "permutation_seed": seed,
                        **delta,
                    })
    return pd.DataFrame(rows)


def summarize_permutation_importance(permutations: pd.DataFrame):
    """Report permutation-repeat mean per fold, then equal-fold mean ± SD."""
    fold = (
        permutations.groupby(
            ["scenario", "fold", "target", "feature_group"], sort=False
        )
        .agg(
            delta_spearman=("delta_spearman", "mean"),
            permutation_sd_spearman=("delta_spearman", "std"),
            delta_mae_log=("delta_mae_log", "mean"),
            permutation_sd_mae_log=("delta_mae_log", "std"),
            n_test=("n_test", "first"),
            n_permutations=("permutation", "nunique"),
        )
        .reset_index()
    )
    summary = (
        fold.groupby(["scenario", "target", "feature_group"], sort=False)
        .agg(
            mean_delta_spearman=("delta_spearman", "mean"),
            sd_across_folds_spearman=("delta_spearman", "std"),
            mean_delta_mae_log=("delta_mae_log", "mean"),
            sd_across_folds_mae_log=("delta_mae_log", "std"),
            mean_permutation_sd_spearman=("permutation_sd_spearman", "mean"),
            mean_permutation_sd_mae_log=("permutation_sd_mae_log", "mean"),
            n_folds=("fold", "nunique"),
            n_test_total=("n_test", "sum"),
            n_permutations=("n_permutations", "min"),
        )
        .reset_index()
    )
    return fold, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="reports/biological_context/validation_dataset.parquet")
    parser.add_argument(
        "--fold-assignments",
        default="reports/biological_context/grouped_fold_assignments.csv",
    )
    parser.add_argument("--output-dir", default="reports/biological_context")
    parser.add_argument("--permutations", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    data = pd.read_parquet(args.data)
    assignments = pd.read_csv(args.fold_assignments)
    permutations = run_permutation_importance(
        data, assignments, args.permutations, args.random_state
    )
    fold, summary = summarize_permutation_importance(permutations)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    permutations.to_csv(output / "grouped_feature_group_permutation.csv", index=False)
    fold.to_csv(output / "grouped_feature_group_permutation_by_fold.csv", index=False)
    summary.to_csv(output / "grouped_feature_group_permutation_summary.csv", index=False)
    metadata = {
        "n_permutations_per_group_fold_target": args.permutations,
        "random_state": args.random_state,
        "feature_groups": FEATURE_GROUPS,
        "permutation_units": PERMUTATION_UNITS,
        "outer_folds": int(fold.groupby("scenario").fold.nunique().max()),
        "aggregation": (
            "Joint permutation within each feature group at its biological unit "
            "(promoter, RBS, or promoter-RBS pair); repeated rows within a unit "
            "receive the same values. "
            "permutation repeats averaged within fold, followed by equal-fold "
            "mean and SD across the fixed outer folds."
        ),
        "interpretation": (
            "Positive delta_spearman and delta_mae_log indicate that permuting "
            "the group degraded predictions. This is model reliance, not causality."
        ),
    }
    (output / "grouped_feature_group_permutation_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(summary.to_string(index=False))
    print(f"\nWrote grouped feature permutation results to {output}")


if __name__ == "__main__":
    main()
