"""Grouped cross-validation and group-bootstrap evaluation for B3/B4/B6/B7.

Run OSTIR first with scripts/run_ostir_baseline.py.  This script can reuse its
validation_dataset.parquet and ostir_predictions.csv artifacts.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.biological_benchmark import ACCESS_FEATURES, BASIC_FEATURES, _metrics, _fit_xgb, prepare_benchmark_data

SCENARIOS = ("random", "promoter_holdout", "rbs_holdout", "double_unseen")
TARGETS = ("log_translation_proxy", "log_prot")
FEATURE_SETS = {
    "B4_sequence_only": BASIC_FEATURES,
    "B3_historical_MFE": BASIC_FEATURES + ["legacy_promoter_tail_mfe"],
    "B6_TSS_accessibility": BASIC_FEATURES + ACCESS_FEATURES,
}
RIDGE_ALPHAS = np.logspace(-5, 7, 25)


def make_inner_cv_splits(train: pd.DataFrame, scenario: str, n_splits: int = 3,
                         random_state: int = 42):
    """Create training-only folds that follow the outer scenario's grouping."""
    n = len(train)
    if n_splits < 2 or n < n_splits:
        raise ValueError("Need at least n_splits training rows")
    if scenario == "random":
        return list(KFold(n_splits, shuffle=True, random_state=random_state).split(np.arange(n)))
    p_ids = np.asarray(sorted(train.promoter_id.dropna().unique()))
    r_ids = np.asarray(sorted(train.rbs_id.dropna().unique()))

    def folds_for(ids):
        if len(ids) < n_splits:
            raise ValueError("Not enough unique groups for grouped inner validation")
        result = {}
        for fold, (_, validation) in enumerate(
            KFold(n_splits, shuffle=True, random_state=random_state).split(ids)
        ):
            for group in ids[validation]:
                result[group] = fold
        return result

    p_fold, r_fold = folds_for(p_ids), folds_for(r_ids)
    row_p = train.promoter_id.map(p_fold).to_numpy()
    row_r = train.rbs_id.map(r_fold).to_numpy()
    splits = []
    for fold in range(n_splits):
        p_held, r_held = row_p == fold, row_r == fold
        if scenario == "promoter_holdout":
            validation, inner_train = p_held, ~p_held
        elif scenario == "rbs_holdout":
            validation, inner_train = r_held, ~r_held
        elif scenario == "double_unseen":
            validation, inner_train = p_held & r_held, ~p_held & ~r_held
        else:
            raise ValueError(f"Unknown split scenario: {scenario}")
        train_idx, validation_idx = np.flatnonzero(inner_train), np.flatnonzero(validation)
        if len(train_idx) and len(validation_idx):
            splits.append((train_idx, validation_idx))
    if len(splits) < 2:
        raise ValueError(f"Not enough non-empty inner folds for {scenario}")
    return splits


def _fit_ridge(train: pd.DataFrame, test: pd.DataFrame, features: list[str],
               target: str, scenario: str, random_state: int):
    """Tune standardized Ridge only inside the current outer training fold."""
    pipeline = make_pipeline(StandardScaler(), Ridge())
    search = GridSearchCV(
        pipeline,
        {"ridge__alpha": RIDGE_ALPHAS},
        scoring="neg_mean_absolute_error",
        cv=make_inner_cv_splits(train, scenario, n_splits=3, random_state=random_state),
        n_jobs=1,
        refit=True,
        error_score="raise",
    )
    search.fit(train[features].astype(float), train[target].astype(float))
    prediction = search.predict(test[features].astype(float))
    best = search.best_estimator_
    coefs = best.named_steps["ridge"].coef_
    return prediction, float(search.best_params_["ridge__alpha"]), coefs


def make_group_cv_splits(df: pd.DataFrame, n_splits: int = 5, random_state: int = 42):
    """Build identical, deterministic outer folds for row and group holdouts.

    In double-unseen folds, training rows touching either held-out group are
    excluded. Test rows are the intersection of the held-out promoter and RBS
    groups; cross-pairs are intentionally unused in that fold.
    """
    if n_splits < 2:
        raise ValueError("n_splits must be at least 2")
    n = len(df)
    if n == 0:
        raise ValueError("Cannot split an empty dataframe")
    scenarios: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {s: [] for s in SCENARIOS}

    row_folds = np.empty(n, dtype=int)
    for fold, (_, test) in enumerate(KFold(n_splits, shuffle=True, random_state=random_state).split(np.arange(n))):
        row_folds[test] = fold
    p_ids = np.asarray(sorted(df["promoter_id"].dropna().unique()))
    r_ids = np.asarray(sorted(df["rbs_id"].dropna().unique()))
    if len(p_ids) < n_splits or len(r_ids) < n_splits:
        raise ValueError("Each grouping variable must have at least n_splits unique IDs")

    def group_fold_map(ids):
        fold_map = {}
        for fold, (_, held) in enumerate(KFold(n_splits, shuffle=True, random_state=random_state).split(ids)):
            for value in ids[held]:
                fold_map[value] = fold
        return fold_map

    p_fold, r_fold = group_fold_map(p_ids), group_fold_map(r_ids)
    row_p = df["promoter_id"].map(p_fold).to_numpy()
    row_r = df["rbs_id"].map(r_fold).to_numpy()
    for fold in range(n_splits):
        random_test = np.flatnonzero(row_folds == fold)
        scenarios["random"].append((np.flatnonzero(row_folds != fold), random_test))

        p_test = row_p == fold
        r_test = row_r == fold
        scenarios["promoter_holdout"].append((np.flatnonzero(~p_test), np.flatnonzero(p_test)))
        scenarios["rbs_holdout"].append((np.flatnonzero(~r_test), np.flatnonzero(r_test)))
        double_test = p_test & r_test
        double_train = ~p_test & ~r_test
        scenarios["double_unseen"].append((np.flatnonzero(double_train), np.flatnonzero(double_test)))

    assignments = []
    source_row_ids = (df["benchmark_row_id"].astype(int).to_numpy()
                      if "benchmark_row_id" in df else np.arange(n))
    for scenario, folds in scenarios.items():
        for fold, (train, test) in enumerate(folds):
            train_set, test_set = set(train), set(test)
            for row_index in range(n):
                role = "train" if row_index in train_set else "test" if row_index in test_set else "unused"
                assignments.append({"scenario": scenario, "fold": fold,
                                    "row_index": row_index,
                                    "benchmark_row_id": int(source_row_ids[row_index]),
                                    "role": role})
    return scenarios, pd.DataFrame(assignments)


def run_grouped_cv(df: pd.DataFrame, ostir: pd.DataFrame, n_splits: int = 5,
                   random_state: int = 42):
    """Fit sequence models inside each fold and evaluate held-out OSTIR too."""
    df = df.reset_index(drop=True).copy()
    ostir = ostir.copy()
    if "row_index" not in ostir or "ostir_expression" not in ostir:
        raise ValueError("OSTIR table must contain row_index and ostir_expression")
    row_ids = (df["benchmark_row_id"].astype(int).to_numpy()
               if "benchmark_row_id" in df else np.arange(len(df)))
    if len(np.unique(row_ids)) != len(row_ids) or ostir["row_index"].duplicated().any():
        raise ValueError("Benchmark and OSTIR row IDs must be unique")
    ostir_by_id = ostir.set_index("row_index")
    available_ids = set(ostir_by_id.index.astype(int))
    ostir_by_id.index = ostir_by_id.index.astype(int)
    keep = np.isin(row_ids, list(available_ids))
    df = df.loc[keep].reset_index(drop=True)
    row_ids = row_ids[keep]
    if len(df) < n_splits:
        raise ValueError("Too few shared scorable constructs for the requested number of folds")
    df["ostir_expression"] = ostir_by_id.loc[row_ids, "ostir_expression"].to_numpy()
    df["log_ostir"] = np.log(df["ostir_expression"].astype(float))
    splits, assignments = make_group_cv_splits(df, n_splits, random_state)
    rows, importance, ridge_coefficients = [], [], []
    for scenario, folds in splits.items():
        for fold, (train_idx, test_idx) in enumerate(folds):
            if not len(train_idx) or not len(test_idx):
                continue
            train, test = df.iloc[train_idx], df.iloc[test_idx]
            for target in TARGETS:
                actual = test[target].to_numpy(float)
                model_predictions: dict[str, np.ndarray] = {}
                model_alpha: dict[str, float] = {}
                for baseline, features in FEATURE_SETS.items():
                    pred, model = _fit_xgb(train, test, features, target)
                    model_predictions[baseline] = pred
                    if baseline == "B6_TSS_accessibility":
                        importance.extend({
                            "scenario": scenario, "fold": fold, "target": target,
                            "feature": feature, "importance": float(value),
                        } for feature, value in zip(features, model.feature_importances_))
                    ridge_name = f"Ridge_{baseline}"
                    ridge_pred, selected_alpha, coefficients = _fit_ridge(
                        train, test, features, target, scenario,
                        random_state + fold,
                    )
                    model_predictions[ridge_name] = ridge_pred
                    model_alpha[ridge_name] = selected_alpha
                    ridge_coefficients.extend({
                        "scenario": scenario, "fold": fold, "target": target,
                        "feature_set": baseline, "feature": feature,
                        "standardized_coefficient": float(coef),
                        "selected_alpha": selected_alpha,
                    } for feature, coef in zip(features, coefficients))

                if target == "log_translation_proxy":
                    # OSTIR is an external mechanistic predictor in different units.
                    # Fit its affine log-scale calibration on training rows only.
                    calibrator = LinearRegression().fit(train[["log_ostir"]], train[target])
                    model_predictions["B7_OSTIR"] = calibrator.predict(test[["log_ostir"]])

                for baseline, predicted in model_predictions.items():
                    for i, row_index in enumerate(test_idx):
                        rows.append({
                            "scenario": scenario, "fold": fold, "target": target,
                            "baseline": baseline, "row_index": int(row_index),
                            "benchmark_row_id": int(test.iloc[i].get("benchmark_row_id", row_index)),
                            "promoter_id": test.iloc[i]["promoter_id"],
                            "rbs_id": test.iloc[i]["rbs_id"],
                            "actual": float(actual[i]), "predicted": float(predicted[i]),
                            "residual": float(actual[i] - predicted[i]),
                            "ostir_expression": float(test.iloc[i]["ostir_expression"]),
                            "selected_alpha": model_alpha.get(baseline, np.nan),
                        })
    coverage = {"n_before_ostir_complete_case_filter": int(len(keep)),
                "n_common_scored": int(len(df)),
                "n_ostir_unscored": int((~keep).sum()),
                "coverage": float(keep.mean())}
    return (pd.DataFrame(rows), assignments, pd.DataFrame(importance),
            pd.DataFrame(ridge_coefficients), coverage)


def _resample_group_rows(frame: pd.DataFrame, scenario: str, rng: np.random.Generator) -> pd.DataFrame:
    if scenario == "promoter_holdout":
        group_col = "promoter_id"
        groups = frame[group_col].drop_duplicates().to_numpy()
        draws = rng.choice(groups, size=len(groups), replace=True)
        return pd.concat([frame.loc[frame[group_col] == group] for group in draws], ignore_index=True)
    if scenario == "rbs_holdout":
        group_col = "rbs_id"
        groups = frame[group_col].drop_duplicates().to_numpy()
        draws = rng.choice(groups, size=len(groups), replace=True)
        return pd.concat([frame.loc[frame[group_col] == group] for group in draws], ignore_index=True)
    if scenario == "double_unseen":
        promoters = frame.promoter_id.drop_duplicates().to_numpy()
        rbs_ids = frame.rbs_id.drop_duplicates().to_numpy()
        p_draw = rng.choice(promoters, size=len(promoters), replace=True)
        r_draw = rng.choice(rbs_ids, size=len(rbs_ids), replace=True)
        p_counts = pd.Series(p_draw).value_counts()
        r_counts = pd.Series(r_draw).value_counts()
        weight = (frame.promoter_id.map(p_counts).fillna(0).to_numpy()
                  * frame.rbs_id.map(r_counts).fillna(0).to_numpy()).astype(int)
        indices = np.repeat(np.arange(len(frame)), weight)
        return frame.iloc[indices].reset_index(drop=True)
    # Row bootstrap is a diagnostic only; grouped scenarios above are primary.
    return frame.iloc[rng.integers(0, len(frame), len(frame))].reset_index(drop=True)


def run_group_bootstrap(predictions: pd.DataFrame, n_boot: int = 1000,
                        random_state: int = 420):
    """Percentile CIs for metrics and paired deltas using biological groups."""
    if n_boot < 1:
        raise ValueError("n_boot must be positive")
    output = []
    rng = np.random.default_rng(random_state)
    for (scenario, target), frame in predictions.groupby(["scenario", "target"], sort=False):
        models = sorted(frame.baseline.unique())
        sampled: dict[str, dict[str, list[float]]] = {
            model: {"spearman": [], "mae_log": []} for model in models
        }
        deltas: dict[tuple[str, str], dict[str, list[float]]] = {}
        pairs = []
        if "B4_sequence_only" in models and "B6_TSS_accessibility" in models:
            pairs.append(("B6_TSS_accessibility", "B4_sequence_only"))
        if "B7_OSTIR" in models and "B6_TSS_accessibility" in models:
            pairs.append(("B7_OSTIR", "B6_TSS_accessibility"))
        for feature_set in FEATURE_SETS:
            xgb_name = feature_set
            ridge_name = f"Ridge_{feature_set}"
            if xgb_name in models and ridge_name in models:
                pairs.append((xgb_name, ridge_name))
        for pair in pairs:
            deltas[pair] = {"spearman": [], "mae_log": []}
        for _ in range(n_boot):
            sample = _resample_group_rows(frame, scenario, rng)
            if sample.empty:
                continue
            sample_metrics = {}
            for model in models:
                subset = sample.loc[sample.baseline == model]
                metric = _metrics(subset.actual.to_numpy(float), subset.predicted.to_numpy(float))
                sample_metrics[model] = metric
                for key in ("spearman", "mae_log"):
                    if np.isfinite(metric[key]):
                        sampled[model][key].append(metric[key])
            for left, right in pairs:
                if left in sample_metrics and right in sample_metrics:
                    for key in ("spearman", "mae_log"):
                        a, b = sample_metrics[left][key], sample_metrics[right][key]
                        if np.isfinite(a) and np.isfinite(b):
                            deltas[(left, right)][key].append(a - b)
        for model, metric_samples in sampled.items():
            for metric, values in metric_samples.items():
                output.append(_ci_record(scenario, target, model, metric, values, n_boot))
        for (left, right), metric_samples in deltas.items():
            for metric, values in metric_samples.items():
                output.append(_ci_record(scenario, target, f"DELTA:{left}-{right}", metric, values, n_boot))
    return pd.DataFrame(output)


def _ci_record(scenario, target, label, metric, values, requested):
    return {
        "scenario": scenario, "target": target, "baseline_or_delta": label,
        "metric": metric, "n_boot_requested": requested, "n_boot_valid": len(values),
        "ci_2_5": float(np.quantile(values, .025)) if values else float("nan"),
        "ci_97_5": float(np.quantile(values, .975)) if values else float("nan"),
        "bootstrap_median": float(np.median(values)) if values else float("nan"),
    }


def summarize_oof(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in predictions.groupby(["scenario", "fold", "target", "baseline"], sort=False):
        metrics = _metrics(group.actual.to_numpy(float), group.predicted.to_numpy(float))
        rows.append(dict(zip(["scenario", "fold", "target", "baseline"], keys),
                         n_test=len(group), **metrics))
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="reports/biological_context/validation_dataset.parquet")
    parser.add_argument("--ostir", default="reports/biological_context/ostir_predictions.csv")
    parser.add_argument("--output-dir", default="reports/biological_context")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--skip-bootstrap", action="store_true",
        help="Skip group bootstrap when only comparing outer split seeds.",
    )
    args = parser.parse_args()
    data_path = Path(args.data)
    required = set(BASIC_FEATURES + ACCESS_FEATURES + ["legacy_promoter_tail_mfe"])
    if data_path.exists():
        df = pd.read_parquet(data_path)
    else:
        df = pd.DataFrame()
    if not required.issubset(df.columns):
        source_path = "data/processed/constructs.parquet"
        df = prepare_benchmark_data(source_path).reset_index(drop=True)
        data_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(data_path, index=False)
    ostir = pd.read_csv(args.ostir)
    predictions, assignments, importance, ridge_coefficients, coverage = run_grouped_cv(
        df, ostir, args.folds, args.random_state
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output / "grouped_oof_predictions.csv", index=False)
    assignments.to_csv(output / "grouped_fold_assignments.csv", index=False)
    summarize_oof(predictions).to_csv(output / "grouped_fold_metrics.csv", index=False)
    importance.to_csv(output / "grouped_b6_importance.csv", index=False)
    ridge_coefficients.to_csv(output / "grouped_ridge_coefficients.csv", index=False)
    if args.skip_bootstrap:
        bootstrap_replicates = 0
    else:
        ci = run_group_bootstrap(predictions, args.bootstrap_replicates, args.random_state + 1)
        ci.to_csv(output / "group_bootstrap_ci.csv", index=False)
        bootstrap_replicates = args.bootstrap_replicates
    (output / "grouped_validation_metadata.json").write_text(
        json.dumps({"folds": args.folds, "random_state": args.random_state,
                    "bootstrap_replicates": bootstrap_replicates,
                    "bootstrap_skipped": bool(args.skip_bootstrap),
                    "bootstrap_seed": None if args.skip_bootstrap else args.random_state + 1,
                    **coverage}, indent=2) + "\n"
    )
    print(summarize_oof(predictions).to_string(index=False))
    if args.skip_bootstrap:
        print(f"\nWrote grouped CV reports (bootstrap skipped) to {output}")
    else:
        print(f"\nWrote grouped CV and bootstrap reports to {output}")


if __name__ == "__main__":
    main()
