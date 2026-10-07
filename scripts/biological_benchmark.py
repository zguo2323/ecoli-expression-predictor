"""Reproducible baselines for biological-context validation.

Run after rebuilding the processed dataset with the current pipeline:
  python scripts/biological_benchmark.py --data data/processed/constructs.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from modules.features.features import compute_mrna_folding_energy
from modules.model.model import FEATURE_COLS, load_and_featurize
from modules.pipeline.pipeline import load_promoter_table, load_rbs_table

BASIC_FEATURES = [
    "gc_promoter", "gc_rbs", "score_minus10", "score_minus35",
    "spacer_length", "spacer_optimal", "score_sd", "sd_spacing", "sd_spacing_optimal",
]
ACCESS_FEATURES = [
    "sd_unpaired_probability", "start_unpaired_probability", "sd_start_opening_energy",
]
TARGETS = {"log_prot": "prot", "log_translation_proxy": None, "log_RNA": "RNA"}


def make_group_splits(df: pd.DataFrame, test_size: float = 0.2, random_state: int = 42):
    """Return row indices for row-random and strict promoter/RBS holdouts."""
    all_idx = np.arange(len(df))
    random_train, random_test = train_test_split(
        all_idx, test_size=test_size, random_state=random_state
    )

    def heldout_ids(column):
        ids = np.array(sorted(df[column].dropna().unique()))
        train_ids, test_ids = train_test_split(ids, test_size=test_size, random_state=random_state)
        return set(train_ids), set(test_ids)

    train_p, test_p = heldout_ids("promoter_id")
    train_r, test_r = heldout_ids("rbs_id")
    p = df["promoter_id"].to_numpy()
    r = df["rbs_id"].to_numpy()
    return {
        "random": (random_train, random_test),
        "promoter_holdout": (all_idx[np.isin(p, list(train_p))], all_idx[np.isin(p, list(test_p))]),
        "rbs_holdout": (all_idx[np.isin(r, list(train_r))], all_idx[np.isin(r, list(test_r))]),
        "double_unseen": (
            all_idx[np.isin(p, list(train_p)) & np.isin(r, list(train_r))],
            all_idx[np.isin(p, list(test_p)) & np.isin(r, list(test_r))],
        ),
    }


def prepare_benchmark_data(path: str) -> pd.DataFrame:
    df = load_and_featurize(path)
    df["benchmark_row_id"] = np.arange(len(df), dtype=int)
    promoters = load_promoter_table("data/raw/sd01.xls")[["promoter_id", "mean_RNA"]]
    rbs = load_rbs_table("data/raw/sd02.xls")[["rbs_id", "mean_xlat"]]
    df = df.merge(promoters, on="promoter_id", how="left").merge(rbs, on="rbs_id", how="left")
    df = df.loc[(df["prot"] > 0) & (df["RNA"] > 0)].copy()
    df["log_prot"] = np.log(df["prot"].astype(float))
    df["log_RNA"] = np.log(df["RNA"].astype(float))
    df["log_translation_proxy"] = df["log_prot"] - df["log_RNA"]
    # Historical feature is retained only as an explicit ablation baseline.
    df["legacy_promoter_tail_mfe"] = [
        compute_mrna_folding_energy(p, r) for p, r in zip(df["promo_seq"], df["rbs_seq"])
    ]
    return df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=FEATURE_COLS + ["log_prot", "log_RNA", "log_translation_proxy", "mean_RNA", "mean_xlat"]
    )


def _fit_xgb(train: pd.DataFrame, test: pd.DataFrame, features: list[str], target: str):
    model = XGBRegressor(
        n_estimators=300, max_depth=5, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, random_state=42,
        eval_metric="rmse", n_jobs=1,
    )
    model.fit(train[features].astype(float), train[target].astype(float))
    return model.predict(test[features].astype(float)), model


def _fit_additive(train: pd.DataFrame, test: pd.DataFrame, target: str):
    pre = ColumnTransformer([
        ("ids", OneHotEncoder(handle_unknown="ignore"), ["promoter_id", "rbs_id"]),
    ])
    model = make_pipeline(pre, Ridge(alpha=10.0))
    model.fit(train[["promoter_id", "rbs_id"]], train[target])
    return model.predict(test[["promoter_id", "rbs_id"]]), model


def _metrics(actual, predicted):
    rank_defined = np.ptp(actual) > 0 and np.ptp(predicted) > 0
    return {
        "spearman": float(spearmanr(actual, predicted).statistic) if rank_defined else float("nan"),
        "mae_log": float(mean_absolute_error(actual, predicted)),
        "r2_log": float(r2_score(actual, predicted)),
    }


def run_benchmark(df: pd.DataFrame, test_size: float = 0.2, random_state: int = 42):
    results, predictions, importances = [], [], []
    splits = make_group_splits(df, test_size, random_state)
    for split_name, (train_idx, test_idx) in splits.items():
        train, test = df.iloc[train_idx], df.iloc[test_idx]
        if not len(train) or not len(test):
            continue
        for target in TARGETS:
            actual = test[target].to_numpy(dtype=float)
            variants = {
                "B0_train_median": np.full(len(test), train[target].median()),
            }
            if target == "log_prot":
                oracle = np.log(test["mean_RNA"].to_numpy()) + np.log(test["mean_xlat"].to_numpy())
            elif target == "log_translation_proxy":
                oracle = np.log(test["mean_xlat"].to_numpy())
            else:
                oracle = np.log(test["mean_RNA"].to_numpy())
            variants["B2_measured_component_oracle"] = oracle
            additive, _ = _fit_additive(train, test, target)
            variants["B1_additive_promoter_rbs"] = additive
            feature_sets = {
                "B4_sequence_only": BASIC_FEATURES,
                "B3_historical_10_feature_model": BASIC_FEATURES + ["legacy_promoter_tail_mfe"],
                "B6_TSS_accessibility": BASIC_FEATURES + ACCESS_FEATURES,
            }
            fitted = {}
            for name, features in feature_sets.items():
                predictions_for_model, model = _fit_xgb(train, test, features, target)
                variants[name] = predictions_for_model
                fitted[name] = model
            for variant, predicted in variants.items():
                metrics = _metrics(actual, predicted)
                results.append({
                    "split": split_name, "target": target, "baseline": variant,
                    "n_train": len(train), "n_test": len(test), **metrics,
                })
                predictions.extend({
                    "split": split_name, "target": target, "baseline": variant,
                    "row_index": int(i), "promoter_id": test.iloc[j]["promoter_id"],
                    "rbs_id": test.iloc[j]["rbs_id"], "actual": actual[j],
                    "predicted": float(predicted[j]), "residual": float(actual[j] - predicted[j]),
                } for j, i in enumerate(test_idx))
            model = fitted["B6_TSS_accessibility"]
            importances.extend({
                "split": split_name, "target": target,
                "feature": name, "importance": float(value),
            } for name, value in zip(BASIC_FEATURES + ACCESS_FEATURES, model.feature_importances_))
    return pd.DataFrame(results), pd.DataFrame(predictions), pd.DataFrame(importances)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/processed/constructs.parquet")
    parser.add_argument("--output-dir", default="reports/biological_context")
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    df = prepare_benchmark_data(args.data)
    metrics, predictions, importance = run_benchmark(df, args.test_size, args.random_state)
    metrics.to_csv(output / "baseline_metrics.csv", index=False)
    predictions.to_csv(output / "heldout_predictions.csv", index=False)
    importance.to_csv(output / "b6_feature_importance.csv", index=False)
    print(metrics.to_string(index=False))
    print(f"\nWrote reports to {output}")


if __name__ == "__main__":
    main()
