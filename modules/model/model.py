import os
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
from modules.features.features import extract_all_features
from modules.pipeline.pipeline import load_rbs_table

FEATURE_COLS = [
    "gc_promoter", "gc_rbs",
    "score_minus10", "score_minus35",
    "spacer_length", "spacer_optimal",
    "score_sd", "sd_spacing", "sd_spacing_optimal",
    "mrna_folding_energy",
]


def load_and_featurize(parquet_path: str) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)
    features = df.apply(
        lambda row: extract_all_features(row["promo_seq"], row["rbs_seq"]),
        axis=1,
        result_type="expand",
    )
    return pd.concat([df.reset_index(drop=True), features], axis=1)


def train_model(
    parquet_path: str,
    target: str = "prot",
    model_output_path: str = "artifacts/model.pkl",
    test_size: float = 0.1,
    val_size: float = 0.1,
    random_state: int = 42,
) -> dict:
    print(f"[1/5] Featurizing dataset from {parquet_path} ...")
    df = load_and_featurize(parquet_path)
    df = df.dropna(subset=FEATURE_COLS + [target])
    print(f"      {len(df)} usable rows, {len(FEATURE_COLS)} features")

    # Stratify by expression bin so all splits see the full expression range
    df["_bin"] = pd.qcut(df[target], q=5, labels=False, duplicates="drop")

    X = df[FEATURE_COLS].astype(float)
    y = df[target].astype(float)
    strat = df["_bin"]

    X_train_val, X_test, y_train_val, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=strat
    )
    strat_train_val = strat.loc[X_train_val.index]
    relative_val = val_size / (1 - test_size)
    X_train, X_val, y_train, y_val = train_test_split(
        X_train_val, y_train_val,
        test_size=relative_val,
        random_state=random_state,
        stratify=strat_train_val,
    )
    print(f"[2/5] Split: {len(X_train)} train / {len(X_val)} val / {len(X_test)} test")

    print(f"[3/5] Training XGBoost (n_estimators=300) ...")
    model = XGBRegressor(
        n_estimators=300,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=random_state,
        eval_metric="rmse",
    )
    model.fit(
        X_train, y_train,
        eval_set=[(X_val, y_val)],
        verbose=False,
    )

    print(f"[4/5] Evaluating on test set ...")
    preds = model.predict(X_test)
    spearman_r = float(spearmanr(y_test, preds).statistic)
    eval_dict = {
        "spearman_r": spearman_r,
        "r2": float(r2_score(y_test, preds)),
        "mae": float(mean_absolute_error(y_test, preds)),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "target": target,
    }

    print(f"[5/5] Saving model to {model_output_path} ...")
    os.makedirs(os.path.dirname(model_output_path), exist_ok=True)
    with open(model_output_path, "wb") as f:
        pickle.dump({
            "model": model,
            "feature_names": FEATURE_COLS,
            "eval": eval_dict,
        }, f)

    print(f"\n=== Training Complete ===")
    print(f"  Spearman r : {eval_dict['spearman_r']:.4f}")
    print(f"  R²         : {eval_dict['r2']:.4f}")
    print(f"  MAE        : {eval_dict['mae']:.2f}")
    print(f"  Train rows : {eval_dict['n_train']}")
    print(f"  Test rows  : {eval_dict['n_test']}")
    return eval_dict


def _load_model(model_path: str) -> dict:
    with open(model_path, "rb") as f:
        return pickle.load(f)


def predict_expression(
    promoter_seq: str,
    rbs_seq: str,
    model_path: str = "artifacts/model.pkl",
) -> dict:
    features = extract_all_features(promoter_seq, rbs_seq)
    artifact = _load_model(model_path)
    model = artifact["model"]
    feature_names = artifact["feature_names"]

    X = pd.DataFrame([{k: float(v) for k, v in features.items()}])[feature_names]
    pred = float(model.predict(X)[0])

    return {
        "predicted_prot": pred,
        "confidence_interval": [pred * 0.75, pred * 1.25],
        "features_used": features,
    }


def rank_rbs_for_promoter(
    promoter_seq: str,
    rbs_library_path: str = "data/raw/sd02.xls",
    model_path: str = "artifacts/model.pkl",
    top_n: int = 5,
) -> list:
    rbs_df = load_rbs_table(rbs_library_path)
    results = []
    for _, row in rbs_df.iterrows():
        pred = predict_expression(promoter_seq, row["sequence"], model_path)
        results.append({
            "rbs_id": row["rbs_id"],
            "rbs_seq": row["sequence"],
            "predicted_prot": pred["predicted_prot"],
        })
    results.sort(key=lambda x: x["predicted_prot"], reverse=True)
    return results[:top_n]


def evaluate_model(
    parquet_path: str,
    model_path: str = "artifacts/model.pkl",
    target: str = "prot",
    test_size: float = 0.1,
    random_state: int = 42,
) -> dict:
    df = load_and_featurize(parquet_path)
    df = df.dropna(subset=FEATURE_COLS + [target])
    df["_bin"] = pd.qcut(df[target], q=5, labels=False, duplicates="drop")

    X = df[FEATURE_COLS].astype(float)
    y = df[target].astype(float)
    strat = df["_bin"]

    _, X_test, _, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=strat
    )

    artifact = _load_model(model_path)
    model = artifact["model"]
    preds = model.predict(X_test)

    return {
        "spearman_r": float(spearmanr(y_test, preds).statistic),
        "r2": float(r2_score(y_test, preds)),
        "mae": float(mean_absolute_error(y_test, preds)),
        "n_test": len(X_test),
    }


class ModelMCP:
    def __init__(self, config):
        self.config = config

    def initiate(self):
        pass

    def run(self, **kwargs):
        mcp_name = self.config.get("execution_details", {}).get("mcp_name")
        if mcp_name == "predict_expression":
            return predict_expression(**kwargs)
        if mcp_name == "rank_rbs_for_promoter":
            return rank_rbs_for_promoter(**kwargs)
        if mcp_name == "evaluate_model":
            return evaluate_model(**kwargs)
        raise ValueError(f"Unknown mcp_name: {mcp_name}")
