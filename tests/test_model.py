import sys
import os
import pytest
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.model.model import (
    train_model,
    predict_expression,
    rank_rbs_for_promoter,
    evaluate_model,
)

PARQUET = "data/processed/constructs.parquet"
MODEL = "artifacts/model.pkl"


@pytest.fixture(scope="session", autouse=True)
def trained_model(tmp_path_factory):
    """Train on a fixed-size sample so integration tests stay fast and isolated."""
    global PARQUET, MODEL
    tmp_path = tmp_path_factory.mktemp("model-integration")
    source = pd.read_parquet(PARQUET)
    sample = source.sample(n=600, random_state=42)
    PARQUET = str(tmp_path / "constructs.parquet")
    MODEL = str(tmp_path / "model.pkl")
    sample.to_parquet(PARQUET, index=False)
    train_model(PARQUET, model_output_path=MODEL)


def test_predict_expression_returns_required_keys():
    result = predict_expression("TTGACATATAATCCGG", "AAAGAGGAGAAA", MODEL)
    assert "predicted_prot" in result
    assert "features_used" in result
    assert "translation_context_available" in result
    assert "confidence_interval" not in result


def test_predict_expression_positive_output():
    result = predict_expression("TTGACATATAATCCGG", "AAAGAGGAGAAA", MODEL)
    assert result["predicted_prot"] > 0


def test_rank_rbs_returns_top_n():
    results = rank_rbs_for_promoter("TTGACATATAATCCGG", model_path=MODEL, top_n=5)
    assert len(results) == 5


def test_rank_rbs_sorted_descending():
    results = rank_rbs_for_promoter("TTGACATATAATCCGG", model_path=MODEL, top_n=10)
    prot_values = [r["predicted_prot"] for r in results]
    assert prot_values == sorted(prot_values, reverse=True)


def test_rank_rbs_context_dependent():
    # Different promoters should produce different rankings
    results1 = rank_rbs_for_promoter("TTGACATATAATCCGG", model_path=MODEL, top_n=5)
    results2 = rank_rbs_for_promoter("GCGCGCGCGCGCGCGC", model_path=MODEL, top_n=5)
    top_ids1 = [r["rbs_id"] for r in results1]
    top_ids2 = [r["rbs_id"] for r in results2]
    assert top_ids1 != top_ids2


def test_evaluate_model_metrics_are_finite():
    result = evaluate_model(PARQUET, model_path=MODEL)
    assert np.isfinite(result["spearman_r"])
    assert np.isfinite(result["mae"])


def test_evaluate_model_returns_required_keys():
    result = evaluate_model(PARQUET, model_path=MODEL)
    for key in ["spearman_r", "r2", "mae", "n_test"]:
        assert key in result
