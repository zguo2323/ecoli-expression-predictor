import sys
import os
import pytest

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
def trained_model():
    """Train once for the whole test session."""
    train_model(PARQUET, model_output_path=MODEL)


def test_predict_expression_returns_required_keys():
    result = predict_expression("TTGACATATAATCCGG", "AAAGAGGAGAAA", MODEL)
    assert "predicted_prot" in result
    assert "confidence_interval" in result
    assert "features_used" in result


def test_predict_expression_positive_output():
    result = predict_expression("TTGACATATAATCCGG", "AAAGAGGAGAAA", MODEL)
    assert result["predicted_prot"] > 0


def test_confidence_interval_ordered():
    result = predict_expression("TTGACATATAATCCGG", "AAAGAGGAGAAA", MODEL)
    assert result["confidence_interval"][0] < result["confidence_interval"][1]


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


def test_evaluate_model_spearman_reasonable():
    result = evaluate_model(PARQUET, model_path=MODEL)
    assert result["spearman_r"] > 0.3


def test_evaluate_model_returns_required_keys():
    result = evaluate_model(PARQUET, model_path=MODEL)
    for key in ["spearman_r", "r2", "mae", "n_test"]:
        assert key in result
