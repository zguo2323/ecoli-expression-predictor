import numpy as np
import pandas as pd
import pytest

from scripts.grouped_oof_error_analysis import (
    add_actual_quantile_bins,
    component_outcome_variance,
    component_error_table,
    paired_component_comparison,
    run_analysis,
    select_extreme_cases,
    summarize_paired_components,
    summarize_actual_bins,
)


def _predictions():
    rows = []
    for index in range(8):
        actual = float(index)
        for baseline, offset in (
            ("B6_TSS_accessibility", 0.1),
            ("Ridge_B6_TSS_accessibility", 0.5),
            ("B7_OSTIR", 1.0),
        ):
            predicted = actual - offset
            rows.append({
                "scenario": "double_unseen",
                "fold": index % 2,
                "target": "log_translation_proxy",
                "baseline": baseline,
                "benchmark_row_id": index,
                "promoter_id": f"p{index // 2}",
                "rbs_id": f"r{index % 4}",
                "actual": actual,
                "predicted": predicted,
                "residual": actual - predicted,
            })
    return pd.DataFrame(rows)


def test_actual_bins_are_shared_across_models_and_have_counts():
    predictions = _predictions()
    binned = add_actual_quantile_bins(predictions, n_bins=4)
    assert binned.groupby("benchmark_row_id").actual_bin.nunique().max() == 1
    summary = summarize_actual_bins(predictions, n_bins=4)
    assert summary.n.sum() == len(predictions)
    assert set(summary.actual_bin) == {1, 2, 3, 4}


def test_component_error_table_includes_double_unseen_pairs():
    components = component_error_table(_predictions())
    pairs = components.loc[components.component_type == "promoter_rbs_pair"]
    assert len(pairs) == 8 * 3
    assert pairs.n.eq(1).all()
    assert set(components.component_type) == {"promoter", "rbs", "promoter_rbs_pair"}


def test_marginal_outcome_variance_is_computed_by_component_id():
    predictions = pd.DataFrame({
        "scenario": ["rbs_holdout"] * 4,
        "target": ["log_translation_proxy"] * 4,
        "baseline": ["B6_TSS_accessibility"] * 4,
        "benchmark_row_id": range(4),
        "actual": [0.0, 0.0, 1.0, 1.0],
        "promoter_id": ["p1", "p1", "p2", "p2"],
        "rbs_id": ["r1", "r2", "r1", "r2"],
    })
    result = component_outcome_variance(predictions).set_index("component_type")
    assert result.loc["promoter", "marginal_eta_squared"] == pytest.approx(1.0)
    assert result.loc["rbs", "marginal_eta_squared"] == pytest.approx(0.0)


def test_paired_component_comparison_uses_same_rows_and_positive_favors_b6():
    comparison = paired_component_comparison(_predictions())
    rbs = comparison.loc[
        (comparison.component_type == "rbs")
        & (comparison.left_model == "B6_TSS_accessibility")
        & (comparison.right_model == "Ridge_B6_TSS_accessibility")
    ]
    assert len(rbs) == 4
    assert np.allclose(rbs.right_minus_left_mae_log, 0.4)
    assert rbs.fraction_left_better.eq(1.0).all()
    summary = summarize_paired_components(rbs)
    assert summary.iloc[0].fraction_components_left_better == 1.0


def test_extreme_cases_are_limited_per_scenario_target_model():
    cases = select_extreme_cases(_predictions(), n_cases=2)
    assert cases.groupby(["scenario", "target", "baseline"]).size().eq(2).all()
    assert set(cases.error_direction) == {"underprediction"}
    with pytest.raises(ValueError, match="positive"):
        select_extreme_cases(_predictions(), n_cases=0)


def test_run_analysis_returns_all_error_diagnostics():
    result = run_analysis(_predictions())
    assert set(result) == {
        "actual_bins", "component_errors", "component_summary",
        "target_component_variance", "paired_component_comparison",
        "paired_component_summary", "extreme_cases",
    }
    assert not result["component_summary"].empty
    assert set(result["target_component_variance"].component_type) == {"promoter", "rbs"}
