import numpy as np
import pandas as pd
import pytest

from scripts.grouped_feature_group_permutation import (
    FEATURE_GROUPS,
    FEATURE_SET,
    permutation_metric_deltas,
    permute_feature_group,
    summarize_permutation_importance,
    validate_feature_groups,
)


def test_feature_groups_partition_all_b6_features():
    validate_feature_groups()
    with pytest.raises(ValueError, match="partition"):
        validate_feature_groups({"bad_group": FEATURE_SET[:-1]}, FEATURE_SET)


def test_feature_block_permutation_keeps_within_block_rows_together():
    unit_values = np.repeat(np.arange(6, dtype=float), 2)
    frame = pd.DataFrame({
        "unit": np.repeat([f"u{i}" for i in range(6)], 2),
        "feature_a": unit_values,
        "feature_b": unit_values * 10,
        "untouched": np.arange(12, dtype=float) + 100,
    })
    permuted = permute_feature_group(
        frame, ["feature_a", "feature_b"], np.random.default_rng(7), ["unit"]
    )
    assert np.array_equal(permuted.feature_b, permuted.feature_a * 10)
    assert np.array_equal(permuted.untouched, frame.untouched)
    assert (permuted.groupby("unit").feature_a.nunique() == 1).all()
    assert not np.array_equal(permuted.drop_duplicates("unit").feature_a.to_numpy(),
                              frame.drop_duplicates("unit").feature_a.to_numpy())


def test_feature_group_permutation_rejects_values_that_vary_within_unit():
    frame = pd.DataFrame({"unit": ["u1", "u1"], "x": [1.0, 2.0]})
    with pytest.raises(ValueError, match="constant within permutation units"):
        permute_feature_group(frame, ["x"], np.random.default_rng(3), ["unit"])


def test_positive_deltas_mean_permuting_group_hurts_prediction():
    actual = np.arange(20, dtype=float)
    baseline = actual.copy()
    permuted = actual[::-1].copy()
    delta = permutation_metric_deltas(actual, baseline, permuted)
    assert delta["delta_spearman"] > 0
    assert delta["delta_mae_log"] > 0


def test_summary_reports_repeat_mean_and_equal_fold_mean_sd():
    rows = []
    for fold, values in ((0, [0.1, 0.3]), (1, [0.5, 0.7])):
        for repeat, value in enumerate(values):
            rows.append({
                "scenario": "promoter_holdout", "fold": fold,
                "target": "log_translation_proxy",
                "feature_group": next(iter(FEATURE_GROUPS)),
                "delta_spearman": value, "delta_mae_log": value * 2,
                "permutation": repeat, "n_test": 10,
            })
    by_fold, summary = summarize_permutation_importance(pd.DataFrame(rows))
    result = summary.iloc[0]
    assert by_fold.delta_spearman.tolist() == [0.2, 0.6]
    assert result.mean_delta_spearman == pytest.approx(0.4)
    assert result.sd_across_folds_spearman == pytest.approx(np.std([0.2, 0.6], ddof=1))
    assert result.mean_delta_mae_log == pytest.approx(0.8)
