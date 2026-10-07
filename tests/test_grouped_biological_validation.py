import numpy as np
import pandas as pd
import pytest

from scripts.grouped_biological_validation import (
    RIDGE_ALPHAS,
    _fit_ridge,
    make_inner_cv_splits,
    make_group_cv_splits,
    run_group_bootstrap,
)
from scripts.run_ostir_baseline import build_ostir_input, parse_ostir_output


def _factorial_data(n_groups=5):
    return pd.DataFrame([
        {"promoter_id": f"p{p}", "rbs_id": f"r{r}"}
        for p in range(n_groups) for r in range(n_groups)
    ])


def test_group_cv_splits_exclude_held_groups_and_save_unused_cross_pairs():
    df = _factorial_data(5)
    splits, assignments = make_group_cv_splits(df, n_splits=5, random_state=13)
    assert len(splits["promoter_holdout"]) == 5
    assert len(splits["rbs_holdout"]) == 5
    for train, test in splits["promoter_holdout"]:
        assert set(df.iloc[train].promoter_id).isdisjoint(df.iloc[test].promoter_id)
    for train, test in splits["rbs_holdout"]:
        assert set(df.iloc[train].rbs_id).isdisjoint(df.iloc[test].rbs_id)
    for train, test in splits["double_unseen"]:
        assert set(df.iloc[train].promoter_id).isdisjoint(df.iloc[test].promoter_id)
        assert set(df.iloc[train].rbs_id).isdisjoint(df.iloc[test].rbs_id)
        assert len(test) == 1
        assert len(train) == 16
    assert set(assignments.role) == {"train", "test", "unused"}
    assert (assignments.row_index == assignments.benchmark_row_id).all()


def test_group_folds_are_reproducible():
    df = _factorial_data(6)
    first, _ = make_group_cv_splits(df, n_splits=3, random_state=91)
    second, _ = make_group_cv_splits(df, n_splits=3, random_state=91)
    for scenario in first:
        for pair_a, pair_b in zip(first[scenario], second[scenario]):
            assert np.array_equal(pair_a[0], pair_b[0])
            assert np.array_equal(pair_a[1], pair_b[1])


@pytest.mark.parametrize("scenario", ["promoter_holdout", "rbs_holdout", "double_unseen"])
def test_ridge_inner_splits_follow_outer_grouping(scenario):
    df = _factorial_data(6)
    splits = make_inner_cv_splits(df, scenario, n_splits=3, random_state=11)
    for train, validation in splits:
        if scenario in {"promoter_holdout", "double_unseen"}:
            assert set(df.iloc[train].promoter_id).isdisjoint(df.iloc[validation].promoter_id)
        if scenario in {"rbs_holdout", "double_unseen"}:
            assert set(df.iloc[train].rbs_id).isdisjoint(df.iloc[validation].rbs_id)
        assert len(train) and len(validation)


def test_ridge_standardizes_features_and_selects_alpha_inside_train():
    df = pd.DataFrame([
        {"promoter_id": f"p{p}", "rbs_id": f"r{r}",
         "x": 1000.0 + p * 10 + r, "y": 2.5 * (p * 10 + r)}
        for p in range(6) for r in range(6)
    ])
    predicted, alpha, coefficients = _fit_ridge(
        df.iloc[:30], df.iloc[30:], ["x"], "y", "random", random_state=3
    )
    assert alpha in set(RIDGE_ALPHAS)
    assert np.isfinite(predicted).all()
    assert np.isfinite(coefficients).all()


def test_ostir_input_uses_exact_downstream_reporter_aug_position():
    context = "GAAUGAUCUUAAUCUAGCCCAGGAACGUUUCAUAUG" + "CGTAAAGGCGAAGAGCTGTTCACTGGTTTCGTCACTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGT"
    result = build_ostir_input(pd.DataFrame({"transcript_context": [context]}))
    row = result.iloc[0]
    assert row["name"] == "0"
    assert int(row["start"]) == len(context) - 89
    assert row["seq"][int(row["start"]) - 1:int(row["start"]) + 2] == "AUG"
    assert row["start"] == row["end"]


def test_ostir_output_requires_one_positive_prediction_for_every_row():
    valid = pd.DataFrame({
        "name": ["0", "1"], "start_position": [25, 27], "expression": [12.5, 2.1],
    })
    parsed = parse_ostir_output(valid, n_rows=2)
    assert parsed.row_index.tolist() == [0, 1]
    assert parsed.ostir_expression.tolist() == [12.5, 2.1]
    with pytest.raises(ValueError, match="incomplete"):
        parse_ostir_output(valid.iloc[:1], n_rows=2)
    partial = parse_ostir_output(valid.iloc[:1], n_rows=2, allow_missing=True)
    assert partial.row_index.tolist() == [0]
    with pytest.raises(ValueError, match="non-positive"):
        parse_ostir_output(pd.DataFrame({
            "name": ["0"], "start_position": [25], "expression": [0.0],
        }), n_rows=1)


def test_group_bootstrap_returns_paired_delta_confidence_intervals():
    rows = []
    for baseline, shift in (("B4_sequence_only", 0.2), ("B6_TSS_accessibility", 0.0), ("B7_OSTIR", 0.3)):
        for p in range(5):
            for r in range(3):
                actual = float(p + r)
                rows.append({
                    "scenario": "promoter_holdout", "target": "log_translation_proxy",
                    "baseline": baseline, "promoter_id": f"p{p}", "rbs_id": f"r{r}",
                    "actual": actual, "predicted": actual + shift,
                })
    ci = run_group_bootstrap(pd.DataFrame(rows), n_boot=30, random_state=4)
    delta = ci.loc[ci.baseline_or_delta == "DELTA:B6_TSS_accessibility-B4_sequence_only"]
    assert set(delta.metric) == {"spearman", "mae_log"}
    assert (delta.n_boot_valid == 30).all()
    assert np.isfinite(delta[["ci_2_5", "ci_97_5"]].to_numpy()).all()
