import pandas as pd

from scripts.summarize_multiseed_grouped_validation import summarize_multiseed


def _fold_metrics():
    rows = []
    for seed, scenario, fold, values in [
        (1, "rbs_holdout", 0, {"B6_TSS_accessibility": (0.8, 0.5),
                              "B4_sequence_only": (0.7, 0.6),
                              "Ridge_B6_TSS_accessibility": (0.6, 0.7)}),
        (1, "rbs_holdout", 1, {"B6_TSS_accessibility": (0.7, 0.6),
                              "B4_sequence_only": (0.6, 0.7),
                              "Ridge_B6_TSS_accessibility": (0.5, 0.8)}),
        (2, "rbs_holdout", 0, {"B6_TSS_accessibility": (0.6, 0.7),
                              "B4_sequence_only": (0.7, 0.6),
                              "Ridge_B6_TSS_accessibility": (0.7, 0.6)}),
        (2, "rbs_holdout", 1, {"B6_TSS_accessibility": (0.5, 0.8),
                              "B4_sequence_only": (0.6, 0.7),
                              "Ridge_B6_TSS_accessibility": (0.6, 0.7)}),
    ]:
        for baseline, (spearman, mae) in values.items():
            rows.append({"seed": seed, "scenario": scenario, "fold": fold,
                         "target": "log_translation_proxy", "baseline": baseline,
                         "spearman": spearman, "mae_log": mae})
    return pd.DataFrame(rows)


def test_multiseed_summary_reports_seed_variation_and_paired_directions():
    reports = summarize_multiseed(_fold_metrics())
    summary = reports["cross_seed_summary"]
    b6 = summary.loc[summary.baseline == "B6_TSS_accessibility"].iloc[0]
    assert b6.n_seeds == 2
    assert b6.mean_seed_spearman == 0.65
    assert b6.mean_seed_mae_log == 0.65

    paired = reports["paired_delta_summary"]
    b6_vs_ridge = paired.loc[
        (paired.better_model == "B6_TSS_accessibility")
        & (paired.comparator_model == "Ridge_B6_TSS_accessibility")
    ].set_index("metric")
    assert b6_vs_ridge.loc["spearman", "seeds_favoring_better_model"] == 1
    assert b6_vs_ridge.loc["mae_log", "seeds_favoring_better_model"] == 1
