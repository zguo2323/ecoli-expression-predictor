"""Summarize grouped-CV stability across outer split seeds."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

METRICS_FILE = "grouped_fold_metrics.csv"
COMPARISONS = (
    ("B6_TSS_accessibility", "B4_sequence_only"),
    ("B6_TSS_accessibility", "Ridge_B6_TSS_accessibility"),
    ("B6_TSS_accessibility", "B7_OSTIR"),
)


def load_seed_runs(run_specs: list[str]) -> pd.DataFrame:
    """Load one fold-metric CSV for every seed, from `seed=directory` specs."""
    frames = []
    seen_seeds = set()
    for spec in run_specs:
        if "=" not in spec:
            raise ValueError(f"Expected seed=directory, received: {spec}")
        seed_text, directory = spec.split("=", 1)
        seed = int(seed_text)
        if seed in seen_seeds:
            raise ValueError(f"Duplicate seed: {seed}")
        seen_seeds.add(seed)
        path = Path(directory) / METRICS_FILE
        if not path.is_file():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = {"scenario", "fold", "target", "baseline", "spearman", "mae_log"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        frame["seed"] = seed
        frames.append(frame)
    if not frames:
        raise ValueError("Provide at least one --run seed=directory")
    return pd.concat(frames, ignore_index=True)


def summarize_multiseed(fold_metrics: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Report per-seed CV averages and paired model-delta stability."""
    keys = ["seed", "scenario", "target", "baseline"]
    seed_level = (fold_metrics.groupby(keys, sort=True)
                  .agg(n_folds=("fold", "nunique"),
                       mean_fold_spearman=("spearman", "mean"),
                       sd_fold_spearman=("spearman", "std"),
                       mean_fold_mae_log=("mae_log", "mean"),
                       sd_fold_mae_log=("mae_log", "std"))
                  .reset_index())
    cross_seed = (seed_level.groupby(["scenario", "target", "baseline"], sort=True)
                  .agg(n_seeds=("seed", "nunique"),
                       mean_seed_spearman=("mean_fold_spearman", "mean"),
                       sd_across_seeds_spearman=("mean_fold_spearman", "std"),
                       min_seed_spearman=("mean_fold_spearman", "min"),
                       max_seed_spearman=("mean_fold_spearman", "max"),
                       mean_seed_mae_log=("mean_fold_mae_log", "mean"),
                       sd_across_seeds_mae_log=("mean_fold_mae_log", "std"),
                       min_seed_mae_log=("mean_fold_mae_log", "min"),
                       max_seed_mae_log=("mean_fold_mae_log", "max"))
                  .reset_index())

    index = ["seed", "scenario", "fold", "target"]
    wide = fold_metrics.pivot(index=index, columns="baseline", values=["spearman", "mae_log"])
    delta_rows = []
    for better, comparator in COMPARISONS:
        for metric in ("spearman", "mae_log"):
            if (metric, better) not in wide or (metric, comparator) not in wide:
                continue
            delta = (wide[(metric, better)] - wide[(metric, comparator)]).dropna()
            for (seed, scenario, fold, target), value in delta.items():
                delta_rows.append({"seed": seed, "scenario": scenario, "fold": fold,
                                   "target": target, "better_model": better,
                                   "comparator_model": comparator, "metric": metric,
                                   "delta_better_minus_comparator": float(value)})
    deltas = pd.DataFrame(delta_rows)
    delta_seed = (deltas.groupby(["seed", "scenario", "target", "better_model",
                                  "comparator_model", "metric"], sort=True)
                  .agg(n_paired_folds=("fold", "nunique"),
                       mean_fold_delta=("delta_better_minus_comparator", "mean"),
                       sd_fold_delta=("delta_better_minus_comparator", "std"))
                  .reset_index())
    delta_summary_rows = []
    summary_keys = ["scenario", "target", "better_model", "comparator_model", "metric"]
    for group_key, group in delta_seed.groupby(summary_keys, sort=True):
        values = group.mean_fold_delta
        metric = group_key[-1]
        favors_better = values > 0 if metric == "spearman" else values < 0
        delta_summary_rows.append(dict(zip(summary_keys, group_key),
                                       n_seeds=group.seed.nunique(),
                                       mean_seed_delta=float(values.mean()),
                                       sd_across_seeds=float(values.std()),
                                       min_seed_delta=float(values.min()),
                                       max_seed_delta=float(values.max()),
                                       seeds_favoring_better_model=int(favors_better.sum())))
    delta_summary = pd.DataFrame(delta_summary_rows)
    return {
        "fold_metrics": fold_metrics,
        "seed_level_metrics": seed_level,
        "cross_seed_summary": cross_seed,
        "paired_deltas_by_fold": deltas,
        "paired_deltas_by_seed": delta_seed,
        "paired_delta_summary": delta_summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True,
                        help="A seed=directory pair; repeat for each seed.")
    parser.add_argument("--output-dir", default="reports/biological_context/multiseed")
    args = parser.parse_args()
    reports = summarize_multiseed(load_seed_runs(args.run))
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for key, frame in reports.items():
        path = output / f"multiseed_{key}.csv"
        frame.to_csv(path, index=False)
        print(f"Wrote {path}: {len(frame)} rows")
    print("\nTranslation-proxy stability (mean ± SD across split seeds):")
    display = reports["cross_seed_summary"].query(
        "target == 'log_translation_proxy' and scenario in "
        "['promoter_holdout', 'rbs_holdout', 'double_unseen'] and "
        "baseline in ['B6_TSS_accessibility', 'Ridge_B6_TSS_accessibility']"
    )
    columns = ["scenario", "baseline", "n_seeds", "mean_seed_spearman",
               "sd_across_seeds_spearman", "mean_seed_mae_log",
               "sd_across_seeds_mae_log"]
    print(display[columns].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
