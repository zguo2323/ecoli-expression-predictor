"""Run OSTIR once for every complete TSS-to-CDS context in the benchmark data.

Example:
  python scripts/run_ostir_baseline.py --conda-env ostir-b7 --threads 4
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

def build_ostir_input(df: pd.DataFrame) -> pd.DataFrame:
    """Create one uniquely named OSTIR row with the exact reporter AUG position."""
    if "transcript_context" not in df:
        raise ValueError("Dataset must include transcript_context")
    rows = []
    row_ids = df["benchmark_row_id"] if "benchmark_row_id" in df else range(len(df))
    for row_index, sequence in zip(row_ids, df["transcript_context"].astype(str)):
        sequence = sequence.upper().replace("T", "U")
        start_1based = len(sequence) - 89  # 90 nt reporter CDS, including AUG
        if start_1based < 1 or sequence[start_1based - 1:start_1based + 2] != "AUG":
            raise ValueError(f"Construct row {row_index} has no reporter AUG at expected position")
        if any(base not in "ACGU" for base in sequence):
            raise ValueError(f"Construct row {row_index} contains ambiguous/non-RNA bases")
        rows.append({"name": str(row_index), "seq": sequence,
                     "start": start_1based, "end": start_1based})
    return pd.DataFrame(rows)


def parse_ostir_output(raw: pd.DataFrame, n_rows: int | list[int],
                       allow_missing: bool = False) -> pd.DataFrame:
    """Validate output; OSTIR may omit sequences without a scorable binding site."""
    if raw.empty:
        raise ValueError("OSTIR returned no binding/initiation predictions")
    required = {"name", "start_position", "expression"}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"OSTIR output is missing columns: {sorted(missing)}")
    out = raw.rename(columns={"expression": "ostir_expression"}).copy()
    out["row_index"] = pd.to_numeric(out["name"], errors="raise").astype(int)
    if out["row_index"].duplicated().any():
        raise ValueError("OSTIR returned multiple start sites for a construct")
    expected_rows = set(range(n_rows)) if isinstance(n_rows, int) else set(map(int, n_rows))
    observed_rows = set(out["row_index"])
    unexpected_rows = observed_rows.difference(expected_rows)
    if unexpected_rows:
        raise ValueError(f"OSTIR returned unexpected row indices: {sorted(unexpected_rows)[:10]}")
    missing_rows = expected_rows.difference(observed_rows)
    if missing_rows and not allow_missing:
        raise ValueError(f"OSTIR predictions incomplete; missing row indices include {sorted(missing_rows)[:10]}")
    if not np.isfinite(out["ostir_expression"].astype(float)).all() or (out["ostir_expression"] <= 0).any():
        raise ValueError("OSTIR returned non-positive or non-finite initiation rates")
    return out.sort_values("row_index").reset_index(drop=True)


def run_ostir(input_csv: Path, output_csv: Path, conda_env: str, threads: int) -> None:
    command = ["conda", "run", "--no-capture-output", "-n", conda_env,
               "ostir", "-i", str(input_csv), "-o", str(output_csv),
               "-t", "csv", "-j", str(threads), "-v", "0"]
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/processed/constructs.parquet")
    parser.add_argument("--output-dir", default="reports/biological_context")
    parser.add_argument("--conda-env", default="ostir-b7")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--limit", type=int, help="Run a prefix only for a smoke test")
    args = parser.parse_args()
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--threads and --limit must be positive")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Preparing TSS-to-CDS contexts from the processed dataset...", flush=True)
    df = pd.read_parquet(args.data).reset_index(drop=True)
    df["benchmark_row_id"] = np.arange(len(df), dtype=int)
    df = df.loc[(df["prot"] > 0) & (df["RNA"] > 0)
                & df["transcript_context"].notna()].copy().reset_index(drop=True)
    if args.limit is not None:
        df = df.iloc[:args.limit].copy()
    input_data = build_ostir_input(df)
    input_path = output_dir / "ostir_input.csv"
    raw_path = output_dir / "ostir_raw_output.csv"
    prediction_path = output_dir / "ostir_predictions.csv"
    dataset_path = output_dir / "validation_dataset.parquet"
    input_data.to_csv(input_path, index=False)
    df.to_parquet(dataset_path, index=False)

    print(f"Running OSTIR on {len(input_data):,} transcript contexts...", flush=True)
    run_ostir(input_path, raw_path, args.conda_env, args.threads)
    raw = pd.read_csv(raw_path)
    expected_ids = input_data["name"].astype(int).tolist()
    predictions = parse_ostir_output(raw, expected_ids, allow_missing=True)
    predictions.to_csv(prediction_path, index=False)
    missing_ids = sorted(set(expected_ids).difference(predictions["row_index"]))
    pd.DataFrame({"benchmark_row_id": missing_ids}).to_csv(
        output_dir / "ostir_unscored_rows.csv", index=False
    )
    metadata = {
        "input_rows": int(len(input_data)),
        "predicted_rows": int(len(predictions)),
        "unscored_rows": int(len(missing_ids)),
        "coverage": float(len(predictions) / len(input_data)),
        "ostir_env": args.conda_env,
        "ostir_version": subprocess.check_output(
            ["conda", "run", "-n", args.conda_env, "ostir", "--version"], text=True
        ).strip(),
        "viennarna_python_version": subprocess.check_output(
            ["conda", "run", "-n", args.conda_env, "python", "-c",
             "import RNA; print(RNA.__version__)"], text=True
        ).strip(),
        "viennarna_cli_version": subprocess.check_output(
            ["conda", "run", "-n", args.conda_env, "RNAfold", "--version"], text=True
        ).strip(),
        "predictions": prediction_path.name,
        "input": input_path.name,
    }
    (output_dir / "ostir_run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Saved {len(predictions):,} OSTIR predictions to {prediction_path}", flush=True)


if __name__ == "__main__":
    main()
