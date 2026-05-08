# Module A — Data Pipeline

## Purpose
Builds the master dataset by loading and joining the three Kosuri et al. supplementary tables into a single parquet file.

## When to use
- Only needed if `data/processed/constructs.parquet` does not exist yet
- Call `build_dataset` once before any model training or evaluation

## Tool: `build_dataset`
- Inputs: paths to `sd01.xls`, `sd02.xls`, `sd03.xls`
- Output: path to the saved `constructs.parquet` file (~11,700 rows)
- This is a one-time setup step, not a per-query tool

## Notes
- Do not call this tool in response to user questions about sequences or expression — it is a data preparation utility, not an analysis tool
- The output parquet file is the input required by `train_model` and `evaluate_model`
