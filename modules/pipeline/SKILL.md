# Module A — Data Pipeline

## Purpose
Builds the master dataset by loading and joining the three Kosuri et al. supplementary tables into a single parquet file.

## MCP exposure
Dataset construction is an offline project setup step and is not exposed as a conversational MCP tool. Use the documented Python setup command when rebuilding the training dataset.

## Notes
- Do not invoke dataset construction in response to user questions about sequences or expression.
- The generated parquet file is used by the offline training and evaluation scripts.
