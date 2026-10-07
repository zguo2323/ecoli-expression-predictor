# E. coli Expression Predictor

An MCP server for **pre-experiment prioritization** of *E. coli* promoter–RBS designs. It uses an XGBoost model trained on the Kosuri et al. 2013 promoter/RBS library and exposes sequence-feature explanations alongside candidate predictions.

The output is the model's prediction of the dataset's `prot` target. It is useful for comparing and prioritizing candidates in the same reporter and experimental context; it is **not an absolute protein concentration, a probability, or a calibrated confidence interval**. Predictions for promoter/RBS identities absent from training are less certain, especially when both are new.

## What the model uses

The model combines promoter features (−10/−35 consensus scores, spacer length and GC content), RBS features (Shine–Dalgarno score and spacing, GC content), and transcript-context accessibility features. The latter estimate Shine–Dalgarno and start-codon accessibility by ViennaRNA folding from a measured transcription start site (TSS) through a provisional 90-nt sfGFP coding context.

For an untested design, a measured `tss_best` is usually unavailable. If it is missing, the assistant asks once whether you have a measured value for a sufficiently matching construct context. If you provide one, it is used; if you do not have one, are unsure, or want to proceed, prediction and ranking continue without TSS-derived accessibility features. The value is an integer offset from the promoter/RBS junction; do not guess one. The sfGFP reference sequence is provisional and has not been verified base-by-base against the Addgene plasmid record. See the [transcript-context data dictionary](docs/transcript_context_data_dictionary.md).

Biological validation compares sequence-only and accessibility models against baselines using grouped cross-validation, group bootstrap, held-out feature-group permutation, and out-of-fold residual analysis. Results and limitations are summarized in [`reports/biological_context/summary.md`](reports/biological_context/summary.md).

## MCP tools

The server exposes nine tools. The feature tools report interpretable sequence descriptors; the two model tools support design decisions.

| Tool | Inputs | Use |
|---|---|---|
| `predict_expression` | `promoter_seq`, `rbs_seq`, optional measured `tss_best` | Predict the dataset `prot` target for one construct; returns the features used and whether transcript-context features were available. |
| `rank_rbs_for_promoter` | `promoter_seq`, `top_n`, optional measured `tss_best` | Rank candidate RBSs from the characterized Kosuri library for a promoter. |
| `extract_all_features` | `promoter_seq`, `rbs_seq`, optional measured `tss_best` | Return promoter/RBS sequence features and, when possible, transcript accessibility estimates. |
| `score_minus10_box` | `promoter_seq` | Score the −10 consensus match. |
| `score_minus35_box` | `promoter_seq` | Score the −35 consensus match. |
| `get_spacer_length` | `promoter_seq` | Estimate the distance between the best −35 and −10 matches. |
| `score_sd_sequence` | `rbs_seq` | Score the Shine–Dalgarno consensus match. |
| `get_sd_spacing` | `rbs_seq` | Estimate SD-to-start spacing. |
| `compute_gc_content` | `seq` | Calculate GC fraction. |

The old promoter-tail MFE calculator is not exposed: promoter sequence upstream of the measured TSS is not a faithful RNA folding context. Dataset building and model evaluation are offline project workflows, not conversational MCP tools. Their Python functions and evaluation scripts remain available in the repository.

## Setup

Requires Python 3.10+ and the dependencies in `requirements.txt`. On macOS, XGBoost may also require OpenMP (`brew install libomp`).

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

The repository includes the Kosuri supplementary files in `data/raw/`. Build the processed dataset and train the model once:

```bash
.venv/bin/python -c "from modules.pipeline.pipeline import build_dataset; build_dataset('data/raw/sd01.xls', 'data/raw/sd02.xls', 'data/raw/sd03.xls')"
.venv/bin/python -c "from modules.model.model import train_model; train_model('data/processed/constructs.parquet')"
```

This creates `data/processed/constructs.parquet` and `artifacts/model.pkl`. The MCP model tools need the trained artifact at the default path.

## Connect an MCP client

Add the server to your MCP client's configuration, replacing the paths with the absolute paths on your machine:

```json
{
  "mcpServers": {
    "ecoli-expression-predictor": {
      "command": "/absolute/path/to/ecoli-expression-predictor/.venv/bin/python",
      "args": ["/absolute/path/to/ecoli-expression-predictor/server.py"],
      "cwd": "/absolute/path/to/ecoli-expression-predictor"
    }
  }
}
```

Restart the client after saving the configuration. `server.py` loads the tool definitions from `modules/` and registers the nine tools above.

To use the optional Gemini demo client instead of another MCP host, copy `.env.example` to `.env`, add a Gemini API key, then run:

```bash
.venv/bin/python client_gemini.py
```

## Example prompts

- “Compare these promoter and RBS candidates using predicted expression and explain which sequence features differ.”
- “Rank the top five RBSs from the characterized library for this promoter.”
- “Show the −10/−35 scores and spacer length for this promoter.”
- “Compare Shine–Dalgarno score and spacing for these RBSs.”

For untested constructs, the assistant should omit `tss_best` unless a matching measured value is available. Treat returned values as prioritization signals and confirm promising designs experimentally.

## Evaluation and project structure

The main assessment uses row-random, promoter-held-out, RBS-held-out, and double-unseen splits so familiar-part performance is not confused with generalization to new parts. It also includes Ridge and OSTIR comparisons. Results, split definitions, feature importance, error analysis, and limitations are in the [biological-context benchmark summary](reports/biological_context/summary.md). Run the project tests with:

```bash
.venv/bin/python -m pytest -p no:capture tests/
```

```text
modules/features/     sequence feature extraction and MCP schemas
modules/model/        model training, prediction, ranking, and MCP schemas
modules/pipeline/     offline data preparation
scripts/              grouped validation and benchmark workflows
reports/              saved metrics, importance, and error analyses
docs/                 stage log, biological context, and project status
server.py             MCP server
client_gemini.py      optional Gemini example client
```

## Data source

Kosuri et al. (2013), “Composability of regulatory sequences controlling transcription and translation in *Escherichia coli*,” *PNAS*. [doi:10.1073/pnas.1301301110](https://doi.org/10.1073/pnas.1301301110).
