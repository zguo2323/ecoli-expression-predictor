# BioE234 Final Project – Predictive Expression Modeling for E. coli
## Spec for Claude Code

---

## Project Goal

Build an MCP server that exposes biologically meaningful, deterministic Python functions
for predicting protein expression levels in E. coli from promoter and RBS sequences.
The system is trained on Kosuri et al. 2013 (~12,500 characterized constructs) and
allows an LLM to reason about construct design interactively.

---

## Repo Structure

```
project/
├── SPEC.md                        # this file
├── README.md                      # to be written last
├── .env                           # GEMINI_API_KEY (never commit)
├── .gitignore                     # include .env, *.pkl, data/raw/
├── requirements.txt
│
├── data/
│   ├── raw/                       # place sd01.xls, sd02.xls, sd03.xls here
│   └── processed/                 # output of pipeline → constructs.parquet
│
├── modules/
│   ├── pipeline/
│   │   ├── pipeline.py            # Module A: data loading and joining
│   │   └── pipeline.json          # MCP wrapper for build_dataset
│   │
│   ├── features/
│   │   ├── features.py            # Module B: all feature extraction functions
│   │   ├── score_minus10_box.json
│   │   ├── score_minus35_box.json
│   │   ├── score_sd_sequence.json
│   │   ├── get_sd_spacing.json
│   │   ├── get_spacer_length.json
│   │   ├── compute_gc_content.json
│   │   ├── extract_all_features.json
│   │   └── compute_mrna_folding_energy.json   # ← NEW
│   │
│   └── model/
│       ├── model.py               # Module C: training, prediction, ranking
│       ├── predict_expression.json
│       ├── rank_rbs_for_promoter.json
│       └── evaluate_model.json
│
├── artifacts/
│   └── model.pkl                  # trained XGBoost model (generated, not committed)
│
├── server.py                      # MCP server entry point (scans modules/)
├── client_gemini.py               # Gemini client (provided by course starter)
│
└── tests/
    ├── test_pipeline.py
    ├── test_features.py
    └── test_model.py
```

---

## Requirements

```
# requirements.txt
pandas>=2.0
pyarrow
xlrd>=2.0.1
xgboost
scikit-learn
scipy
numpy
pytest
python-dotenv
google-generativeai
ViennaRNA        # ← NEW: mRNA secondary structure prediction
```

---

## Module A – Data Pipeline (`modules/pipeline/pipeline.py`)

### Purpose
Load and join the three Kosuri supplementary tables into a single clean DataFrame.
This DataFrame is the contract input for Module B.

### Data files (in `data/raw/`)
- `sd01.xls` – ~114 rows. Promoter sequences and their measured RNA/protein levels.
- `sd02.xls` – ~111 rows. RBS sequences and their measured RNA/protein/translation levels.
- `sd03.xls` – 12,655 rows, sheet name "Constructs". All promoter×RBS combinations
  with measured expression. Key columns: Promoter, RBS, RNA, prot, deltaG, bad.prot, bad.RNA.

### Column name cleaning
Strip surrounding quotes from string values in Promoter and RBS columns.
sd01/sd02 values are stored as `"pFAB124"` (with quotes) – strip them.

### Functions to implement

```python
def load_promoter_table(path: str) -> pd.DataFrame:
    """
    Load sd01.xls. Clean column names (strip quotes from string values).
    Return DataFrame with columns:
        promoter_id (str), sequence (str), mean_RNA (float),
        mean_prot (float), TSS_best (int)
    Rename: Promoter->promoter_id, Sequence->sequence,
            mean.RNA->mean_RNA, mean.prot->mean_prot, TSS.best->TSS_best
    Drop rows where sequence is null.
    """

def load_rbs_table(path: str) -> pd.DataFrame:
    """
    Load sd02.xls. Clean column names.
    Return DataFrame with columns:
        rbs_id (str), sequence (str), mean_RNA (float),
        mean_prot (float), mean_xlat (float)
    Rename: RBS->rbs_id, Sequence->sequence,
            mean.RNA->mean_RNA, mean.prot->mean_prot, mean.xlat->mean_xlat
    Drop rows where sequence is null.
    """

def load_construct_table(path: str) -> pd.DataFrame:
    """
    Load sd03.xls sheet "Constructs". Clean string values.
    Filter: keep only rows where bad.prot == False AND bad.RNA == False.
    Return DataFrame with columns:
        promoter_id (str), rbs_id (str), RNA (float),
        prot (float), deltaG (float)
    Rename: Promoter->promoter_id, RBS->rbs_id
    Drop rows where prot or RNA is null.
    """

def build_dataset(
    sd01_path: str,
    sd02_path: str,
    sd03_path: str,
    output_path: str = "data/processed/constructs.parquet"
) -> pd.DataFrame:
    """
    Join all three tables.
    Steps:
      1. load_promoter_table -> promoters_df (keep promoter_id, sequence as promo_seq)
      2. load_rbs_table -> rbs_df (keep rbs_id, sequence as rbs_seq)
      3. load_construct_table -> constructs_df
      4. merge constructs_df with promoters_df on promoter_id (left join)
      5. merge result with rbs_df on rbs_id (left join)
      6. drop rows where promo_seq or rbs_seq is null
      7. save to output_path as parquet
      8. return final DataFrame

    Final columns must include at minimum:
        promoter_id, rbs_id, promo_seq, rbs_seq,
        RNA, prot, deltaG
    """
```

### MCP wrapper (`pipeline.json`)
```json
{
  "id": "org.bioe234.build_dataset",
  "name": "Build Expression Dataset",
  "description": "Loads and joins the three Kosuri et al. supplementary tables into a single clean dataset of promoter-RBS constructs with measured expression levels. Returns path to saved parquet file.",
  "type": "function",
  "inputs": [
    {"name": "sd01_path", "type": "string"},
    {"name": "sd02_path", "type": "string"},
    {"name": "sd03_path", "type": "string"}
  ],
  "outputs": [{"name": "output_path", "type": "string"}],
  "execution_details": {
    "execution": "run",
    "mcp_name": "build_dataset",
    "seq_params": ["sd01_path", "sd02_path", "sd03_path"]
  }
}
```

---

## Module B – Feature Extraction (`modules/features/features.py`)

### Purpose
Translate raw DNA sequences into biologically meaningful numerical features.
Every function is pure and deterministic: same input → same output, always.
No external API calls. No randomness.

### Biological background (for implementation guidance)

**Promoter features:**
- The −35 element consensus in E. coli is `TTGACA`. RNA polymerase sigma factor
  recognizes this hexamer approximately 35 nt upstream of the transcription start site.
- The −10 element consensus is `TATAAT`. This is the most conserved element and
  most predictive of promoter strength.
- Optimal spacer between −10 and −35 elements is 17 nt (range 15–21 is functional).
- Score both elements using position weight matrix (PWM) scoring against consensus.

**RBS features:**
- The Shine-Dalgarno (SD) consensus is `AGGAGG`. It base-pairs with 16S rRNA.
- SD sequence typically sits 5–10 nt upstream of the AUG start codon.
- Score by finding the best-matching hexamer in the RBS sequence (sliding window).
- Spacing to start codon critically affects translation efficiency.

<mark>**mRNA secondary structure (NEW):**</mark>
<mark>- The same RBS sequence can be accessible or buried depending on what promoter is</mark>
<mark>  upstream of it. The junction between the promoter 3' end and the RBS folds into</mark>
<mark>  a specific mRNA secondary structure that gates ribosome access.</mark>
<mark>- Without modeling this, `rank_rbs_for_promoter` always returns the same ranking</mark>
<mark>  regardless of promoter context, because RBS features are computed in isolation.</mark>
<mark>- Use ViennaRNA (`RNA.fold`) on the last 30 nt of the promoter + full RBS sequence</mark>
<mark>  to compute the minimum free energy (MFE) of the junction region.</mark>
<mark>- More negative MFE = more tightly folded = less ribosome access = lower translation.</mark>

### Functions to implement

```python
def compute_gc_content(seq: str) -> float:
    """
    Returns GC fraction of seq (0.0 to 1.0).
    seq: DNA string, any case. Convert to uppercase internally.
    Return 0.0 for empty string.
    """

def score_pwm(seq: str, consensus: str) -> float:
    """
    Helper. Score best alignment of consensus against seq using sliding window.
    For each position in consensus, score = fraction of matching characters
    across the window, weighted by conservation.

    Simple implementation:
      - Slide a window of len(consensus) across seq
      - For each window, compute match score:
          score = sum(1 for i where window[i] == consensus[i]) / len(consensus)
      - Return max score across all windows
    Return 0.0 if len(seq) < len(consensus).
    """

def score_minus10_box(promoter_seq: str) -> float:
    """
    Score the −10 element (consensus: TATAAT) in promoter_seq.
    Uses score_pwm(promoter_seq, "TATAAT").
    Returns float 0.0–1.0.
    promoter_seq should be the full promoter sequence (~40–60 nt).
    """

def score_minus35_box(promoter_seq: str) -> float:
    """
    Score the −35 element (consensus: TTGACA) in promoter_seq.
    Uses score_pwm(promoter_seq, "TTGACA").
    Returns float 0.0–1.0.
    """

def get_spacer_length(promoter_seq: str) -> int:
    """
    Find the positions of the best −35 match and best −10 match,
    return the distance in nt between the end of the −35 element
    and the start of the −10 element.

    Implementation:
      - Find start index of best TTGACA match (pos_35)
      - Find start index of best TATAAT match (pos_10), requiring pos_10 > pos_35
      - spacer = pos_10 - (pos_35 + 6)
      - Return spacer. If either element not found or spacer < 0, return -1.
    """

def score_sd_sequence(rbs_seq: str) -> float:
    """
    Score the Shine-Dalgarno element (consensus: AGGAGG) in rbs_seq.
    Uses score_pwm(rbs_seq, "AGGAGG").
    Returns float 0.0–1.0.
    rbs_seq: the RBS sequence (~15–25 nt upstream of start codon).
    """

def get_sd_spacing(rbs_seq: str) -> int:
    """
    Return the distance in nt between the end of the best SD match
    and the end of rbs_seq (proxy for distance to start codon).

    Implementation:
      - Find start index of best AGGAGG match (pos_sd)
      - spacing = len(rbs_seq) - (pos_sd + 6)
      - Return spacing. If no match found, return -1.
    """
```

```python
# ---- NEW FUNCTION ----
def compute_mrna_folding_energy(promoter_seq: str, rbs_seq: str) -> float:
    """
    Compute the minimum free energy (MFE) of the mRNA junction region formed
    by the 3' end of the promoter and the RBS sequence.

    This captures context-dependent ribosome accessibility: the same RBS can
    be accessible or occluded depending on what promoter it is paired with.

    Implementation:
      - junction = last 30 nt of promoter_seq + full rbs_seq
      - Convert T -> U (DNA to RNA)
      - Call RNA.fold(junction) from the ViennaRNA package
      - Return the MFE (float, kcal/mol). More negative = more folded = less
        accessible. Return 0.0 if ViennaRNA is not available.
    """
```

```python
def extract_all_features(promoter_seq: str, rbs_seq: str) -> dict:
    """
    Calls all above functions and returns a single feature dict.
    This is the main tool the LLM will call for a full feature breakdown.

    Returns:
    {
        "gc_promoter": float,          # GC content of promoter
        "gc_rbs": float,               # GC content of RBS
        "score_minus10": float,        # −10 box PWM score
        "score_minus35": float,        # −35 box PWM score
        "spacer_length": int,          # nt between −35 and −10
        "spacer_optimal": bool,        # True if spacer is 15–21 nt
        "score_sd": float,             # Shine-Dalgarno PWM score
        "sd_spacing": int,             # nt from SD end to sequence end
        "sd_spacing_optimal": bool,    # True if sd_spacing is 5–10 nt
        "mrna_folding_energy": float,  # ← NEW: MFE of promoter+RBS junction (kcal/mol)
    }
    """
```

### MCP wrappers
Create one `.json` file per public function. Existing wrappers unchanged.

<mark>New wrapper to add:</mark>

```json
// compute_mrna_folding_energy.json  ← NEW
{
  "id": "org.bioe234.compute_mrna_folding_energy",
  "name": "mRNA Junction Folding Energy Calculator",
  "description": "Computes the minimum free energy (MFE) of the mRNA secondary structure formed at the junction between a promoter's 3' end and an RBS sequence. A more negative value means the RBS is more likely to be occluded by folding, reducing ribosome access and translation. Use this to evaluate whether a promoter-RBS pair is structurally compatible, not just sequence-compatible.",
  "type": "function",
  "inputs": [
    {"name": "promoter_seq", "type": "string"},
    {"name": "rbs_seq", "type": "string"}
  ],
  "outputs": [{"name": "mrna_folding_energy", "type": "number"}],
  "examples": [{
    "input": {
      "promoter_seq": "TTGACATATAATCCGG",
      "rbs_seq": "AAAGAGGAGAAATTTA"
    },
    "output": -3.2
  }],
  "execution_details": {
    "execution": "run",
    "mcp_name": "compute_mrna_folding_energy",
    "seq_params": ["promoter_seq", "rbs_seq"]
  }
}
```

---

## Module C – Model (`modules/model/model.py`)

### Purpose
Train a gradient boosted model on the feature-engineered dataset,
evaluate it rigorously, and expose prediction and ranking as MCP tools.

### Functions to implement

```python
def load_and_featurize(parquet_path: str) -> pd.DataFrame:
    """
    Load constructs.parquet, run extract_all_features on each row,
    append feature columns to the DataFrame.
    Also include deltaG from sd03 as a feature (it's precomputed, valid to use).
    Return augmented DataFrame.
    """

def train_model(
    parquet_path: str,
    target: str = "prot",
    model_output_path: str = "artifacts/model.pkl",
    test_size: float = 0.1,
    val_size: float = 0.1,
    random_state: int = 42
) -> dict:
    """
    Full training pipeline:
      1. load_and_featurize(parquet_path)
      2. Split into train/val/test (80/10/10) → stratify by expression bin
         (create 5 equal-frequency bins of target, use as stratification key)
      3. Train XGBoost regressor on train set
         Suggested params: n_estimators=300, max_depth=5, learning_rate=0.05,
                           subsample=0.8, colsample_bytree=0.8
      4. Evaluate on test set → compute spearman_r, r2, mae
      5. Save model to model_output_path with pickle
      6. Return eval dict:
         {"spearman_r": float, "r2": float, "mae": float,
          "n_train": int, "n_test": int, "target": str}
    """

def predict_expression(
    promoter_seq: str,
    rbs_seq: str,
    model_path: str = "artifacts/model.pkl"
) -> dict:
    """
    Pure prediction function – the main MCP tool.
    Steps:
      1. extract_all_features(promoter_seq, rbs_seq)
         [includes mrna_folding_energy, computed fresh for novel sequences]
      2. Load model from model_path
      3. Predict on feature vector
      4. Compute a simple confidence interval:
         CI = [prediction * 0.75, prediction * 1.25]  (placeholder ±25%)
    Returns:
    {
        "predicted_prot": float,
        "confidence_interval": [float, float],
        "features_used": dict   # output of extract_all_features for transparency
    }
    Note: deltaG from sd03 is NOT available for novel sequences.
    Use mrna_folding_energy (computed via ViennaRNA) as the structural
    feature for novel inputs instead of the sd03 deltaG placeholder.
    """

def rank_rbs_for_promoter(
    promoter_seq: str,
    rbs_library_path: str = "data/raw/sd02.xls",
    model_path: str = "artifacts/model.pkl",
    top_n: int = 5
) -> list[dict]:
    """
    Score every RBS in the characterized library against a given promoter.
    Steps:
      1. Load rbs_library from rbs_library_path (use load_rbs_table)
      2. For each RBS, call predict_expression(promoter_seq, rbs_seq)
         [mrna_folding_energy is now computed per promoter+RBS pair,
          making the ranking context-dependent]
      3. Sort by predicted_prot descending
      4. Return top_n as list of dicts:
         [{"rbs_id": str, "rbs_seq": str, "predicted_prot": float}, ...]
    """

def evaluate_model(
    parquet_path: str,
    model_path: str = "artifacts/model.pkl"
) -> dict:
    """
    Rerun evaluation on the held-out test set.
    Use the same random_state=42 split as train_model to recover the test set.
    Return {"spearman_r": float, "r2": float, "mae": float, "n_test": int}
    """
```

### MCP wrappers
Unchanged from original spec.

---

## Tests (`tests/`)

### `test_pipeline.py` — unchanged

### `test_features.py`

```python
# All original tests unchanged, plus:

def test_compute_mrna_folding_energy_returns_float():          # ← NEW
    energy = compute_mrna_folding_energy("TTGACATATAATCCGG", "AAAGAGGAGAAATTTA")
    assert isinstance(energy, float)

def test_compute_mrna_folding_energy_context_dependent():      # ← NEW
    # Two different promoters paired with the same RBS should give different MFEs
    energy1 = compute_mrna_folding_energy("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    energy2 = compute_mrna_folding_energy("GCGCGCGCGCGCGCGC", "AAAGAGGAGAAA")
    assert energy1 != energy2

def test_extract_all_features_keys():
    result = extract_all_features("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    expected_keys = {
        "gc_promoter", "gc_rbs", "score_minus10", "score_minus35",
        "spacer_length", "spacer_optimal", "score_sd",
        "sd_spacing", "sd_spacing_optimal",
        "mrna_folding_energy",    # ← NEW
    }
    assert expected_keys == set(result.keys())
```

### `test_model.py` — unchanged

---

## Implementation Notes for Claude Code

1. **Quote stripping**: The Promoter and RBS values in the xls files are stored
   as `"pFAB124"` (with surrounding double quotes as part of the string).
   Strip them everywhere: `str(val).strip('"')`.

2. **deltaG in model**: The `deltaG` column in sd03 is precomputed RBS–ribosome
   ΔG. Include it as a training feature. <mark>For novel sequence prediction,</mark>
   <mark>use `mrna_folding_energy` (ViennaRNA) instead of the sd03 deltaG mean placeholder.</mark>
   <mark>The model should be trained with both deltaG (for training rows) and</mark>
   <mark>mrna_folding_energy (computed for all rows) as separate features.</mark>

3. **Module boundaries**: Modules import from each other in one direction only:
   `model.py` imports from `features.py`. `features.py` imports nothing from
   the other modules. `pipeline.py` imports nothing from features or model.

4. **MCP class structure**: Each module's functions should be wrapped in a class
   with `initiate(self)` and `run(self, **kwargs)` methods, following the
   course's gc_content.py pattern. The `run` method dispatches to the
   appropriate function based on the json wrapper's `mcp_name`.

5. **No randomness in feature functions**: `features.py` must be fully
   deterministic. Never use random seeds, sampling, or external calls inside
   any function in that file. <mark>ViennaRNA's `RNA.fold` is deterministic — safe to use.</mark>

6. **Parquet over CSV**: Always save processed data as `.parquet` (via pandas).
   It preserves dtypes and is faster to load.

7. **Model persistence**: Save the trained model as a dict containing both the
   model object and the training metadata:
   ```python
   pickle.dump({
       "model": xgb_model,
       "feature_names": feature_cols,
       "deltaG_mean": float,
       "eval": eval_dict
   }, f)
   ```

---

## Definition of Done

- [ ] `build_dataset` runs end-to-end and produces `constructs.parquet` with >10,000 rows
- [ ] All feature functions pass their unit tests
- [ ] `train_model` completes and saves `artifacts/model.pkl`
- [ ] `evaluate_model` returns spearman_r > 0.3
- [ ] All 8 MCP wrappers (json files) are present and valid  ← NEW count
- [ ] `pytest tests/` passes with no failures
- [ ] MCP server starts without errors (`python server.py`)
- [ ] All 3 cross-tool demo prompts produce a response that calls at least 2 tools
- [ ] <mark>`rank_rbs_for_promoter` returns a different ranking for different promoters</mark>
  <mark>(validates context-dependence via mrna_folding_energy)</mark>
