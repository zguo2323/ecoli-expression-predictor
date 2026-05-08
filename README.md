#    BioE234 Final Project — Predictive Expression Modeling for E. coli

An MCP server that predicts protein expression levels in *E. coli* from promoter and RBS sequences, trained on Kosuri et al. 2013 (~12,500 characterized constructs). A connected LLM (Gemini) can call the tools interactively to reason about construct design.

---

## Modules

### Module A — Data Pipeline (`modules/pipeline/`)

Loads and joins the three Kosuri et al. supplementary tables (`sd01.xls`, `sd02.xls`, `sd03.xls`) into a single clean dataset. Strips quoted string values, filters out flagged low-quality constructs, and merges promoter and RBS sequences with their measured expression levels. Outputs `data/processed/constructs.parquet` (~11,700 rows).

### Module B — Feature Extraction (`modules/features/`)

Translates raw DNA sequences into biologically meaningful numerical features used by the model. All functions are pure and deterministic. Features extracted:

- **−10 and −35 box scores** — sliding-window PWM match against the *E. coli* sigma-70 consensus elements (`TATAAT` and `TTGACA`), which RNA polymerase recognizes to initiate transcription
- **Spacer length** — distance between the −35 and −10 elements (optimal: 17 nt)
- **Shine-Dalgarno score and spacing** — PWM match against the ribosome binding consensus `AGGAGG` and its distance to the start codon (optimal: 5–10 nt)
- **GC content** — for both promoter and RBS sequences
- **mRNA junction folding energy** — minimum free energy (MFE) of the mRNA secondary structure formed at the junction between the promoter's 3′ end and the RBS, computed using ViennaRNA

> **Feedback from Demo Day:** Our professor noted that without modeling mRNA secondary structure, the RBS ranking would always return the same winner regardless of promoter context, because RBS features were computed in isolation. We addressed this by adding `compute_mrna_folding_energy`, which folds the specific promoter+RBS junction sequence using ViennaRNA. The same RBS can have different accessibility depending on what's upstream of it — a GC-rich promoter tail can form a stem-loop that occludes the Shine-Dalgarno sequence, suppressing translation even if the SD sequence itself is strong. This makes `rank_rbs_for_promoter` genuinely context-dependent: the top-ranked RBS will differ for different promoters.

### Module C — Model (`modules/model/`)

Trains a gradient-boosted model on the featurized dataset and exposes prediction and ranking as MCP tools.

**What is XGBoost and why do we use it?**
XGBoost (Extreme Gradient Boosting) is a machine learning algorithm that builds an ensemble of decision trees sequentially, where each new tree corrects the errors of the previous ones. It is well-suited for this task because:
- The relationship between sequence features and expression is highly non-linear and cannot be captured by a linear model
- XGBoost handles the mixed feature types we use (continuous scores, integer spacer lengths, boolean flags) without requiring normalization
- It performs well on tabular datasets of this size (~11,700 rows, 10 features) without overfitting
- It is fast to train and produces feature importance scores that help interpret which biological signals matter most

We train on an 80/10/10 train/val/test split stratified by expression level bin to ensure all parts of the expression range are represented. The trained model achieves a **Spearman r = 0.83** on the held-out test set, meaning it ranks constructs by expression level with high accuracy.

At inference time, the model uses the same 10 computed sequence features for any novel promoter+RBS pair — no precomputed data required.

---

## Setup

**Install dependencies:**
```bash
pip install -r requirements.txt
brew install libomp   # required for XGBoost on macOS
```

**Add your Gemini API key to `.env`:**
```
GEMINI_API_KEY=your_key_here
```

The raw data files are included in the repository under `data/raw/`:
- `sd01.xls` — promoter sequences and measured expression
- `sd02.xls` — RBS sequences and measured expression
- `sd03.xls` — all promoter×RBS construct combinations

---

## Running the Project

**Step 1 — Build the dataset** (produces `data/processed/constructs.parquet`):
```bash
python3 -c "
from modules.pipeline.pipeline import build_dataset
build_dataset('data/raw/sd01.xls', 'data/raw/sd02.xls', 'data/raw/sd03.xls')
"
```

**Step 2 — Train the model** (produces `artifacts/model.pkl`):
```bash
python3 -c "
from modules.model.model import train_model
train_model('data/processed/constructs.parquet')
"
```

**Step 3 — Run tests:**
```bash
pytest tests/
```

**Step 4 — Start the MCP server + Gemini client:**
```bash
python3 client_gemini.py
```

---

## Connecting to MCP Clients

The MCP server is client-agnostic — it can be connected to any MCP-compatible LLM client. Add the server to your client's configuration using the following pattern:

```json
{
  "mcpServers": {
    "ecoli-expression-predictor": {
      "command": "/path/to/your/.venv/bin/python3",
      "args": ["/path/to/bioe234-final-project/server.py"],
      "cwd": "/path/to/bioe234-final-project"
    }
  }
}
```

| Client | Config file location |
|---|---|
| **Claude Desktop** | `~/Library/Application Support/Claude/claude_desktop_config.json` |
| **Cursor** | `.cursor/mcp.json` in your project root |
| **VS Code (Copilot)** | `.vscode/mcp.json` in your project root |

After adding the config, restart the client and the 12 tools will be available automatically.

---

## Example Prompts

Once connected, try these:

```
"I have a promoter sequence TTGACAAATATTTGTATAAT. How strong is its -10 box? What about its -35 box?"

"Predict the protein expression level for promoter TTGACATATAATCC paired with RBS AAAGAGGAGAAA."

"Rank the top 5 RBS sequences for maximum protein output with promoter TTGACATATAATCC."

"Compare these two constructs: promoter TTGACATATAATCC + RBS AAAGAGGAGAAA versus promoter TTTACATATAATCC + RBS AGGATGT. Which will express more highly and why?"
```

---

## Sample Agent Run

The following is a real interaction with Claude using the MCP server, demonstrating cross-tool reasoning — Claude automatically chains multiple tool calls to answer a single question:

**You:** I designed a new synthetic promoter: TTGACATATAATCCGGAATTC. Is it likely to be strong? What's the best RBS to pair it with?

**Tool calls made by Claude:**
```
→ score_minus10_box(promoter_seq="TTGACATATAATCCGGAATTC")  →  1.0
→ score_minus35_box(promoter_seq="TTGACATATAATCCGGAATTC")  →  1.0
→ get_spacer_length(promoter_seq="TTGACATATAATCCGGAATTC")  →  0
→ rank_rbs_for_promoter(promoter_seq="TTGACATATAATCCGGAATTC", top_n=1)
  → [{"rbs_id": "B0030_RBS", "predicted_prot": 231757.36}]
```

**Claude:** While your promoter has perfect scores for both the −10 and −35 boxes (1.0 each), which suggests strong RNA polymerase recognition and transcription initiation signals, the spacer length between these elements is 0 nucleotides. An optimal spacer length is typically 17 nucleotides (functional range: 15–21 nt), so this significantly suboptimal spacing might reduce the overall strength of your promoter.

The best RBS to pair with your promoter for maximum protein output is **B0030_RBS** (`ATTAAAGAGGAGAAATTA CATATG`). This combination is predicted to yield a protein output of **231,757 units**.

---

## Project Structure

```
project/
├── data/
│   ├── raw/              # sd01.xls, sd02.xls, sd03.xls (Kosuri et al. supplementary data)
│   └── processed/        # constructs.parquet (generated by Step 1)
├── modules/
│   ├── pipeline/         # Module A: data loading and joining
│   ├── features/         # Module B: biological feature extraction
│   └── model/            # Module C: XGBoost training and prediction
├── artifacts/            # model.pkl (generated by Step 2, not committed)
├── tests/                # 29 unit tests across all three modules
├── server.py             # MCP server — scans modules/ and registers all tools
├── client_gemini.py      # Gemini LLM client
└── SPEC.md               # full project specification
```

---

## Data Source

Kosuri et al. (2013). *Composability of regulatory sequences controlling transcription and translation in Escherichia coli.* PNAS. [doi:10.1073/pnas.1301301110](https://doi.org/10.1073/pnas.1301301110)
