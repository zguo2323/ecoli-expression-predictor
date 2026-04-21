# BioE234 Final Project — Predictive Expression Modeling for E. coli

An MCP server that predicts protein expression levels in *E. coli* from promoter and RBS sequences. Trained on Kosuri et al. 2013 (~12,500 characterized constructs).

## Setup

```bash
pip install -r requirements.txt
```

Copy your Gemini API key into `.env`:
```
GEMINI_API_KEY=your_key_here
```

Place the raw data files in `data/raw/`:
- `sd01.xls` — promoter sequences and measured expression
- `sd02.xls` — RBS sequences and measured expression
- `sd03.xls` — all promoter×RBS construct combinations

## Usage

**Build the dataset:**
```bash
python -c "from modules.pipeline.pipeline import build_dataset; build_dataset('data/raw/sd01.xls', 'data/raw/sd02.xls', 'data/raw/sd03.xls')"
```

**Train the model:**
```bash
python -c "from modules.model.model import train_model; print(train_model('data/processed/constructs.parquet'))"
```

**Start the MCP server:**
```bash
python server.py
```

**Run tests:**
```bash
pytest tests/
```

## Project Structure

```
project/
├── data/
│   ├── raw/          # place sd01.xls, sd02.xls, sd03.xls here (not committed)
│   └── processed/    # constructs.parquet (generated)
├── modules/
│   ├── pipeline/     # Module A: data loading and joining
│   ├── features/     # Module B: biological feature extraction
│   └── model/        # Module C: XGBoost training and prediction
├── artifacts/        # trained model.pkl (generated, not committed)
├── tests/
├── server.py         # MCP server entry point
└── client_gemini.py  # Gemini LLM client
```

## Data Source

Kosuri et al. (2013). *Composability of regulatory sequences controlling transcription and translation in Escherichia coli.* PNAS.
