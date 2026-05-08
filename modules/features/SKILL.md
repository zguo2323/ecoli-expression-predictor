# Module B — Biological Feature Extraction

## Purpose
Translates raw DNA sequences into biologically meaningful numerical features. All functions are pure and deterministic — no randomness, no external calls.

## When to use
- Use `extract_all_features` when the user wants a complete breakdown of a promoter+RBS pair
- Use individual tools (`score_minus10_box`, `score_minus35_box`, `score_sd_sequence`, etc.) when the user asks about a specific element
- Use `compute_mrna_folding_energy` when the user asks about structural compatibility between a specific promoter and RBS

## Tools at a glance
| Tool | Input | What it tells you |
|---|---|---|
| `score_minus10_box` | promoter_seq | Strength of the −10 transcription element (0–1) |
| `score_minus35_box` | promoter_seq | Strength of the −35 transcription element (0–1) |
| `get_spacer_length` | promoter_seq | Distance between −35 and −10 (optimal: 17 nt) |
| `score_sd_sequence` | rbs_seq | Strength of the Shine-Dalgarno element (0–1) |
| `get_sd_spacing` | rbs_seq | Distance from SD end to start codon proxy (optimal: 5–10 nt) |
| `compute_gc_content` | seq | GC fraction (0–1) |
| `compute_mrna_folding_energy` | promoter_seq, rbs_seq | MFE of the mRNA junction (kcal/mol) — more negative means RBS is more occluded |
| `extract_all_features` | promoter_seq, rbs_seq | All of the above in one call |

## Notes
- Sequences can be any length and any case — functions handle uppercase conversion internally
- `compute_mrna_folding_energy` is context-dependent: the same RBS will score differently paired with different promoters
- A more negative `mrna_folding_energy` is generally worse for translation — the RBS is buried in secondary structure
- `spacer_optimal` is True for 15–21 nt; `sd_spacing_optimal` is True for 5–10 nt
- Prefer `extract_all_features` when doing construct comparison — it gives the full picture in one tool call
