# Module C — Expression Model

## Purpose
Predicts protein expression levels for novel promoter+RBS pairs and ranks RBS candidates for a given promoter. Powered by an XGBoost model trained on ~11,700 Kosuri et al. constructs (Spearman r = 0.83 on held-out test set).

## When to use
- Use `predict_expression` when the user provides a specific promoter+RBS pair and wants a predicted expression level
- Use `rank_rbs_for_promoter` when the user has a fixed promoter and wants to find the best RBS candidates from the characterized library
- Use `evaluate_model` when the user asks about model reliability or performance metrics
- For construct comparison questions, call `predict_expression` for each construct and compare

## Tools at a glance
| Tool | Key inputs | Output |
|---|---|---|
| `predict_expression` | promoter_seq, rbs_seq | predicted_prot, confidence_interval, features_used |
| `rank_rbs_for_promoter` | promoter_seq, top_n | ranked list of RBS sequences by predicted expression |
| `evaluate_model` | parquet_path | spearman_r, r2, mae, n_test |

## Notes
- `predict_expression` works on any novel sequence — it does not require the sequence to be in the training library
- `rank_rbs_for_promoter` scores all ~110 RBS sequences from the Kosuri library against the given promoter; rankings are context-dependent because mRNA folding energy is computed per promoter+RBS pair
- The confidence interval is ±25% of the predicted value — treat it as an order-of-magnitude guide, not a precise bound
- Predicted protein values are in the same relative units as the Kosuri et al. dataset; use them for ranking and comparison, not as absolute expression levels
- When a user asks which construct will "express more," call `predict_expression` for each and compare `predicted_prot`; also call `extract_all_features` on both to explain the sequence-level reasons
