# Module C — Expression Model

## Purpose
Predicts protein expression levels for novel promoter+RBS pairs and ranks RBS candidates for a given promoter. Powered by an XGBoost model trained on ~11,700 Kosuri et al. constructs (Spearman r = 0.83 on held-out test set).

## When to use
- Use `predict_expression` when the user provides a specific promoter+RBS pair and wants a predicted expression level
- Use `rank_rbs_for_promoter` when the user has a fixed promoter and wants to find the best RBS candidates from the characterized library
- For construct comparison questions, call `predict_expression` for each construct and compare

## TSS follow-up
- Before calling `predict_expression` or `rank_rbs_for_promoter`, if `tss_best` is missing and the user has not already answered for this construct context, ask once: “Do you have an experimentally measured `tss_best` for a matching construct context? It is an integer offset from the promoter/RBS junction. If you have it, please provide it; otherwise I can continue without TSS-based accessibility features.”
- If the user supplies a measured value, pass it through. If they say they do not have one, are unsure, or ask to proceed, omit it and continue; do not infer or invent a measured TSS.
- Do not repeat the question for the same construct context after the user has answered.

## Tools at a glance
| Tool | Key inputs | Output |
|---|---|---|
| `predict_expression` | promoter_seq, rbs_seq, optional tss_best | predicted_prot, features_used, translation_context_available |
| `rank_rbs_for_promoter` | promoter_seq, top_n | ranked list of RBS sequences by predicted expression |

## Notes
- `predict_expression` accepts candidate sequences, but its reliability on novel promoter/RBS identities is lower than on familiar library parts; consult the grouped holdout results documented in the project reports
- `rank_rbs_for_promoter` scores the characterized RBS sequences from the Kosuri library against the given promoter
- Pass `tss_best` only when it is experimentally measured for a sufficiently matching construct context. If unknown, omit it; the tool marks translation-context features unavailable. Do not infer or invent a measured TSS
- Predictions are in the Kosuri dataset's `prot` target scale. They are not absolute protein concentrations and do not include a calibrated confidence interval; use them to prioritize and compare candidates
- When a user asks which construct will "express more," call `predict_expression` for each and compare `predicted_prot`; also call `extract_all_features` on both to explain the sequence-level reasons
