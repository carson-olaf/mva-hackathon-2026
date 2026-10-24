# Public-data repurposing engine: BUB1B

No direction-supported rescue hypothesis with a ChEMBL phase-4 proxy was established. Graph proximity or a broad disease-label match alone is insufficient.

| Rank | Drug with ChEMBL phase-4 proxy | Score | Grade | Evidence status |
|---|---|---:|---|---|
| 1 | ABEMACICLIB | 15 | E | pathway_context_only |
| 2 | ALECTINIB HYDROCHLORIDE | 15 | E | pathway_context_only |
| 3 | BORTEZOMIB | 15 | E | pathway_context_only |
| 4 | BORTEZOMIB D-MANNITOL | 15 | E | pathway_context_only |
| 5 | BRENTUXIMAB VEDOTIN | 15 | E | pathway_context_only |
| 6 | CABAZITAXEL | 15 | E | pathway_context_only |
| 7 | CARFILZOMIB | 15 | E | pathway_context_only |
| 8 | CERITINIB | 15 | E | pathway_context_only |
| 9 | COLCHICINE | 15 | E | pathway_context_only |
| 10 | CRIZOTINIB | 15 | E | pathway_context_only |

Grades: L = known disease-label match, with route/manifestation/age unverified; C = supplied literature bridge; D = direct target with compatible direction; E = graph context only; X = incompatible direct action.

Scores are review priorities. They are not efficacy probabilities. Grade E drugs are not rescue candidates. Grade L is a retrieval benchmark, not a novel prediction or automatic rescue recommendation. Equal scores are sorted alphabetically; their displayed order is not biological evidence.

## PGx integration

Input rows: 0; matched rows: 0. Unmatched rows are explicit in JSON and do not establish normal metabolism.

## Limitations

- No patient genome is required or uploaded; input is a small named-variant list and HPO terms.
- Phase-4 ChEMBL status is a market-approval proxy; jurisdiction and current availability require regulator verification.
- Exact disease-label matches are retrieval benchmarks, not novel predictions; route, formulation, manifestation and age may differ (for example topical sirolimus for TSC facial angiofibromas).
- With --gene supplied, named variants are provenance only; HPO labels/context do not drive scores.
- Shared Reactome membership does not establish causal direction, tissue efficacy, or genotype rescue.
- HPO terms are resolved for provenance, but do not assign unvalidated phenotype-to-drug scores.
- Missing PGx matches are unknown, never proof of normal metabolism; any matched PGx needs manual review.
- Pediatric use is unknown unless a supplied evidence bridge explicitly documents it.
- Literature bridge citations are retrieved; causal assertions remain the curator's responsibility.
- Evidence scores are transparent review priorities, not calibrated treatment probabilities.
