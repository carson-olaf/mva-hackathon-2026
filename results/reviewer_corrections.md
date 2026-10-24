# External-review corrections — completed locally

All eight items were applied as separate commits. No push or submission. The report, summary and pitch use the revised framing.

| Item | Local commit | Result |
|---|---|---|
| 1 | `8c1f7b3` | Requested title and opening retained verbatim. |
| 2 | `3a1beb7` | Both biological transfers and corrected-control limitation are on PDF page 1. |
| 3 | `f5ca628` | Exposure values unchanged; post-peak period, baseline-only sensitivity and reference concentration qualified. |
| 4 | `a89692d` | Inconclusive routes deprioritized; benchmark RMSE is not a variant CI; 18 runs are robustness checks. |
| 5 | `dc4c3d6` | Pandey is the existing DOI, cited once in the main bibliography; opposing findings and metabolic measurements retained. |
| 6 | `e8ddc01` | Burrell checkpoint comparison and NCIH508 nonresponse verified in primary main text. |
| 7 | `96bdafa` | Qualified novelty statement added; AI disclosure and chronology unchanged. |
| 8 | `3e08bb7` | Prior routes and genome detail moved to Appendices A–B. |

**Verification:** 28 tests pass in a freshly synchronized offline environment; the CLI runs and appropriately abstains from a supported BUB1B rescue claim. Exposure replay is byte-identical and every pre-review exposure value is unchanged. The PDF has 11 pages and the summary has one. The renderer now ends the copied acknowledgement at the next section, so appendices are not copied into the Track 1 methods. See [validation manifest](submission_validation.json).

**Literature checks:** [Pandey primary full paper, DOI 10.1093/nar/gkaf1035](https://www.diva-portal.org/smash/get/diva2:2010300/FULLTEXT01.pdf), Figures 1 and 3; [wording correction, DOI 10.1093/nar/gkag532](https://pmc.ncbi.nlm.nih.gov/articles/PMC13187833/); [Burrell primary manuscript, PMID 23446422 / DOI 10.1038/nature11935](https://pmc.ncbi.nlm.nih.gov/articles/PMC4636055/), main text before Figures 1 and 4.

**Limits:** no requested literature correction remains unverified. The exact earlier phrase “a connection of published facts, not a discovery” was absent from the current report; the requested contribution statement was added. Uridine-alone sufficiency, a nucleotide-responsive component in this genotype, patient exposure and long-term cancer safety remain unproved. The Burrell dose-assignment/component questions remain unresolved as recorded in the supplement audit. No clinician signoff or patient experiment is claimed.

**Attribution:** DGIdb bug identification and repair remain the maintainers' work. Our contribution is historical quantification, rare-disease exposure analysis and reproducible regression checks. The 11 remaining discrepancies stay unresolved.

The final independent review found no blocking contradiction. Its two precision fixes were applied: a failed model experiment demonstrates absence of response only under the tested conditions, and the summary labels the 11 discrepancies unresolved.
