# Reproduce the public-data engine

The engine uses no private genome. It accepts one gene, up to ten named HGVS variants, and HPO identifiers. If the gene is omitted, Ensembl VEP resolves the small named-variant list. When `--gene` is provided, variants remain local. The complete public HPO ontology is downloaded once; all phenotype lookup is local. Neither HPO inputs nor disease labels are sent to an external API.

From the repository root:

```bash
uv sync --locked
uv run mva-repurpose --gene BUB1B \
  --variant 'ENST00000287598.11:c.2210T>G' \
  --variant 'ENST00000287598.11:c.3006T>G' \
  --hpo HP:0002859 --hpo HP:0000121 --hpo HP:0004322 \
  --hpo HP:0001508 --hpo HP:0003202 --hpo HP:0001622 \
  --hpo HP:0001518 \
  --mechanism loss_of_function \
  --disease 'mosaic variegated aneuploidy syndrome' \
  --cache examples/public_api_cache --offline \
  --output work/repurpose_bub1b_reproduced

uv run mva-repurpose --gene TSC2 --hpo HP:0001250 \
  --mechanism loss_of_function --disease 'tuberous sclerosis' \
  --cache examples/public_api_cache --offline \
  --output work/repurpose_positive_control_reproduced
```

For a fresh live run, omit `--offline` and use a new writable cache, such as `--cache work/new-public-cache`. New genes require internet access. `--output` is a prefix for JSON, CSV and Markdown artifacts. `python -m mva_repurpose` exposes the same CLI. The repository lock file fixes dependencies. Each response is stored with its request, retrieval timestamp, HTTP status and SHA-256 digest. Offline mode raises an error on any missing or altered cache entry.

## What the engine does

1. Resolve the gene with Ensembl and obtain its leaf Reactome pathways from Open Targets.
2. Retrieve public pathway proteins. Retain at most 20 leaf pathways and 200 proteins by default. Rank discovery priority by shared pathway count, then accession. Follow one explicit human protein-complex cofactor edge, capped at 50 extra proteins. Counts and any truncation are reported in JSON.
3. Retrieve human ChEMBL target mechanisms. This uses curated mechanism records, not arbitrary binding-assay hits. Complex membership is explicit because ChEMBL may attach a drug to its binding cofactor rather than the downstream enzyme.
4. Keep only molecules with `max_phase=4` and no withdrawn flag. This is the specified ChEMBL phase-4 proxy; it does not prove current availability in a particular country. Retrieve phase-4 indications separately.
5. Exclude direct target inhibition for a loss-of-function disease, and direct activation for a gain-of-function disease. Shared pathways do not establish rescue direction. A same-pathway drug remains grade E unless a direction-compatible direct action or supplied literature bridge adds evidence. An exact phase-4 disease-label match receives grade L for retrieval; its route, formulation, disease manifestation and age applicability remain unknown. It is not promoted to a rescue hypothesis.
6. Score the evidence components below and preserve all provenance. Abstain from a rescue recommendation when only graph proximity exists.

HPO terms are resolved and preserved as context. They do not receive a speculative phenotype-to-drug score. This is a limitation: the engine does not yet reason over phenotype mechanisms or estimate efficacy. It makes those missing links visible.

## Scores and grades

| Component | Points | Meaning |
|---|---:|---|
| Approval | +10 | ChEMBL phase 4, withdrawn flag false |
| Mechanism link | +20 / +5 / 0 | Compatible direct action / shared graph / no graph link |
| Evidence | +30 / 0–25 / 0 | Exact phase-4 disease-label match / supplied literature priority / absent |
| Pediatric use | +5 / 0 | Explicit supplied evidence / unknown |
| PGx | −10 / −3 / 0 | Actionable adjustment / metabolizer review flag / no supported penalty |
| Direct direction conflict | −100 | Action would worsen the supplied gene mechanism |
| Negative evidence | −15 per supplied negative bridge, capped at −50 | Conflicting supplied literature |

These weights are transparent review priorities, not calibrated treatment probabilities. A normal-metabolizer row receives no PGx penalty. A missing row is unknown. A matched PGx row is never a dosing instruction. Grade L means a known disease-label match with route, manifestation and age unverified; C means a supplied complete functional-rescue bridge, H means an unverified laboratory hypothesis with missing evidence, D means compatible direct target direction, E means context only, and X means a direct direction conflict. Equal scores sort alphabetically.

The frozen BUB1B run finds 33 compounds with a ChEMBL phase-4 proxy in its pathway context and zero supported rescue hypotheses. Many pathway neighbors are antimitotics. Their appearance is a warning against using pathway overlap as a treatment recommendation. The engine does not manufacture a top ten of treatments from these records.

The TSC2 retrieval benchmark finds everolimus first among 48 graph compounds with a ChEMBL phase-4 proxy without being given its name. The phase-4 TSC disease label supplies the retrieval signal. This is recovery of a known association, not a novel prediction. Both everolimus and sirolimus receive grade L, and neither becomes an automatically supported rescue hypothesis. Broad labels lose route and manifestation: [topical sirolimus (Hyftor)](https://www.ema.europa.eu/en/medicines/human/EPAR/hyftor) for TSC-associated facial angiofibromas does not establish a systemic seizure indication. EXIST-3 is a published randomized trial in TSC-associated treatment-resistant focal seizures: French et al., 2016, **PMID 27613521**, **DOI [10.1016/S0140-6736(16)31419-2](https://doi.org/10.1016/S0140-6736(16)31419-2)**. This verifies recovery of an established gene/pathway-to-approved-drug case. It does not validate efficacy for a new disease, formulation, manifestation or age group. Current jurisdiction-specific approval and pediatric suitability remain unverified by the engine. The negative control changes the supplied disease to MVA; everolimus drops to grade E and the engine abstains.

## Add reviewed literature or PGx

Use `--pgx path/to/your_pgx_table.csv` for a small table with `gene`, `phenotype`, `candidate_drugs_affected` and optional `clinical_action`, `recommendation` and `status`. The legacy `candidate_drugs` column is also accepted. Complete drug names are matched inside slash lists and explanatory prose using explicit ChEMBL synonyms and parent-molecule links, so vincristine can match vincristine sulfate. Unmatched input rows are reported in JSON. Drug classes are not expanded into guessed members. JSON rows or a TSV table are also accepted. No genotype-scale file should be supplied.

Use `--evidence path/to/bridges.json` to add explicit, reviewed bridges. Each record requires `gene`, `drug`, `direction` (`rescue`, `against`, `unknown`, or `hypothesis`), `pmids` and an integer `evidence_score` from 0–25. Optional `rationale`, `limitations`, `evidence_level` and `pediatric_use` fields remain in the artifact. The engine retrieves citation metadata through Europe PMC and checks identity. It cannot prove that the curator's causal interpretation is sound. An indirect hypothesis or incomplete rescue bridge stays grade H. Grade C requires an explicit supplied functional-rescue claim with no listed missing bridges; this remains a curator assertion, not automatic clinical validation. A PMID establishes neither efficacy nor clinical eligibility.

An empty template is in `evidence-schema.json`. No BUB1B treatment assertion is embedded in the engine or tests.

## Tests and source records

`uv run pytest tests/test_repurpose.py` runs four behavioral tests using cached public responses and synthetic boundary inputs. They forbid network requests. They test BUB1B abstention despite known antimitotic graph neighbors, approval filtering, eight local HPO lookups, everolimus discovery in TSC2/TSC, the wrong-disease control, and PGx review scoring. The boundary tests also check inverse agonists, ambiguity of partial agonists, per-edge causal flags, the clinical PGx column schema with a salt/parent synonym, and the curated-evidence caps and clinical gates.

`public_api_cache/` contains only public responses; no VCF, sample identifier, patient genotype table or patient-specific phenotype request is included. `freeze_public_cache.py` recreates the subset from a live local cache and the two result manifests. It does not fabricate responses. `results/repurpose_references.json` records opened bibliographic metadata.

Scientific methods references, opened and checked against Europe PMC core records:

- Reactome: Milacic et al., 2024. **PMID 37941124**, **DOI [10.1093/nar/gkad1025](https://doi.org/10.1093/nar/gkad1025)**.
- ChEMBL: Zdrazil et al., 2024. **PMID 37933841**, **DOI [10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004)**.
- HPO: Gargano et al., 2024. **PMID 37953324**, **DOI [10.1093/nar/gkad1005](https://doi.org/10.1093/nar/gkad1005)**.


## Regulatory hold and input limitations

Ataluren is removed from the eligible/context list despite its cached ChEMBL phase-4 value. The [EMA Translarna EPAR](https://www.ema.europa.eu/en/medicines/human/EPAR/translarna), opened on 2026-10-24 UTC, records that the European Commission did not renew its conditional authorization on **28 March 2025**. An alternative current eligible jurisdiction has not been verified. The engine records this explicit hold; it does not claim worldwide withdrawal. This regulatory source has no PMID or DOI.

With `--gene` provided, named variants are provenance only: their functional effects do not alter drug scores. HPO terms supply local labels/context only. Pediatric suitability is currently unknown unless explicit supplied evidence establishes it. Phase-4 data remain a proxy for every other graph compound, not verification of current market authorization.

Patient-specific integration rerun:

```bash
uv run mva-repurpose --gene BUB1B \
  --variant 'ENST00000287598.11:c.2210T>G' \
  --variant 'ENST00000287598.11:c.3006T>G' \
  --hpo HP:0002859 --hpo HP:0000121 --hpo HP:0004322 \
  --hpo HP:0001508 --hpo HP:0003202 --hpo HP:0001622 \
  --hpo HP:0001518 \
  --mechanism loss_of_function \
  --disease 'mosaic variegated aneuploidy syndrome' \
  --pgx path/to/your_pgx_table.csv \
  --cache examples/public_api_cache --offline \
  --output work/repurpose_bub1b_pgx_reproduced
```
