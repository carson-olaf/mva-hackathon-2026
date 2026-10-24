# Curated BUB1B laboratory hypotheses

This branch combines the external ChatGPT Pro mechanism review with the actual public-data CLI and the local marker-level PGx table. It does not replace the automatic graph-context run. Arimoclomol, azithromycin and prescription nicotinic acid are **indirect laboratory hypotheses**, with no demonstrated rescue of this genotype. The structural and LINCS analyses add zero positive rescue lines.

The evidence bundle names every missing bridge. A source that demonstrates one step in another protein, species or disease does not complete a drug-to-functional-BubR1 path. All three candidates have E capped at 1. No functional-rescue or clinical gate is passed by their ordinal priority. Expected clinical effects remain unknown.

## Reproduce the actual run

```bash
uv run mva-repurpose --gene BUB1B \
  --variant 'ENST00000287598.11:c.2210T>G' \
  --variant 'ENST00000287598.11:c.3006T>G' \
  --hpo HP:0002859 --hpo HP:0000121 --hpo HP:0004322 \
  --hpo HP:0001508 --hpo HP:0003202 --hpo HP:0001622 \
  --hpo HP:0001518 \
  --mechanism loss_of_function \
  --disease 'mosaic variegated aneuploidy syndrome' \
  --evidence examples/bub1b_curated_evidence.json \
  --approval-audit results/repurpose_curated_approval_audit.json \
  --pgx path/to/your_pgx_table.csv \
  --cache examples/bub1b_curated_cache \
  --cache-fallback examples/public_api_cache --offline \
  --output work/replay/repurpose_curated_bub1b
```

The additional public cache contains bibliographic metadata and small ChEMBL drug responses. The fallback reuses the existing reference cache. No patient genome, HPO combination or PGx table is sent externally. Named variants stay local because a gene is supplied. Their identity must match the variant-specific evidence bundle; that check does not confirm their function or phase. Complete copyrighted abstracts are not redistributed; the source ledger records access, hashes and links. For a fresh public run, omit `--offline` and use a writable cache.

## What the score means

The external Pro report supplied ordinal G/E/S/F/N judgments. The engine computes `6G + 5E + 4S + 3F + 2N`; it does not learn or discover these judgments. It enforces applicability caps and reports one-point sensitivity. The supplied inputs give arimoclomol 49, azithromycin 48 and niacin 43, with fragile ordering. A one-point change in G, E or S can reverse the co-leads. These are priorities for a discriminating experiment, not efficacy estimates or the hackathon judging rubric.

| Component | External judgment being scored | Interpretation retained by the engine |
|---|---|---|
| G | Conditional genotype/allele-class fit | Phase and allele function are unproved. The two stability estimates do not establish substantial N1002K destabilization. More mutant protein may not restore function. |
| E | Relevant efficacy evidence | Every proposed drug-to-MVA path has an untested bridge; E is capped at 1, including if a supplied score requests 5. Clinical efficacy in another disease does not remove this cap. |
| S | Pediatric/cancer suitability | This is an external caution judgment. An approved pediatric use in another indication does not pass MVA chronic-safety, oncology or renal gates. |
| F | Feasibility of the first test | Proposed assay feasibility is not a completed experiment or an established exposure bridge. |
| N | Novelty | Subjective research priority; no novelty point demonstrates efficacy. |

The original Pro documents' SHA-256 digests are stored in the input bundle. The result separately hashes the actual evidence bundle, PGx table and approval audit used. Current official-product approval comes from the independent verification audit. A missing, mismatched, stale or unverified official-product audit leaves the approval gate `UNKNOWN`. After seven days it must be refreshed; historical numerical priorities remain reproducible, while a stale approval must not be claimed current.

## Three distinct causal gaps

- **Arimoclomol:** analogous opsin proteostasis experiments and unstable BubR1 mutant precedent do not demonstrate arimoclomol rescue of useful N1002K BubR1. Context-dependent negative NPC1 results and HSF1 tumor-support biology remain visible.
- **Azithromycin:** analogous premature-stop reporter/protein signals do not establish endogenous Leu737Ter readthrough after NMD, a useful inserted residue, or functional rescue at feasible free exposure. Chronic antibacterial and organ-safety questions remain open.
- **Prescription nicotinic acid:** adult niacin-to-NAD data and separate NMN/SIRT2-to-BubR1 data leave the critical combined edge untested. K668/K250 regulation and RMS survival biology make direction uncertain. Niacin is not interchangeable with nicotinamide, NR or NMN.

The actual PGx input is applied to all resolved molecules. Explicitly unmatched rows remain in the JSON. No selected mapping for a lead means the current small marker table supplies no genotype-based advice for it; it does not mean normal metabolism or safety. No dose is selected.

## Outputs and gates

The JSON separates `curated_hypotheses` from automatic `ranking` and `supported_hypotheses`. The curated CSV labels priorities as external judgments. Grade H is an indirect laboratory hypothesis. The market-approved laboratory list requires an independent current official-product record, while every lead remains `clinical_eligible: false` because exact-genotype rescue, exposure, phase and clinical safety are unknown. Product approval for NPC, infection or lipid treatment is never MVA approval.

`results/repurpose_curated_sources.json` records the ten supporting/contrary paper records actually opened. Some publisher pages were inaccessible or incomplete; those cases rely only on the verified primary abstract and do not support figure-specific concentration claims. `results/repurpose_curated_bub1b.*` records the actual CLI run. The source ledger and the graph retain contrary evidence and missing bridges rather than counting multiple papers on separate steps as independent drug-rescue replications.

Clinical gate completeness is a supplied curator assessment, never software authorization to treat. Automatic target-action direction conflicts remain hard exclusions from the laboratory shortlist regardless of numerical priority. Case-specific scenario and selection language live in this evidence bundle; unrelated bundles receive generic unknown-applicability defaults.
