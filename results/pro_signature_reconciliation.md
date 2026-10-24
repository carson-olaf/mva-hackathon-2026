# Reconciliation of the Pro and Codex LINCS diagnostics

**Conclusion: the two runs agree that transcriptional evidence earns no drug-ranking bonus. The Codex cross-modality/drug calculation is now classified `uninformative_quality_failure`, not a validated negative rescue result.** The Pro report prompted a provider-checksum audit that the initial Codex analysis had omitted. Local hashing and deterministic replay did not establish agreement with the provider's advertised source version. That omission is corrected in code and in the saved outputs.

## Keep the cohorts separate

| Item | Pro diagnostic, as reported | Codex diagnostic |
|---|---|---|
| Genetic perturbation | 15 BUB1B shRNA profiles | 21 BUB1B shRNA and 8 CRISPR profiles |
| Cells | A549, MCF7, PC3 | Those three plus A375 |
| Knockdown time | 96 h | 96 h |
| Drug data | 125 individual 2021 chemical profiles; selected drugs, cells, doses and 24-h exposure | 1,262 analyzed 2022 chemical-consensus profiles, pooled across contexts; phase-4/non-withdrawn filter |
| Genes | 959 common measured landmarks; 887 after E2F/G2M removal | 977 measured landmarks after excluding BUB1B |
| Reported issue | 125 chemical SHA mismatches; 67 chemical **file byte counts** mismatched | All 8 CRISPR files contradict provider hashes; drug-consensus metadata has **no usable advertised hashes** |

The Pro statement concerns **67 files with mismatched byte counts**, not 67-byte error payloads. The supplied Pro report and specification were read. Their companion `lincs/` manifests and raw files were not present among the supplied downloads, so the complete 125-file claim was not independently re-executed here. A small fresh named-drug audit below confirms the same class of upstream inconsistency. The cohorts are not pooled. Shared cell/condition subsets may overlap; agreement between the reports is not independent biological replication.

## Exact source-integrity audit

The cached public metadata and signature bytes were audited without any patient queries.

| Source | Files audited | Provider SHA/MD5 result | Provider size result | Meaning |
|---|---:|---|---|---|
| All cached 2021 shRNA profiles | 1,380 | 1,380 SHA and MD5 matches | 1,380 matches | Provider-verified bytes; this includes unrelated controls and exploratory profiles |
| Selected BUB1B 96-h shRNA profiles | 21 | All match | All match | Valid input to a shRNA-only reproducibility diagnostic |
| Selected BUB1B 96-h CRISPR profiles | 8 | **0/8 SHA matches; 0/8 MD5 matches** | 7 match; 1 differs by 16 bytes | Contradictory source provenance; blocks cross-modality inference |
| Analyzed 2022 chemical-consensus profiles | 1,262 | Metadata SHA is blank; no MD5 supplied | Metadata size fields are zero | Upstream checksum provenance unavailable; this is distinct from a contradictory hash |
| Fresh sample of 2021 chemical profiles | 14 | **0/14 SHA matches; 0/14 MD5 matches** | 10 match; 4 differ | Confirms the upstream problem described by Pro, in a separate small sample |

The default CLI's current selection contains 1,094 provider-verified shRNA vectors (21 BUB1B plus 1,073 selected unrelated controls), 8 contradictory CRISPR vectors, 1,262 drug vectors with no advertised checksums, and one unavailable drug file (HTTP 403). The larger 1,380-file audit includes earlier exploratory cached shRNA profiles and is not the inferential cohort size.

A concrete example is BUB1B CRISPR signature `dd01e09e-0af0-51bd-bc43-f98f029228c9` in A375 at 96 h:

- Provider SHA-256: `e639b727ac4964fdbd104449a124ad0731a86256961ac1912361936159ae2468`.
- Downloaded SHA-256: `83927be539e6b6d070c5edc934a1585ea4ec7881a0d21d0516ad9dbc5a746af1`.
- Both advertised and observed byte counts: **217,477**. Matching length does not repair the hash conflict.
- [Public source object](https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/cd/xpr/L1000_LINCS_DCIC_XPR015_A375.311_96H_N22_BUB1B.tsv).

Per-file advertised/actual hashes and sizes are retained in `signature_BUB1B_provider_hash_audit.csv`. Chemical-consensus audit rows are in `signature_analyzed_chemical_hash_audit.csv`; absent advertised values mean unavailable verification, never a pass.

## What the files actually contain

All 1,262 analyzed chemical-consensus objects are valid gzip files and parse as one numeric vector with **12,327 unique gene symbols**, finite values, and a drug column whose normalized name matches the associated metadata. No two analyzed drug vectors are identical under a gene-sorted numeric serialization. Every cached file's SHA and byte count agrees with its local download sidecar. These checks exclude obvious HTML/error-body contamination, malformed gzip payloads, label/header mismatch and duplicate-vector substitution. They do not prove the biological identity or correct generation of every vector.

A deterministic, dispersed sample of nine chemical-consensus files was fetched again, along with one provider-verified shRNA file and one contradictory CRISPR file. **All 11 fresh HTTP responses were status 200 and byte-identical to the cached files. All 11 S3 ETags matched the actual payload MD5.** The nine chemical responses had `binary/octet-stream` content type and valid gzip contents; the other two were tab-separated text. See `signature_transport_refetch_audit.json`.

Thus the cache and network transport are reproducible. The repeated fetch does **not** resolve a disagreement between the source object's current bytes and a separate provider manifest. It also does not manufacture a missing advertised checksum for the 2022 library.

For a fresh sample of two profiles each for azithromycin, everolimus, hydroxychloroquine, metformin, niacin, phenylbutyrate and sirolimus, all 14 public 2021 objects returned HTTP 200 and proper TSV content types. Their observed payload MD5s matched S3 ETags, while their advertised metadata hashes did not. These records are in `signature_2021_chemical_provider_hash_audit.json`. File lengths were approximately 217 KB, not tiny error pages.

Possible whitespace, newline, header, row-order and float-format explanations were tested on one contradictory CRISPR object: **356 plausible serializations, no advertised SHA match**. This is a diagnostic search, not proof that every possible transformation was excluded. The official retrieval tutorial also documents alternate `LINCS-data-2020/L1000/...` objects. Three small alternate objects were fetched; they contain multi-replicate expression matrices rather than the one-column characteristic-direction signatures, and their decompressed hashes do not equal the advertised signature hashes. They must not be silently substituted. See `signature_alternate_2020_object_audit.json`.

**The cause remains unresolved.** Stale metadata or regenerated objects are plausible explanations. Neither was proved, and numerical equivalence to the intended source version was not established. The audit does not justify calling every public table corrupt, nor does it permit treating failed hashes as verified data.

## The shRNA-only result that remains valid

The 21 selected BUB1B shRNA files and all shRNA null profiles used here pass provider SHA-256 and byte-count checks. Within each cell, we calculated the median pairwise Spearman correlation on the 977 measured non-target landmarks. For a cell/time matched null, we drew equally sized sets of complete unrelated shRNA profile vectors without replacement, with **999 resamples** and seed **20261024**. Whole-profile sampling retains covariance between genes. This diagnostic uses no chemical or CRISPR vectors.

| Cell | BUB1B profiles | Median pairwise rho | Competitive null p | BH q across 4 cells |
|---|---:|---:|---:|---:|
| A375 | 6 | 0.033119 | 0.049 | 0.098 |
| A549 | 6 | 0.006467 | 0.379 | 0.379 |
| MCF7 | 3 | 0.067042 | 0.042 | 0.098 |
| PC3 | 6 | 0.024573 | 0.077 | 0.102667 |

These are weak correlations, with no cell passing the four-test BH threshold of 0.05. The null is competitive reproducibility relative to the cached same-cell unrelated shRNA profiles; it is not a treatment-efficacy null. Batch and off-target effects remain possible, and repeated profiles from one cell line are not independent human samples. This finding concerns this assay's reproducibility, not whether BUB1B biology exists or any approved drug can work.

The differing Pro medians (0.008, 0.059, 0.025) are close in scale but must remain attributed to its 959-gene analysis. We used 977 non-target landmarks. Without the Pro companion manifests, exact signature and gene-set identity cannot be confirmed. No numerical values are averaged between runs.

## Effect on the earlier drug screen

The historical 4,999-resample calculation and its raw scores remain available for audit, explicitly labeled unvalidated. Its reported cross-modality rho −0.03377, zero adjusted discoveries, and best raw minimum reversal rho 0.07749 are **not validated biological conclusions**. They cannot establish that drugs fail to reverse BUB1B loss or that indoramin is preferable to the mechanism-led candidates.

All rows in `signature_ranked_approved_drugs.csv` and the diagnostics now carry `evidence_status = uninformative_quality_failure`, `ranking_bonus = 0`, and a failed source-gate flag. The figure is labeled as a historical unvalidated diagnostic. `signature_reversal.md` states the correction prominently.

The code now distinguishes `verified`, `mismatch`, `not_advertised` and `unavailable` provenance. It writes `signature_source_qc.json` and **stops before cross-modality correlation or drug permutations** when the source gate fails. A small HTTP re-fetch sample supports local reproducibility of the 2022 chemical library, but does not by itself justify validated drug ranking. No new drug-null calculation was run after the gate failed.

A verified shRNA-only diagnostic remains runnable:

```bash
UV_CACHE_DIR=work/.uv-cache uv run python -m mva_hackathon_2026.signature \
  --offline --shrna-qc --cache work/signature --output results
```

The ordinary offline command now reproduces the failed source gate, rather than repeating invalid inference. Resolving the conflict requires an authoritative manifest/object reconciliation or regeneration from an independently versioned source with documented transformations. Even then, the weak shRNA reproducibility must be addressed before a reversal score becomes useful.

## Sources and audit scope

The Pro documents read were `BUB1B_MVA_Track2.md` and `PIPELINE_SPEC.md`; their reported 125-file result remains a separately attributed result pending the companion bundle. Public-source methods were checked against the [official SigCom retrieval tutorial](https://maayanlab.github.io/dex-benchmark/L1000_data/Retrieving_L1000_Data.html) and SigCom's paper: **PMID 35524556; DOI 10.1093/nar/gkac328**, [opened primary full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9252724/fullTextXML). The characteristic-direction method source was also opened: **PMID 24650281; DOI 10.1186/1471-2105-15-79**, [primary full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4000056/fullTextXML).

All audit requests used public gene/drug names or public signature IDs. No patient data were queried or transmitted, and no large dataset was downloaded. Exact public URLs, returned byte hashes and HTTP metadata accompany the audit artifacts.
