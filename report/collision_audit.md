# When target 701 became Gene 701: a public database collision census

24 October 2026. Public data only. No patient, genome or gated dataset was accessed for this audit.

**This was not two isolated errors. The December 2024 DGIdb export contains 13,146 distinct misassigned gene–GtoP-ligand pairs across 1,403 genes.** These collapse to **12,994 gene–normalized-drug pairs**. The strict count is 80.33% of that export's 16,365 Guide to Pharmacology rows. It is not 80.33% of DGIdb as a whole. A separate implementation reproduced the December count. [Computed census and source hashes](../results/collision_audit/summary.json), [independent checker](../scripts/independent_collision_check.py).

**The importer defect is already fixed upstream.** DGIdb merged PR624 on 14 April 2026. Both the June and October 2026 exports have zero confirmed instances of this exact numeric substitution. The October export and live API still contain **11 other target-identity conflicts**, whose cause is not established here. The actionable contribution is the historical scope, the rare-disease exposure map, and a reproducible list of residual conflicts—not a claim to have discovered an unfixed importer bug. [DGIdb PR624](https://github.com/dgidb/dgidb-v5/pull/624), [current conflicts](../results/collision_audit/current_identity_conflicts.csv), [live evidence](../results/collision_audit/live_api_evidence.json).

**What was counted**

An equal integer in two databases is only a collision opportunity. We count an observed erroneous association only when all these conditions hold:

1. A DGIdb GuideToPharmacology row carries a source claim `NCBIGENE:n`.
2. GtoP target `n` resolves to a different HGNC gene than the row's assigned gene.
3. The independent human NCBI Gene file maps Entrez `n` to the row's assigned HGNC gene.

Unresolved IDs and other disagreements are separate categories. The primary unit is a unique `(assigned HGNC gene, GtoP ligand ID)` pair. A second count uses DGIdb's normalized drug concept IDs; it is not a fresh chemical-structure normalization. All raw inputs, request bodies, dates and SHA-256 hashes are pinned in the [input manifest](../results/collision_audit/input_manifest.json). The census reads complete export subsets, not a sample, so no sampling confidence interval is appropriate.

| DGIdb snapshot | All source rows | GtoP rows | Exact numeric-collision rows | Distinct wrong source pairs | Other identity conflicts |
|---|---:|---:|---:|---:|---:|
| 2024-06 | 98,239 | 16,070 | 12,888 | 12,841 | 1,206 |
| 2024-12 | 98,920 | 16,365 | 13,146 | 13,146 | 1,158 |
| 2026-06 | 99,834 | 16,365 | 0 | 0 | 11 |
| 2026-10 | 122,904 | 17,105 | 0 | 0 | 11 |

Source: [recomputed snapshot summaries](../results/collision_audit/summary.json). The snapshots embed GtoP versions 2024.1, 2024.3, 2024.4 and 2026.2 respectively. The independent reference files are GtoP **2026.3.1**, published 12 October 2026. Historical GitHub release upload dates are not treated as original dataset dates.

Of the 13,146 December collisions, **13,112** still have the same human target–ligand pair in the current primary GtoP export; **9,196** also retain the same interaction-type wording. The type distinction matters: vocabulary changes can break a string match without changing target identity. The remaining pairs are retained in the historical count because the observed wrong-gene assignment and NCBI mapping meet the stated criterion. The 1,158 other December conflicts are excluded from that count. [Pair-level evidence](../results/collision_audit/identity_conflicts.csv).

**The cause is visible in code**

The old DGIdb importer constructed a gene claim as:

```ruby
create_gene_claim("NCBIGENE:#{line['Target ID']}", GeneNomenclature::NCBI_ID)
```

`Target ID` was a GtoP identifier. The importer had already built a target-to-Entrez dictionary, but this line did not use it. It also attached target aliases to that claim. The repaired importer creates a namespaced GtoP claim and retrieves it by the GtoP target ID. The patch also removes an alias path that labeled Ensembl IDs as NCBI IDs. [Opened upstream patch](https://github.com/dgidb/dgidb-v5/pull/624/files); its exact API response is included in the pinned public cache.

For example, `IUPHAR.RECEPTOR:2760` means **PDCD1**, whereas `NCBIGENE:2760` means **GM2A**. A valid PD-1 antibody association therefore became an invalid GM2A target association. This is a semantic identity failure; changing a drug-ranking weight cannot repair it.

**Named cases worth checking first**

These cases combine large misleading candidate sets or a substantial change in the apparent mechanism. “Worst” here is an audit priority, not a measured clinical-harm score. Every row is an observed December 2024 collision; every representative pair appears under its correct target in both tested 2026 exports. The full evidence, IDs, counts, current-source checks and opened pharmacology references are in [named_cases.csv](../results/collision_audit/named_cases.csv).

| Rare-disease search | Wrong historical lead | Actual target / colliding integer | Wrong source pairs for the gene | Opened pharmacology evidence |
|---|---|---|---:|---|
| GM2A; GM2 gangliosidosis AB variant | Nivolumab | PDCD1 / **2760** | 16 | PD-1 competition assays use nivolumab as comparator; PMID [32106752](https://pubmed.ncbi.nlm.nih.gov/32106752/), DOI 10.1080/19420862.2020.1724751 |
| GATM; primary Fanconi renotubular syndrome | Rituximab | MS4A1/CD20 / **2628** | 14 | Anti-CD20 binding and cytotoxicity study; PMID [15102696](https://pubmed.ncbi.nlm.nih.gov/15102696/), DOI 10.1158/1078-0432.ccr-03-0493 |
| ABCD1; X-linked cerebral adrenoleukodystrophy | (−)-Sulpiride | DRD2 / **215** | 78 | Dopamine D2/D3 receptor binding, including sulpiride; PMID [8301582](https://pubmed.ncbi.nlm.nih.gov/8301582/) |
| AARS1; Charcot–Marie–Tooth type 2N | Atropine | CHRM4 / **16** | 91 | Muscarinic receptor binding with atropine comparator; PMID [27085897](https://pubmed.ncbi.nlm.nih.gov/27085897/), DOI 10.1016/j.ejphar.2016.04.022 |
| ABCA4; Stargardt disease | Phenylephrine | ADRA1D / **24** | 49 | Human alpha-1 receptor functional assays; PMID [34355529](https://pubmed.ncbi.nlm.nih.gov/34355529/), DOI 10.1002/prp2.799 |
| CPT2; carnitine palmitoyltransferase II deficiency | Indomethacin | PTGS2/COX-2 / **1376** | 44 | COX-2 study with indomethacin comparator; PMID [29031075](https://pubmed.ncbi.nlm.nih.gov/29031075/), DOI 10.1016/j.ejmech.2017.09.056 |
| BUB1B; mosaic variegated aneuploidy | Niflumic acid | CLCNKB / **701** | 3 | Kidney chloride-current assay; PMID [16244177](https://pubmed.ncbi.nlm.nih.gov/16244177/), DOI 10.1124/mol.105.017384 |

The disease labels come from the assessed Orphanet gene associations, not from drug experiments. The GM2A and GATM disease links also have opened primary reports: PMID [8900233](https://europepmc.org/article/MED/8900233) and PMID [29654216](https://pubmed.ncbi.nlm.nih.gov/29654216/), DOI 10.1681/ASN.2017111179. The case-publication abstracts are preserved through Europe PMC; no unviewed full-text results are claimed. For nivolumab, the PubMed Figure 1 legend was also inspected. For niflumic acid, the pharmacology paper was traced separately because the current GtoP interaction's PMID field is empty.

An erroneous target assignment does not prove that the drug can never influence the wrong gene indirectly. It means **this source record cannot support that target mechanism**. No source found in this audit establishes the false direct-target bridges as rare-disease treatments.

**Which rare diseases are most exposed?**

Using Orphanet's 23 June 2026 gene map, we require an exact HGNC match, “Assessed” status and “Disease-causing germline mutation(s) in.” The December error set intersects **373 disease-associated genes and 512 ORPHA entities: 372 disorders plus 140 subtypes**. Of those entities, **207** inherit at least one ligand marked approved and not withdrawn in current GtoP. This measures misleading search leads, not affected patients, independent diseases or demonstrated harm. [Orphanet source](https://github.com/Orphanet/Orphadata_aggregated), [full ranking](../results/collision_audit/rare_disease_exposure.csv), [gene-association evidence](../results/collision_audit/rare_disease_gene_evidence.csv).

Ranking uses distinct normalized drug concepts with those flags; source-ligand counts remain available separately. Database approval flags are **not independently verified regulator approvals**.

| ORPHA entity | Affected gene(s) | Distinct approved-flag normalized drug concepts inherited from wrong links |
|---|---|---:|
| X-linked cerebral adrenoleukodystrophy / adrenomyeloneuropathy | ABCD1 | 41 each; overlapping gene and drug set |
| Trichothiodystrophy | AARS1, CARS1, ERCC2, ERCC3 | 29 |
| Marfanoid habitus–facial dysmorphism–skeletal abnormality–heart defect syndrome | ABL1 | 27 |
| Charcot–Marie–Tooth disease type 2N | AARS1 | 26 |
| Stargardt disease | ABCA4 | 22 |
| CPT II deficiency: myopathic, severe infantile or neonatal entity | CPT2 | 21 each; overlapping gene and drug set |
| GM2 gangliosidosis AB variant | GM2A | 12 |
| Primary Fanconi renotubular syndrome | GATM | 11 |

The unedited ranking also places Hirschsprung disease at 41 because the source labels ABCD1, ATP7A, ERBB2 and ERBB3 as assessed germline associations. We retain that result in the CSV, but do not equate those broad annotations with equally established monogenic causality. Retinitis pigmentosa and cone-rod dystrophy likewise aggregate several genes. This is why named single-gene examples are more useful than a bare “most exposed disease” leaderboard.

**Could a team encounter this today? Yes, through historical derivatives and current residual conflicts.**

A currently accessible 2025 epilepsy paper explicitly reports DGIdb-derived GM2A interactions with **nivolumab, cemiplimab and pembrolizumab** in Table 1. That is direct evidence of propagation beyond a database export. It does not demonstrate that those drugs were administered, nor invalidate every other analysis in the paper. [Opened full text, section 3.4 and Table 1; PMID 41165052, DOI 10.1002/epi4.70175](https://pmc.ncbi.nlm.nih.gov/articles/PMC12716302/).

The earlier [BUB1B audit](target_identity_audit.md) documented a second paper's checkpoint-gene/channel-drug mismatch. Its PMC page returned a browser challenge during this census, so that earlier inspection is not counted as a fresh full-text verification here. The independently reopened GM2A paper is sufficient to establish propagation. We did not survey every downstream website or count publications, users or patients affected.

The latest DGIdb export has **eight PDE10A associations assigned to LINC00473** and **three MCL1 associations assigned to ADAMTSL4-AS1**. All 11 were matched to live API responses on 24 October. Examples include papaverine and mardepodect in the first group, and venetoclax, AZD5991 and obatoclax in the second. These are identity conflicts, not proof that every listed ligand is a potent or selective inhibitor. Their source claims now use the correct GtoP namespace; they do not meet the demonstrated numeric-swap criterion. Root cause remains unknown. Other sources also contribute interactions to these normalized concepts, but those extra edges were not adjudicated. [Complete 11-row evidence](../results/collision_audit/current_identity_conflicts.csv), [saved live requests and responses](../results/collision_audit/live_api_evidence.json).

**What about our engine's actual sources?**

`mva_repurpose` currently uses Ensembl, Open Targets, Reactome and ChEMBL; it does not ingest DGIdb. We therefore tested the actual path separately. The complete ChEMBL 37 mechanism collection contains **7,561 records and 1,518 referenced targets**. Of the historical wrong pairs, **7,800 distinct wrong gene–ChEMBL-compound pairs** can be formed through GtoP ligand cross-references. **None overlaps the current ChEMBL mechanism collection**, checking both molecule and parent-molecule IDs against human target-component HGNC IDs. There are also **1,008 shared GtoP–ChEMBL protein identity checks with zero conflicts**. [Engine-source audit](../results/collision_audit/engine_source_audit.json).

HGNC target identities were available for **5,734 of 6,065 human mechanism records**; the other **331** are not certified by this HGNC comparison. As a positive control, the same join recovers **1,613 correctly assigned source counterpart gene–compound pairs** in ChEMBL. The zero wrong-pair result therefore does not come from a join that returns nothing.

For the 1,222 affected genes with one Ensembl cross-reference in the human NCBI file, Open Targets returned all 1,222 requested IDs with no ID or symbol discrepancy. A further 244 shared Open Targets–ChEMBL protein checks had no symbol discrepancy. The other 181 affected genes did not qualify for that request set. Reactome supplies pathway membership rather than direct drug–gene assertions; no exhaustive Reactome content audit is claimed. These bounded zeroes do not certify all database tables or all types of error.

**What maintainers can do now**

- Keep the upstream namespaced-ID fix. Reject any join that changes `IUPHAR.RECEPTOR:n` into `NCBIGENE:n` by retaining only the integer.
- Validate normalized identity after import: source human HGNC/Entrez/UniProt cross-references must agree with the normalized concept. Quarantine conflicts and unresolved cases before scoring.
- Investigate the MCL1 and PDE10A concept clusters. Rebuild aliases and derived associations from clean source assertions; editing a display label alone is insufficient.
- Publish affected-version and corrected-version guidance with a machine-readable invalidation list. Historical downloads, cached knowledge graphs and literature-derived candidate lists can outlive an importer repair.
- Run the sentinel cases plus the full identity reconciliation on each release. Preserve evidence types so a biomarker or cell-viability association cannot become a direct-target claim.

The [maintainer note](collision_maintainer_note.md) contains exact live queries, the existing fix, a reference validation rule and acceptance checks. It is a local draft; no issue or message was sent.

**Reproduce and falsify the census**

From the repository root, with the preserved public cache:

```sh
python3 scripts/replay_collision_audit.py
UV_CACHE_DIR=work/.uv-audit-cache uv run --offline pytest tests/test_collision_audit.py
```

The replay verifies **52 public input files** and reproduces **eight output files byte for byte**. It uses only Python's standard library. Explicit `--download` can fetch missing inputs using `requests`; any changed source hash stops reproduction. Mutable API responses and GtoP/NCBI download URLs may change later, so an exact public-only cache archive is also prepared locally. See [reproduction instructions](../results/collision_audit/README.md).

Six behavioral tests pass, including an end-to-end synthetic CLI run, a tampered-input rejection and removal of stale findings when a result is empty. A separate checker reproduces the main historical count and Orphanet join without importing the production audit. No genome tests or patient pipeline were run for this task. Total cached response bodies received: **171,629,134 bytes**, within the approved 250,000,000-byte limit.

The main limitation is temporal: the independent GtoP and NCBI identities are current, not archived 2024 reference files. The repaired exports, source code and named primary experiments corroborate the interpretation, but we do not claim to have re-executed the historical Rails importer. No clinical harm or global database error rate was measured. Source-specific licenses and attribution remain in [NOTICE.md](../results/collision_audit/NOTICE.md).

AI assistance: OpenAI; Pro; model training and data sharing disabled, as Carson confirmed. Computational findings were checked against saved public responses, a second implementation and explicit regression tests.
