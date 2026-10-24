# Uridine triacetate: one lead, one decisive rescue/harm test

**We propose Xuriden to test whether supplying uridine reduces replication-associated chromosome errors in BUB1B-MVA.** The FDA-approved oral product has chronic pediatric precedent in hereditary orotic aciduria. Its mean plasma uridine peak is **88.7 µM (n=4)**: **2.96×** the **30 µM nominal comparator**, not an established uridine-only effective concentration or an efficacy margin. Burrell's **45–55%** error reduction belongs to a four-nucleoside mixture. The retrieved supplement contains no component-specific arm. Uracil-DNA damage reported at **10 or 30 µM** makes benefit versus harm the decisive question. [Full report and sources](../report/track2_report.md).

**Two unproved transfers:** BUB1B errors must include a nucleotide-responsive component, and uridine alone must reproduce enough mixture benefit. Genetic correction establishes BUB1B dependence, not the correctable mechanism.

**Reject progression if uridine alone fails at measured, PK-informed exposure, or increases genomic damage or RMS fitness.** Compare exact-genotype and corrected cells, uridine alone and the mixture, with RMS counter-screens. Proposed laboratory gates: **≥30% fewer errors per completed division**, 95% interval excluding no benefit, **≥90% viable proliferative output**, and no reproducible harm. These are design thresholds, not expected effects or clinical eligibility. [Protocol](../report/nucleoside_decision.md).

**Track 1:** phenotype-first discovery ranks BUB1B Leu737Ter/Asn1002Lys first. NMD predicted from a **746-nt** stop-to-final-junction gap; phase unproven. CSV uses VCF sample name **WGS_EX2312012**; production scorer acceptance untested. [CSV](track1_carcango_hpo_recessive.csv).

| Routes deprioritized under available evidence | Verified result / boundary |
|---|---|
| Stability | DynaMut2 −0.051 / KORPM −0.573 kcal/mol; RMSE 1.14/1.33. Benchmark RMSE is not a variant CI; no convincing large effect, not proven neutrality or exclusion of stabilization rescue. |
| Pockets | 30 pockets; maximum 0.199, nearest 0.001/0.005. No modeled pocket passed; other conformations remain untested. |
| Expression | 8/8 CRISPR hash failures; verified shRNA rho 0.00647–0.06704, all four BH q≥0.098. Uninformative for efficacy. |
| Bulk aneuploidy | 1,779,707 SNVs, 60 tests, zero outliers. Low-pass 0.917×/0.801×, 61 primary tests, zero outliers. Conditional gain/loss limits 7.07–12.40% / 7.77–12.48%; modeled clonal events only, not total MVA burden. 18 runs are robustness checks, not biological replications. chr8p unconfirmed. |
| Modifiers / secondary | 49 genes; all 31 noncausal coding calls AF>0.001. Zero exact qualifying matches across 84 ACMG SF genes. Not negative clinical tests. |
| PGx | 36 markers; seven observed genotypes, no full diplotypes or actionable earlier-candidate dose rule. No uridine dose follows. |
| Mouse allele | Mouse L1002P is human L1012P, not N1002K. No exact-genotype in-vivo model verified. |
| Earlier shortlist | Arimoclomol, azithromycin and niacin have no independent positive BUB1B rescue line. Historical hypotheses only. |

**Methods finding:** **13,146** historical wrong gene–ligand pairs across **1,403** genes. The maintainers fixed the exact numeric importer bug; **11 unresolved discrepancies** remain. This census improves evidence checking, not uridine efficacy. [Census](../report/collision_audit.md).

**Contribution:** a potentially novel repurposing hypothesis, an audit of supporting/opposing evidence, and a mechanism-testing experiment. **Scalability:** public API engine, cached replay, explicit unknowns and TSC2–everolimus retrieval control. Automatic BUB1B output supplies context, not rescue evidence; uridine is a curated literature hypothesis. [Engine output and limits](../examples/README.md).

**Verified:** source identities/access scope, actual supplement, label arithmetic and byte-identical replay; **28 tests pass**. **Unverified:** uridine-only rescue, continuous effective exposure, patient PK, long-term cancer safety, phase/RNA/protein, total disease-cell burden and exact dataset citation. No clinician signoff, patient experiment, push or submission. AI: **OpenAI Pro; training/data sharing disabled**, as reported by Carson.
