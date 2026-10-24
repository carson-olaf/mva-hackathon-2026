# Independent lead-drug verification

All three leads satisfy the market-approved ingredient requirement through verified current US prescription products. **None has demonstrated rescue of this BUB1B genotype, and none is clinically eligible on this evidence.** This audit checks current official products, labels and specific scientific claims; it does not recommend a dose or a treatment.

## Official product and formulation checks

| Lead | Verified current US product | Label boundary and practical caveats |
|---|---|---|
| Arimoclomol | MIPLYFFA, NDA214927; approved 2024-09-20; oral capsules 47/62/93/124 mg base | NPC neurological manifestations **with miglustat**, age ≥2. Capsules can be opened per instructions. GI/weight effects, hypersensitivity and renal exposure matter; OCT2 inhibition can raise creatinine without reduced GFR. This does not establish MVA monotherapy or oncology safety. |
| Azithromycin | ZITHROMAX oral suspension, NDA050710; approved 1995-10-19; current Rx 100/200 mg per 5 mL; tablets separately current | Specific bacterial indications; selected pediatric indications ≥6 months, pharyngitis ≥2 years. QT/torsades, hepatic, hearing/GI and resistance concerns matter. Infection-course precedent does not establish chronic MVA safety. Ingredient initial approval1991 is distinct from this suspension application. |
| Prescription nicotinic acid | NIACIN ER, ANDA203899; approved 2017-06-16; current Rx 500/1000 mg tablets | Lipid indications. Safety/effectiveness and PK at age ≤16 are unestablished. ER tablets must remain whole; liver, glucose, urate/flushing and renal issues matter. Nicotinamide, NR, NMN and supplements are not interchangeable evidence for this prescription formulation. |

Official labels were opened: [MIPLYFFA March 2026](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=5feffc0e-453d-47fa-91dd-38d4952309bc), [ZITHROMAX July 2026 FDA label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/050710s054lbl.pdf), [niacinER](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=249e27ac-e1e7-42dc-84f6-cd32fd1d9dee). FDA product records, source hashes and retrieval times are in pro_drug_verification.json. Official regulatory records do not have PMID/DOI identifiers; none are fabricated. Market status does not establish UK licensing, local supply or MVA benefit.

## Exposure checks

The arimoclomol label reports NPC pediatric mean serum trough/peak 206/523 ng/mL (population SD 60/194). A calculated free-base molecular weight 313.78 g/mol converts these to 0.657/1.667 µM total; using approximate 10% plasma binding gives 0.591/1.500 µM unbound. **Pro's arithmetic is correct.** These are population exposure anchors using an approximate binding assumption, not this child's free tissue PK. The roughly 4-hour half-life makes a constant 24-hour bath an imperfect comparator. Glucuronidation does not establish UGT1A1 genotype actionability.

The azithromycin label's pediatric day 5 mean total-plasma peaks 0.216/0.383 µg/mL are confirmed. The 100 µg/mL reporter condition is approximately 463/261 times those means. Nominal overlap of a 0.1 µg/mL cell bath with total plasma is not proof of sustained free cytosolic exposure, readthrough or safety; tissue/lysosomal accumulation cannot supply the missing measurement. No pediatric niacin exposure bridge was established.

## Scientific claim checks

- **PMID [26620677](https://pubmed.ncbi.nlm.nih.gov/26620677/), DOI [10.1007/s00109-015-1364-1](https://doi.org/10.1007/s00109-015-1364-1) — public publisher captions and figures.** Fig 2B: APC(R1450X) reporter 100 microgram/mL for 48 hours. Fig 5A: RTT MeCP2(R294X) protein signal at 0.1 microgram/mL for 7 days. Nominal concentration ratio 1000. Different systems/endpoints. SMN1-deleted/SMN2-copy Fig 5C is not a defined nonsense-allele rescue assay; its low-dose labels differ 0.1/0.01. No BUB1B readthrough or functional rescue shown.

- **PMID [24853414](https://pubmed.ncbi.nlm.nih.gov/24853414/), DOI [10.1038/cddis.2014.214](https://doi.org/10.1038/cddis.2014.214) — fulltext methods results Fig1 discussion.** Arimoclomol 1 micromolar for 24 hours in mutant-opsin SK-N-SH cells increased Hsp70 and reduced aggregation; discussion explicitly reports no improved mutant trafficking. Analogous misfolded protein; lower aggregation can reflect disposal. Does not establish useful mutant BubR1 stabilization or HSP90-dependent rescue.

- **PMID [34481829](https://pubmed.ncbi.nlm.nih.gov/34481829/), DOI [10.1016/j.jlr.2021.100114](https://doi.org/10.1016/j.jlr.2021.100114) — fulltext results Fig7.** In NPC1-I1061T fibroblasts, arimoclomol did not correct cholesterol storage across tested exposures; body describes 1.6–500 micromolar for 3/5 days. HSP70 blots at 50/400 micromolar were negative. Cell/context-specific counterevidence. Does not refute NPC clinical approval or prove lack of effect in BUB1B cells.

- **PMID [32386566](https://pubmed.ncbi.nlm.nih.gov/32386566/), DOI [10.1016/j.cmet.2020.04.008](https://doi.org/10.1016/j.cmet.2020.04.008) — PubMed abstract only.** Adult-onset mitochondrial-myopathy niacin study reports increased blood/muscle NAD and muscle endpoints. This is a niacin-to-NAD bridge in another disease. No BUB1B/SIRT2 functional rescue or pediatric MVA regimen. Full formulation/methods not independently verified here. Linked erratum 32640244 exists; its text remains inaccessible and its effect is unverified.

- **PMID [24825348](https://pubmed.ncbi.nlm.nih.gov/24825348/), DOI [10.15252/embj.201386907](https://doi.org/10.15252/embj.201386907) — PubMed and PMC abstract only.** Abstract separates SIRT2 overexpression-mediated hypomorph lifespan benefit from NMN-mediated BubR1 abundance. K668 deacetylation is the proposed stability link. Neither niacin nor NMN lifespan treatment efficacy is established by this abstract. PMC API returned front/abstract without body; no figure-level quantitative lifespan claim independently verified.

- **PMID [25285631](https://pubmed.ncbi.nlm.nih.gov/25285631/), DOI [10.1016/j.bbrc.2014.09.128](https://doi.org/10.1016/j.bbrc.2014.09.128) — PubMed abstract only.** SIRT2 deacetylates BubR1 K250; knockdown increased acetylation but did not substantially raise total BubR1 or delay prometaphase-to-anaphase transition in that study. Do not assume K250 deacetylation changes total protein or that NAD/SIRT2 enhancement has an invariably favorable checkpoint effect.

- **PMID [41533580](https://pubmed.ncbi.nlm.nih.gov/41533580/), DOI [10.1093/nar/gkaf1517](https://doi.org/10.1093/nar/gkaf1517) — fulltext abstract results discussion.** Aurora B-dependent phosphorylation supports BubR1 K250acetylation, MCC/fibrous-corona maintenance and CENP-E-linked checkpoint control; acetylation-mimetic rescue is reported. This is a separate pathway constraint, not a niacin/SIRT2 drug experiment. Measure both K668andK250 and segregation function.

- **PMID [25341037](https://pubmed.ncbi.nlm.nih.gov/25341037/), DOI [10.1038/cddis.2014.385](https://doi.org/10.1038/cddis.2014.385) — fulltext results Fig7 Fig8.** SIRT1/SIRT2 knockdown impaired RD rhabdomyosarcoma growth; tenovin-6 also reduced growth in the reported xenograft experiment. Does not demonstrate niacin-induced recurrence or imply the reverse intervention necessarily promotes tumors. It warrants an RMS counter-screen.

The azithromycin quantitative distinction is independently supported by publicly accessible publisher [Figure2](https://link.springer.com/article/10.1007/s00109-015-1364-1/figures/2) captions and [Figure5](https://link.springer.com/article/10.1007/s00109-015-1364-1/figures/5) image/caption. The subscription article body was not accessed. Use **MeCP2 Figure 5A** for the exact 0.1 µg/mL low-dose claim; do not merge the SMN result into allele-specific nonsense rescue.

## Correction and verification gaps

Niacin has linked erratum PMID [32640244](https://pubmed.ncbi.nlm.nih.gov/32640244/), DOI [10.1016/j.cmet.2020.05.020](https://doi.org/10.1016/j.cmet.2020.05.020). Identity/link were verified, but publisher correction text returned403 and the PubMed record lacks an abstract. Its contents and effect remain unverified. Retain this warning rather than assume the correction is immaterial.

The Caspi correction PMID [26897011](https://pubmed.ncbi.nlm.nih.gov/26897011/), DOI [10.1007/s00109-016-1394-3](https://doi.org/10.1007/s00109-016-1394-3), was opened as a publisher preview/record. Author spellings differ, but the correction body is subscription content; this independent audit cannot certify that names are its only changes. The separate 2023 Apc-mouse corrigendum is outside this bounded figure audit and remains covered by the main reference ledger.

Four complete PMC article bodies were inspected (24853414,34481829,41533580,25341037). The PMC export for 24825348 contains front matter/abstract only; it was not silently counted as full-text review. No numerical lifespan effect or niacin formulation equivalence is independently claimed here.

## Pharmacogenomics

The current CPIC drug index returned no matches for the queried arimoclomol/azithromycin/niacin/nicotinic names. No actionable CPIC dose rule was established for these leads from the seven existing marker genotypes. The marker-level table is withheld from this public repository to protect the child's privacy. Absence of a guideline is not proof of absence of genetic variability or safety; organ function, interactions, formulation and measured exposure remain decisive gates.
