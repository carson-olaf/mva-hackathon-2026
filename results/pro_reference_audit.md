# Pro reference and mechanism audit

**All 51 title–PMID–DOI identities match original PubMed XML and Europe PMC metadata.** This verifies citations, not treatment efficacy. The two 2026 references are real records, and R07's full text was opened. No indexed retraction flags or links were found in the retrieved records; this is not an exhaustive research-integrity assessment.

## Access and verification scope

Actual PubMed HTML and DOI pages were attempted for every reference. All 51 PubMed HTML responses were status 203 without the expected citation metadata. Of 51 DOI attempts, 21 landing pages exposed matching citation metadata; 16 were JavaScript redirect pages and 14 returned 403. Original PubMed EFetch XML successfully supplied all 51 records, avoiding false verification from blank pages. Titles were compared after punctuation, capitalization and HTML-format normalization; DOI case was ignored.

The machine-readable ledger records every requested URL, HTTP result, source scope and cache hash. The original Pro opening grades are preserved separately from this audit. Other drug-claim and regulatory-label checks belong to the separate drug audit.

## Changes required when using the evidence

- R09 is the corrigendum for R08, not another allele-series experiment. It corrects duplicated Figure 2E sample images and the related pHH3 Methods description. Use the updated original article.
- R20 and R22 are also correction records. They must not add independent supporting lines.
- Additional corrections were found for R16 (PMID 32640244, DOI 10.1016/j.cmet.2020.05.020) and R38 (DOI 10.1007/s10577-014-9441-9). Their correction bodies were not verified here; the R38 publisher preview is paywalled. Do not invent their effect on conclusions.
- R05's broad KARD/PP2A mechanism is supported by its primary abstract, but full-text/figure access was not obtained. R06's relevant full primary-paper sections were read through its CC-BY mirror. Do not label both as newly audited journal full texts.
- The Pro 103-sequence conservation result, side-chain SASA and N1002→mouse N992 mapping remain unreplicated because its accompanying analysis files were unavailable. R08 independently confirms the different human L1012→mouse L1002 mapping.

## Selected mechanism claims

| Ref | Verified scope and finding | Required boundary |
|---|---|---|
| [R03](https://www.hubrecht.eu/app/uploads/2015/08/Kops_Research_2010_Suijkerbuijk_Molecular-Causes-for-BUBR1-Dysfunction-in-the-Human-Cancer-Predisposition-Syndrome-Mosaic-Variegated-Aneuploidy.pdf) | Full author manuscript, Figures 2–4: low protein, HSP90 dependence and forced-expression checkpoint rescue. | Other alleles/cultured cells; no N1002K or drug rescue. |
| [R05](https://pubmed.ncbi.nlm.nih.gov/23079597/) | Primary abstract: PLK1–KARD–PP2A-B56 regulates Aurora B-dependent attachment stability. | Full-text figures were inaccessible; exact-genotype effect untested. |
| [R06](https://www.researchgate.net/publication/347785117_BUBR1_Pseudokinase_Domain_Promotes_Kinetochore_PP2A-B56_Recruitment_Spindle_Checkpoint_Silencing_and_Chromosome_Alignment) | Full primary article at CC-BY mirror, Figures 3–4: normalized/titrated BubR1 can remain functionally defective. | Measure useful function separately from abundance. |
| [R07](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13385977/fullTextXML) | Full PMC XML and Figure 1B verify the Pro motif/domain map and TPR/MCC mechanism. | Annotation boundaries and engineered mutants do not establish either case allele's activity. |
| [R08](https://www.jci.org/articles/view/126863) | Publisher full text: allele/tissue-specific effects and human L1012→mouse L1002. | Mouse biology; apply R09 correction; no exact-genotype mTOR treatment evidence. |
| [R09](https://www.jci.org/articles/view/144781) | Publisher full corrigendum: Figure 2E images and related Methods corrected. | Correction provenance only; no independent biological experiment. |
| [R11](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3707109/fullTextXML) | Full PMC text: high WT transgene improves selected mouse fidelity, cancer and lifespan outcomes. | Genetic WT overexpression is not a mutant-rescue drug trial or human safety result. |

No selected source demonstrates rescue of N1002K, Leu737Ter or this genotype by any proposed drug. The mechanistic rationale remains a falsifiable experimental proposal.

## Complete identity ledger

Every row below matched both title and DOI in the original PubMed XML and Europe PMC record. ‘DOI page’ means a landing page also exposed matching citation metadata; it does not imply full-text access.

| Ref | PMID | DOI | Record note |
|---|---|---|---|
| R01 | [15475955](https://pubmed.ncbi.nlm.nih.gov/15475955/) | [10.1038/ng1449](https://doi.org/10.1038/ng1449) | Identity verified; DOI page |
| R02 | [9710612](https://pubmed.ncbi.nlm.nih.gov/9710612/) | [10.1128/mcb.18.9.5272](https://doi.org/10.1128/mcb.18.9.5272) | Identity verified |
| R03 | [20516114](https://pubmed.ncbi.nlm.nih.gov/20516114/) | [10.1158/0008-5472.can-09-4319](https://doi.org/10.1158/0008-5472.can-09-4319) | Identity verified |
| R04 | [22698286](https://pubmed.ncbi.nlm.nih.gov/22698286/) | [10.1016/j.devcel.2012.03.009](https://doi.org/10.1016/j.devcel.2012.03.009) | Identity verified |
| R05 | [23079597](https://pubmed.ncbi.nlm.nih.gov/23079597/) | [10.1016/j.devcel.2012.09.005](https://doi.org/10.1016/j.devcel.2012.09.005) | Identity verified |
| R06 | [33207204](https://pubmed.ncbi.nlm.nih.gov/33207204/) | [10.1016/j.celrep.2020.108397](https://doi.org/10.1016/j.celrep.2020.108397) | Identity verified |
| R07 | [42477419](https://pubmed.ncbi.nlm.nih.gov/42477419/) | [10.1038/s42003-026-10384-9](https://doi.org/10.1038/s42003-026-10384-9) | Identity verified; DOI page |
| R08 | [31738183](https://pubmed.ncbi.nlm.nih.gov/31738183/) | [10.1172/jci126863](https://doi.org/10.1172/jci126863) | Identity verified; DOI page; linked correction |
| R09 | [33136097](https://pubmed.ncbi.nlm.nih.gov/33136097/) | [10.1172/jci144781](https://doi.org/10.1172/jci144781) | Correction record; DOI page |
| R10 | [15208629](https://pubmed.ncbi.nlm.nih.gov/15208629/) | [10.1038/ng1382](https://doi.org/10.1038/ng1382) | Identity verified; DOI page |
| R11 | [23242215](https://pubmed.ncbi.nlm.nih.gov/23242215/) | [10.1038/ncb2643](https://doi.org/10.1038/ncb2643) | Identity verified; DOI page |
| R12 | [22048312](https://pubmed.ncbi.nlm.nih.gov/22048312/) | [10.1038/nature10600](https://doi.org/10.1038/nature10600) | Identity verified; DOI page |
| R13 | [24825348](https://pubmed.ncbi.nlm.nih.gov/24825348/) | [10.15252/embj.201386907](https://doi.org/10.15252/embj.201386907) | Identity verified; DOI page |
| R14 | [25285631](https://pubmed.ncbi.nlm.nih.gov/25285631/) | [10.1016/j.bbrc.2014.09.128](https://doi.org/10.1016/j.bbrc.2014.09.128) | Identity verified |
| R15 | [41533580](https://pubmed.ncbi.nlm.nih.gov/41533580/) | [10.1093/nar/gkaf1517](https://doi.org/10.1093/nar/gkaf1517) | Identity verified |
| R16 | [32386566](https://pubmed.ncbi.nlm.nih.gov/32386566/) | [10.1016/j.cmet.2020.04.008](https://doi.org/10.1016/j.cmet.2020.04.008) | Identity verified; linked correction |
| R17 | [24853414](https://pubmed.ncbi.nlm.nih.gov/24853414/) | [10.1038/cddis.2014.214](https://doi.org/10.1038/cddis.2014.214) | Identity verified; DOI page |
| R18 | [17889646](https://pubmed.ncbi.nlm.nih.gov/17889646/) | [10.1016/j.cell.2007.07.020](https://doi.org/10.1016/j.cell.2007.07.020) | Identity verified |
| R19 | [26620677](https://pubmed.ncbi.nlm.nih.gov/26620677/) | [10.1007/s00109-015-1364-1](https://doi.org/10.1007/s00109-015-1364-1) | Identity verified; DOI page; linked correction |
| R20 | [26897011](https://pubmed.ncbi.nlm.nih.gov/26897011/) | [10.1007/s00109-016-1394-3](https://doi.org/10.1007/s00109-016-1394-3) | Correction record; DOI page |
| R21 | [37276642](https://pubmed.ncbi.nlm.nih.gov/37276642/) | [10.1016/j.biopha.2023.114968](https://doi.org/10.1016/j.biopha.2023.114968) | Identity verified; linked correction |
| R22 | [37385919](https://pubmed.ncbi.nlm.nih.gov/37385919/) | [10.1016/j.biopha.2023.115082](https://doi.org/10.1016/j.biopha.2023.115082) | Correction record |
| R23 | [34066907](https://pubmed.ncbi.nlm.nih.gov/34066907/) | [10.3390/ijms22094923](https://doi.org/10.3390/ijms22094923) | Identity verified |
| R24 | [31382608](https://pubmed.ncbi.nlm.nih.gov/31382608/) | [10.3390/antibiotics8030110](https://doi.org/10.3390/antibiotics8030110) | Identity verified |
| R25 | [23680642](https://pubmed.ncbi.nlm.nih.gov/23680642/) | [10.1016/j.jcf.2013.04.006](https://doi.org/10.1016/j.jcf.2013.04.006) | Identity verified |
| R26 | [30266908](https://pubmed.ncbi.nlm.nih.gov/30266908/) | [10.1038/s41467-018-06507-4](https://doi.org/10.1038/s41467-018-06507-4) | Identity verified; DOI page |
| R27 | [37094517](https://pubmed.ncbi.nlm.nih.gov/37094517/) | [10.1016/j.redox.2023.102701](https://doi.org/10.1016/j.redox.2023.102701) | Identity verified |
| R28 | [23602569](https://pubmed.ncbi.nlm.nih.gov/23602569/) | [10.1016/j.celrep.2013.03.028](https://doi.org/10.1016/j.celrep.2013.03.028) | Identity verified |
| R29 | [18519516](https://pubmed.ncbi.nlm.nih.gov/18519516/) | [10.3324/haematol.12793](https://doi.org/10.3324/haematol.12793) | Identity verified; DOI page |
| R30 | [17616695](https://pubmed.ncbi.nlm.nih.gov/17616695/) | [10.1158/0008-5472.can-06-3012](https://doi.org/10.1158/0008-5472.can-06-3012) | Identity verified |
| R31 | [16138009](https://pubmed.ncbi.nlm.nih.gov/16138009/) | [10.4161/cc.4.10.2061](https://doi.org/10.4161/cc.4.10.2061) | Identity verified |
| R32 | [19855003](https://pubmed.ncbi.nlm.nih.gov/19855003/) | [10.1073/pnas.0904343106](https://doi.org/10.1073/pnas.0904343106) | Identity verified |
| R33 | [25341037](https://pubmed.ncbi.nlm.nih.gov/25341037/) | [10.1038/cddis.2014.385](https://doi.org/10.1038/cddis.2014.385) | Identity verified; DOI page |
| R34 | [28633018](https://pubmed.ncbi.nlm.nih.gov/28633018/) | [10.1016/j.devcel.2017.05.022](https://doi.org/10.1016/j.devcel.2017.05.022) | Identity verified |
| R35 | [26404941](https://pubmed.ncbi.nlm.nih.gov/26404941/) | [10.1101/gad.269118.115](https://doi.org/10.1101/gad.269118.115) | Identity verified |
| R36 | [28738408](https://pubmed.ncbi.nlm.nih.gov/28738408/) | [10.1038/nature23449](https://doi.org/10.1038/nature23449) | Identity verified; DOI page |
| R37 | [29342134](https://pubmed.ncbi.nlm.nih.gov/29342134/) | [10.1038/nature25432](https://doi.org/10.1038/nature25432) | Identity verified; DOI page |
| R38 | [24981203](https://pubmed.ncbi.nlm.nih.gov/24981203/) | [10.1007/s10577-014-9432-x](https://doi.org/10.1007/s10577-014-9432-x) | Identity verified; DOI page; linked correction |
| R39 | [42434306](https://pubmed.ncbi.nlm.nih.gov/42434306/) | [10.3389/fendo.2026.1838559](https://doi.org/10.3389/fendo.2026.1838559) | Identity verified; DOI page |
| R40 | [21552266](https://pubmed.ncbi.nlm.nih.gov/21552266/) | [10.1038/ng.822](https://doi.org/10.1038/ng.822) | Identity verified; DOI page |
| R41 | [28553959](https://pubmed.ncbi.nlm.nih.gov/28553959/) | [10.1038/ng.3883](https://doi.org/10.1038/ng.3883) | Identity verified; DOI page |
| R42 | [35044816](https://pubmed.ncbi.nlm.nih.gov/35044816/) | [10.1126/sciadv.abk0114](https://doi.org/10.1126/sciadv.abk0114) | Identity verified |
| R43 | [26621532](https://pubmed.ncbi.nlm.nih.gov/26621532/) | [10.1007/s00439-015-1619-5](https://doi.org/10.1007/s00439-015-1619-5) | Identity verified; DOI page |
| R44 | [29195078](https://pubmed.ncbi.nlm.nih.gov/29195078/) | [10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049) | Identity verified |
| R45 | [35524556](https://pubmed.ncbi.nlm.nih.gov/35524556/) | [10.1093/nar/gkac328](https://doi.org/10.1093/nar/gkac328) | Identity verified |
| R46 | [26771021](https://pubmed.ncbi.nlm.nih.gov/26771021/) | [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) | Identity verified |
| R47 | [34791371](https://pubmed.ncbi.nlm.nih.gov/34791371/) | [10.1093/nar/gkab1061](https://doi.org/10.1093/nar/gkab1061) | Identity verified |
| R48 | [31513481](https://pubmed.ncbi.nlm.nih.gov/31513481/) | [10.1200/jco.19.00576](https://doi.org/10.1200/jco.19.00576) | Identity verified |
| R49 | [38936378](https://pubmed.ncbi.nlm.nih.gov/38936378/) | [10.1016/s1470-2045(24)00255-9](https://doi.org/10.1016/s1470-2045(24)00255-9) | Identity verified |
| R50 | [39264246](https://pubmed.ncbi.nlm.nih.gov/39264246/) | [10.1158/1078-0432.ccr-24-1098](https://doi.org/10.1158/1078-0432.ccr-24-1098) | Identity verified |
| R51 | [34481829](https://pubmed.ncbi.nlm.nih.gov/34481829/) | [10.1016/j.jlr.2021.100114](https://doi.org/10.1016/j.jlr.2021.100114) | Identity verified |

## Reproduction and provenance

The public-only fetch/parse inputs are cached in ignored `work/pro-reference-cache/`. Run `uv run python work/pro-reference-build-audit.py` to regenerate this audit from those caches. The audit JSON retains input/source SHA-256 digests and the exact access distinctions. No case genotype table or private clinical details are present.
