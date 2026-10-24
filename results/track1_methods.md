# Phenotype-first discovery

No causal gene is supplied. Public HPO annotations rank genes by summed inverse prevalence of the eight terms; the top 40 genes are scanned locally.

QC: PASS, QUAL>=30, GQ>=30, DP>=10, biallelic calls; heterozygous balance 0.2–0.8. Canonical CDS references establish coding effect locally with reference matching and strand-aware translation. Score: HPO prior/3 + consequence (6 truncating; 2 missense) + inheritance compatibility (4 unphased recessive pair; 3 homozygous; 0 dominant). Up to 12 named shortlisted variants go to VEP; AF>1% excludes returned alleles.

- Phenotype-first shortlist, not exhaustive genome-wide consequence annotation.
- Canonical coding transcripts only; structural variants and most splice variants not assessed.
- Single-sample INFO/AF is not population frequency and is never used.
- Rarity confirmed only where VEP returns frequencies; missing frequency is unknown, not zero.
- Non-shortlisted variants use impact and phenotype heuristics without a rarity claim.
- Recessive pairs unphased; dominant models merely de novo compatible without parental data.
- All supplied terms included; one family-history term is withheld from the public copy.
- EPCR is subjective, not a statistical output.

Counts and NMD arithmetic: track1_evidence.json. NMD 50–55 nucleotide rule: PMID:35052820 (https://pubmed.ncbi.nlm.nih.gov/35052820/).
