# Prospective low-pass analysis choices

Recorded 2026-10-24T01:59:15.671492+00:00, before the first chromosome-depth inference. This is a local analysis plan, not an external preregistration. The VCF screen and public-data results were already known. Downloads and paired subsampling had completed; the first whole-reference index attempt stopped at the 8-GB swap guard. No CN result had been accepted.

- Primary contrast: patient versus the public normal, matched on mapped CIGAR bases; 100-kb reference bins and 1-Mb spatial blocks.
- Primary family: all eligible whole-autosome and arm contrasts within that run. Bonferroni familywise threshold 0.05. Acrocentric short arms and insufficient blocks are excluded; actual tested count will be reported.
- Sensitivities: 50/250-kb bins, 5-Mb blocks, same-normal disjoint-template self-split, and patient versus half-normal at comparable mapped depth. These are controls/sensitivities, not 18 independent confirmatory runs or an across-run FWER guarantee.
- A block must retain at least half of its expected bins before equal block weighting. GC fits, bin masks and exact filter counts will be retained.
- Compare actual usable MAPQ30 autosomal depth as well as mapped CIGAR-base matching; the latter includes other contigs and is not identical to analysis depth.
- Self-split measures stochastic same-library sampling noise, not between-library, tissue or culture variation. Full/half contrasts reuse one biological normal and are not independent replication.
- The empirical technical null uses other-chromosome offsets and spatial blocks. At most 21 comparator offsets exist; resampling does not create new normal genomes. Confidence/detection limits remain conditional.
- Mosaic-fraction limits assume a single one-copy gain/loss against a mostly diploid baseline. Post-correction shift injection does not validate signal recovery through GC fitting. Opposing clones, genome-wide changes and tissue-restricted disease can be missed. No total MVA-cell fraction or clinical exclusion is inferred.
- No result selection by the most favorable bin/block combination. Compare primary findings with all sensitivities, control behavior and the VCF screen; discrepancies remain unresolved until explained.
- Restarted stages must use atomic outputs and verify complete primary-autosome coordinates plus input/parameter/code provenance. Partial files are not completed stages.

The final methods will identify the resource-compatible alignment command, software versions, realized coverage, controls, source checksums, runtime and limitations. This plan contains no patient-scale data.

## Pre-GC signal-recovery addendum

Recorded 2026-10-24T02:59:21.983942+00:00 while the normal alignment was still running, before any patient/control CN inference. After the primary analysis, inject one-copy gain/loss fractions 10%, 20% and 30% separately into each normal-selfsplit chromosome by multiplying raw covered bases by 1 ± f/2. Re-run the same masks and GC fit, and compare the target block statistic with the unspiked selfsplit baseline. Report recovery relative to the known log2 depth shift and retained-bin changes. This checks correction/filter attenuation at the depth-table level. It does not simulate reads, prove clinical sensitivity, supply independent biological controls, or automatically replace the conditional detection limits.

## Pairing-QC sensitivity addendum

Recorded 2026-10-24T03:15:34.577456+00:00 after normal alignment QC but before any patient depth inference. The normal has 15,975,720 usable Q30 primary records; 12,978,429 (81.2385%) have the proper-pair flag. The prespecified primary remains unchanged. Add one separate 100-kb/1-Mb contrast after requiring proper pairs in both samples and re-matching mapped CIGAR bases. This is a QC sensitivity outside the 18-run grid, not independent replication or an across-run familywise guarantee. Pairing filters can themselves change library/mappability selection.
