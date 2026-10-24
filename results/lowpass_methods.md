# Low-pass chromosome and arm depth results

No eligible contrast in the prespecified primary analysis exceeded its nominal familywise threshold. An additional proper-pair sensitivity flagged chr8p; this small, method-sensitive contrast is unconfirmed. Neither outcome establishes or excludes mosaic aneuploidy or total MVA burden.

Primary: full patient versus public normal; 100-kb bins, 1-Mb blocks. The other 17 runs are sensitivities/controls, not independent confirmation or an across-18-run FWER claim. The most favorable sensitivity is not selected.

Computational methods and executable public-module commands: [low-pass methods](../report/lowpass_methods.md).

Validated execution summary SHA256: `17d3147a0bf7f2bf80d9aca5ee9d4a1c285962e48cb73277f219c7843ed7066d`.

## Realized input and mapping

- Patient: 11,412,959 selected pairs, 3,401,061,782 selected bases (1.101284x nominal input); 3,257,969,076 mapped CIGAR bases before matching; mapped fraction 0.986498, proper-pair fraction 0.845466.
- Normal: 11,256,953 selected pairs, 3,399,599,806 selected bases (1.100810x nominal input); 2,809,346,909 mapped CIGAR bases before matching; mapped fraction 0.972579, proper-pair fraction 0.709060.

| Contrast | Side | Selection fraction | Actual matched CIGAR bases | Nominal mapped coverage |
|---|---|---:|---:|---:|
| Primary | left | 0.862300053642 | 2,809,194,471 | 0.909634x |
| Primary | right | 1 | 2,809,346,909 | 0.909683x |
| Normal self-split | left | 0.999591478175 | 1,404,368,384 | 0.454743x |
| Normal self-split | right | 1 | 1,404,386,476 | 0.454749x |
| Patient/half-normal | left | 0.431056388578 | 1,404,075,127 | 0.454648x |
| Patient/half-normal | right | 1 | 1,404,368,384 | 0.454743x |

Usable primary MAPQ30 autosomal-bin mean depth: patient **0.916896x**, normal **0.801107x**. CIGAR matching includes other contigs and is not identical to this depth. The depth difference must be considered before interpreting conditional power.

## Primary whole-autosome results

| Chromosome | Relative CN | 95% empirical-null interval | Adjusted p | Conditional gain 90% limit | Conditional loss 90% limit |
|---|---:|---|---:|---:|---:|
| chr1 | 1.9974 | 1.9593–2.0256 | 1 | 7.29% | 7.98% |
| chr2 | 1.9994 | 1.9617–2.0268 | 1 | 7.07% | 7.77% |
| chr3 | 1.9920 | 1.9534–2.0209 | 1 | 7.45% | 8.12% |
| chr4 | 1.9845 | 1.9464–2.0132 | 1 | 7.38% | 8.07% |
| chr5 | 1.9978 | 1.9586–2.0267 | 1 | 7.55% | 8.24% |
| chr6 | 1.9970 | 1.9576–2.0264 | 1 | 7.72% | 8.41% |
| chr7 | 2.0055 | 1.9649–2.0349 | 1 | 7.88% | 8.63% |
| chr8 | 1.9955 | 1.9549–2.0261 | 1 | 8.00% | 8.61% |
| chr9 | 2.0001 | 1.9592–2.0314 | 1 | 8.28% | 8.91% |
| chr10 | 2.0174 | 1.9753–2.0469 | 1 | 8.33% | 9.17% |
| chr11 | 1.9935 | 1.9525–2.0241 | 1 | 8.06% | 8.72% |
| chr12 | 1.9816 | 1.9420–2.0127 | 1 | 7.80% | 8.38% |
| chr13 | 2.0179 | 1.9732–2.0499 | 1 | 8.59% | 9.32% |
| chr14 | 1.9743 | 1.9314–2.0039 | 1 | 8.50% | 9.29% |
| chr15 | 1.9810 | 1.9368–2.0143 | 1 | 9.09% | 9.75% |
| chr16 | 2.0388 | 1.9964–2.0725 | 1 | 8.57% | 9.00% |
| chr17 | 1.9969 | 1.9519–2.0314 | 1 | 9.25% | 9.77% |
| chr18 | 2.0257 | 1.9799–2.0602 | 1 | 9.49% | 9.97% |
| chr19 | 2.0270 | 1.9781–2.0647 | 1 | 9.99% | 10.45% |
| chr20 | 2.0161 | 1.9674–2.0526 | 1 | 9.77% | 10.30% |
| chr21 | 2.0291 | 1.9738–2.0740 | 1 | 12.40% | 12.48% |
| chr22 | 2.0206 | 1.9649–2.0659 | 1 | 11.72% | 11.96% |

Full arm results, within-region bootstrap intervals and individual GC-normalized patient/normal profiles remain in the corresponding CSV/JSON. The 95% intervals are pointwise, not simultaneous familywise bands. The smallest primary adjusted p-value is 0.3192 at chr8p (relative CN 2.0696); a pointwise interval excluding 2 does not establish a clone. Primary arm limits are 7.76–15.93% gain and 8.46–15.62% loss under the conditional single-copy model.

## Every sensitivity and control

| Run | Eligible tests | Nominal per-run outliers |
|---|---:|---|
| case_100kb_1mb | 61 | none |
| case_100kb_5mb | 59 | none |
| case_250kb_1mb | 61 | none |
| case_250kb_5mb | 59 | none |
| case_50kb_1mb | 61 | none |
| case_50kb_5mb | 59 | none |
| case_half_100kb_1mb | 61 | none |
| case_half_100kb_5mb | 59 | none |
| case_half_250kb_1mb | 61 | none |
| case_half_250kb_5mb | 59 | none |
| case_half_50kb_1mb | 61 | none |
| case_half_50kb_5mb | 59 | none |
| normal_selfsplit_100kb_1mb | 61 | none |
| normal_selfsplit_100kb_5mb | 59 | none |
| normal_selfsplit_250kb_1mb | 61 | none |
| normal_selfsplit_250kb_5mb | 59 | none |
| normal_selfsplit_50kb_1mb | 61 | none |
| normal_selfsplit_50kb_5mb | 59 | none |

The normal self-split controls stochastic same-library sampling noise. It does not control between-person, tissue, culture or library differences. All full/half comparisons reuse one biological normal and are not independent rescue evidence. GC fitting can absorb chromosome signal; the limits shift already-corrected statistics and assume signal preservation through fitting/filtering. They are not end-to-end validated sensitivity or a patient mosaicism upper bound. Opposing clones and tissue-restricted disease may be missed.

## Integrity and resource record

Original source files verified unchanged by SHA256, size and mtime: 3. Downloaded patient FASTQs removed after validated BAM QC: True. Downloaded normal FASTQs removed after validated BAM QC: True.

Among 704 recorded samples across attempts, minimum observed free disk was 33.817GB and maximum sampled swap was 6.040GB. These are sampled observations, not continuous extrema. The initial full-index attempt crossed the 8GB swap trigger and was stopped; the full-read partition attempt was stopped for projected scratch space. Neither attempt supplied an accepted BAM. Final alignment uses 750,000-pair batches with globally merged reference partitions and global duplicate removal.

## Pre-GC correction recovery

The separately planned control perturbation completed all 132 chromosome/direction/fraction cases (22 autosomes, gains/losses, 10/20/30% single-copy fractions). It multiplied covered-base depth before refiltering and GC fitting, then compared the corrected statistic with the unspiked self-split. Recovered signal ranged 88.4051–102.2996%, median 99.5673%; the minimum was chr19 loss at 10%. This tests correction/filter attenuation in one control, not read alignment, clinical sensitivity or smaller clone fractions. No conditional detection limit was automatically adjusted. [Recovery aggregates](lowpass_normal_selfsplit_pre_gc_recovery.csv) and [manifest](lowpass_normal_selfsplit_pre_gc_summary.json) retain all cases and input/code/output hashes.

## Additional proper-pair sensitivity and unresolved chr8p contrast

This check was planned before patient depth inference and is outside the 18-run grid. It requires proper pairs, re-matches mapped CIGAR bases, and uses the same 100-kb/1-Mb estimator. The realized mapped-base ratio is 1.000135; usable patient/normal depths are 0.740494×/0.677605×. Of 61 eligible tests, chr8p crosses the conditional within-run threshold. No across-run correction or independent biological replication is claimed.

| chr8p contrast | Relative CN | Pointwise 95% block-bootstrap interval | Within-run adjusted p |
|---|---:|---|---:|
| Prespecified primary | 2.0696 | 2.0389–2.0990 | 0.3192 |
| Additional proper-pair filter | 2.0875 | 2.0040–2.1169 | 0.03457 |
| Normal self-split, primary settings | 2.0276 | 1.9978–2.0593 | 1.0000 |

The proper-pair empirical-null 95% interval is 2.0322–2.1359, distinct from the within-region interval above. Across the six full-depth grid settings, chr8p relative CN is 2.0516–2.0814 with adjusted p≥0.1179. All six half-depth settings remain nonsignificant. These reuse the same patient and comparator; significance is method-sensitive. Pair filtering changes library/mapping selection, not only measurement noise, and a smaller p-value does not establish greater accuracy.

The existing strict VCF screen has 34,426 chr8p heterozygous SNVs, mean alternate fraction 0.5014 and adjusted p=1. Its conditional gain detection limit is 10.04%; absence of an allele-balance outlier does not exclude a smaller clone. There is no confirmed patient clone or defensible clone-fraction estimate. A matched uncultured reference and an orthogonal cell-level or targeted copy-number assay are appropriate follow-up, before assigning biological origin. The extra result is retained in [pairing provenance](lowpass_pairing_sensitivity.json) and [all regions](lowpass_proper_case_100kb_1mb_regions.csv).
