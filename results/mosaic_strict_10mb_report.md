# VCF-only mosaic aneuploidy screen

No chromosome or arm exceeded the conservative familywise empirical-null threshold. This does not exclude MVA or mosaic aneuploidy.

## Observed data

- 1,779,707 autosomal high-quality heterozygous SNVs; mean allele-depth 44.82; mean alternate allele fraction 0.50058.
- 287 eligible 10-Mb spatial blocks; 49 chromosome/arm tests; 30,000 bootstrap replicates; seed 20261024.
- Same-genome pseudo-controls: 1 / 21 whole chromosomes have unadjusted one-sided p<0.05. They are not known negative controls.
- Whole-chromosome 90%-power conditional gain fraction range: 8.8–10.3%; loss: 8.1–9.3%. Each region has its own limit in the CSV. These are model-based planning estimates, not patient-specific upper confidence bounds.

## Estimator and null

We use biallelic PASS SNVs with QUAL>=30, MQ>=40, GQ>=30, DP15–100, diploid heterozygous GT, and usable ref+alt depth15–100 and >=80% of DP. Both alleles must have at least one read. We do not filter on a narrow allele-balance range. Centromeres and 1-Mb flanks are excluded; blocks require >=500 SNVs; regions require >=5 blocks. X is descriptive only because sex/ploidy/PAR and tissue baselines are not established. Acrocentric short arms without sufficient markers are untestable.

Exact filter for this run: `FILTER="PASS" && TYPE="snp" && N_ALT=1 && QUAL>=50 && INFO/MQ>=60 && INFO/QD>=5 && INFO/FS<=20 && INFO/SOR<=2 && INFO/DB=1 && FMT/GQ>=50 && FMT/DP>=20 && FMT/DP<=80`. The optional strict sensitivity also requires a dbSNP flag, MQ60, QD5, FS<=20, SOR<=2, QUAL50, GQ50 and DP20–80. Genotype and allele-depth safeguards still apply.

For alternate count a of n reads, q=((a−n/2)²−n/4)/(n(n−1)) is unbiased for (p−0.5)² under independent binomial reads. We subtract a leave-one-chromosome, exact-depth-stratified autosomal baseline. The region statistic weights spatial blocks equally. The null resamples within-chromosome centered 10-Mb block residuals and adds a randomly sampled comparator chromosome offset. This includes a conservative floor for between-chromosome technical variation. One-sided p-values test excess broadening; Bonferroni controls the nominal 0.05 family across eligible chromosomes and arms. Empirical support is only 21 comparator chromosomes. Confidence intervals are observed effect minus the 2.5/97.5% empirical-null quantiles; calibration relies on exchangeability.

## Conditional detection model

For a simple clonal single-homolog gain in fraction f, the true unphased bands are 0.5 ± f/[2(2+f)]; for loss, 0.5 ± f/[2(2−f)]. Expected broadening is delta². We inject this expected broadening into each empirical-null statistic. The reported 90%-power threshold is the familywise upper null quantile minus the lower 10% null quantile. This is an aggregate signal injection, not a genotype re-calling experiment. It assumes the same measurement noise, retained heterozygotes, and one altered homolog throughout the region. Loss/gain allele identity cannot be resolved without phasing or depth. Opposite-homolog events in different cells, balanced gains/losses, heterogeneous cell clones, and tissue-restricted events can cancel. The genome baseline could itself contain diffuse MVA signal; the screen detects regional contrasts and cannot measure total aneuploid cell burden.

## Depth and limitations

Depth is available only at selected variant loci. The reported depth ratios describe the VCF and cannot serve as GC-corrected read-depth CN calls: missing invariant sites, heterozygote ascertainment, mapping quality, GC/mappability and depth truncation all bias them. Ref/alt bias and local overdispersion can mimic BAF broadening. No matched public normal was processed for this VCF-only scan. The conditional detection estimates are not validated clinical sensitivity and must not be used to rule out low-level, nonclonal, or tissue-specific MVA. Read-depth alignment, a matched processed normal, and/or single-cell karyotyping are the next measurements.

## Reproduction and public sources

`UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz`

Public GRCh38 arm boundaries: [UCSC cytoBand](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz); cached under work/, SHA-256 in JSON. The source was downloaded and parsed; no patient information was sent.
Methods context, opened: Jakubek et al., *Mosaic chromosomal alterations in blood across ancestries using whole-genome sequencing*. PMID [37904051](https://pubmed.ncbi.nlm.nih.gov/37904051/), DOI [10.1038/s41588-023-01553-1](https://doi.org/10.1038/s41588-023-01553-1). This exploratory implementation is not MoChA and does not claim its validated performance.
