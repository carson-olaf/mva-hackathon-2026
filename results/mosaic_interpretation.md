# Mosaicism: integrated interpretation

The VCF does not establish a chromosome- or arm-level aneuploid clone. The high-quality refinement analyzed **1,779,707 autosomal heterozygous SNVs**, mean AD depth **44.82**, with **60 eligible chromosome/arm tests**. No region passed nominal familywise 0.05 correction. This is a negative exploratory screen with a conditional sensitivity estimate; it is not evidence against the molecular MVA diagnosis.

## Quality-control finding

The broad PASS/MQ40/GQ30 set showed alternate-allele depletion on chr22 (mean BAF **0.47455**, excess variance **0.004057**). Requiring dbSNP, MQ60, QD5, FS<=20, SOR<=2, QUAL50, GQ50 and DP20–80 gave mean BAF **0.49997** and excess variance **0.00000932** on chr22. The disappearance of the deviation under stringent quality filtering supports a technical/ascertainment explanation rather than a robust aneuploid signal. This refinement was prompted by observed bias; it was not a prospectively registered filter. Both result sets are retained.

## Conditional detection limits

The limits below are the injected expected broadening needed for 90% power under the observed spatial-block/chromosome null and the single-homolog clonal gain/loss model. They are not clinical sensitivity and are not upper confidence bounds on this child's total aneuploid cell fraction.

| Read-quality set | Spatial blocks | Eligible tests | Whole-chromosome gain | Whole-chromosome loss | Arm gain | Arm loss |
|---|---:|---:|---|---|---|---|
| Broad PASS | 5 Mb | 60 | 37.8–47.1% | 27.4–32.0% | 39.5–56.6% | 28.3–36.1% |
| Broad PASS | 10 Mb | 51 | 47.7–56.2% | 32.3–36.0% | 49.8–60.6% | 33.2–37.7% |
| Stringent | 5 Mb | 60 | 7.8–10.3% | 7.3–9.3% | 8.2–11.5% | 7.6–10.3% |
| Stringent | 10 Mb | 49 | 8.8–10.3% | 8.1–9.3% | 9.0–10.4% | 8.3–9.4% |

The stringent 10-Mb sensitivity run also has no familywise outliers. It tests 21 whole chromosomes and 28 arms; chr21 and several small arms have too few independent blocks for that coarse resolution. The stringent 5-Mb run tests all 22 autosomes and 38 arms. All X results are descriptive; the VCF lacks an established sex/PAR/ploidy baseline for this analysis. Acrocentric short arms and chr18p are untestable with the 5-Mb block criterion.

The difference between broad and stringent conditional limits demonstrates why a pure binomial read-count limit would be misleading: low-quality marker classes can dominate regional background variance. The remaining stringent baseline is still only one genome, with 21 comparator chromosomes, and no independently processed normal.

## What was verified

- bcftools streamed the actual VCF; only aggregate depth/allele-count histograms were stored under ignored work/.
- Depth-corrected excess allelic variance is calibrated against leave-one-chromosome empirical spatial blocks plus a between-chromosome technical-variance floor.
- The null has 30,000 replicates per eligible region, fixed seed 20261024, one-sided tests and Bonferroni correction.
- The binomial-noise correction and inversion of the gain/loss model passed two behavioral tests, including simulated reads with a known allelic shift.
- All four complete scans ran; the high-quality result survives changing block size from 5 to 10 Mb. Figures were generated and visually checked.

## What was not verified

- No BAM, GC-corrected depth, invariant sites, matched normal or independent cytogenetics were available.
- The sensitivity calculation injects expected aggregate delta-squared into null statistics; it does not simulate a full read-generation and genotype-recalling pipeline. Genotype dropout, cell-type composition and tissue-specific sensitivity remain unvalidated.
- The tissue source is not established by the VCF. Even if this is blood, blood may miss disease-relevant clones in other tissues.
- Unphased broadening cannot distinguish gain, loss and copy-neutral allelic imbalance. Depth at called variant sites is ascertainment-biased and cannot resolve copy number.
- MVA contains multiple unstable clones. Gains/losses in opposing homologs or opposing directions can cancel in bulk sequencing. A common genome-wide burden also enters the within-genome baseline and is invisible to a regional-contrast test. **No total MVA burden is estimated.**

## Exact reruns

```sh
UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz
UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --block-mb 10 --prefix results/mosaic_10mb
UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --qc strict --prefix results/mosaic_strict
UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --qc strict --block-mb 10 --prefix results/mosaic_strict_10mb
uv run pytest -q tests/test_mosaic.py
```

Use `--reuse-cache` only when the private aggregate histograms exist and their path/size/mtime/filter/cytoband identities match. Public arm annotation is downloaded automatically if absent; no patient data are sent. Per-run JSON files retain annotation SHA-256, filters and QC counts; region CSV files contain aggregates only.

References opened: Jakubek et al., *Mosaic chromosomal alterations in blood across ancestries using whole-genome sequencing*, PMID [37904051](https://pubmed.ncbi.nlm.nih.gov/37904051/), DOI [10.1038/s41588-023-01553-1](https://doi.org/10.1038/s41588-023-01553-1). This work uses a separate exploratory estimator and does not claim MoChA performance. GRCh38 boundaries came from the downloaded [UCSC cytoBand table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz).
