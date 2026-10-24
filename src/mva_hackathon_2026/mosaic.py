"""Conservative, unphased chromosome/arm screen from a single-sample VCF.

Only aggregate statistics leave work/. This is a research screen, not a CNV caller.
Run with python -m mva_hackathon_2026.mosaic --help.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

AUTOSOMES = tuple(str(i) for i in range(1, 23))
FILTER = ('FILTER="PASS" && TYPE="snp" && N_ALT=1 && QUAL>=30 '
          '&& INFO/MQ>=40 && FMT/GQ>=30 && FMT/DP>=15 && FMT/DP<=100')
STRICT_FILTER = ('FILTER="PASS" && TYPE="snp" && N_ALT=1 && QUAL>=50 '
                 '&& INFO/MQ>=60 && INFO/QD>=5 && INFO/FS<=20 && INFO/SOR<=2 '
                 '&& INFO/DB=1 && FMT/GQ>=50 && FMT/DP>=20 && FMT/DP<=80')
CYTO_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"


def excess_moment(alt, depth):
    """Unbiased estimate of (true allelic probability - 0.5)^2 for binomial reads."""
    return ((alt - depth / 2) ** 2 - depth / 4) / (depth * (depth - 1))


def fraction_from_broadening(broadening: float, event: str) -> float | None:
    """Invert delta²; one haplotype gained/lost in fraction f of sampled cells."""
    delta = math.sqrt(max(0.0, broadening))
    denominator = 1 - 2 * delta if event == "gain" else 1 + 2 * delta
    fraction = 4 * delta / denominator if denominator > 0 else math.inf
    return fraction if fraction <= 1 else None


def broadening_from_fraction(fraction: float, event: str) -> float:
    denominator = 2 * (2 + fraction if event == "gain" else 2 - fraction)
    return (fraction / denominator) ** 2


def boundaries(path: Path) -> dict[str, tuple[int, int, int]]:
    """Return acen interval and chromosome end; UCSC coordinates are zero-based."""
    rows = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, start, end, _, stain = line.strip().split("\t")
            chrom = chrom.removeprefix("chr")
            if chrom not in (*AUTOSOMES, "X"):
                continue
            row = rows.setdefault(chrom, [10**10, 0, 0])
            row[2] = max(row[2], int(end))
            if stain == "acen":
                row[0], row[1] = min(row[0], int(start)), max(row[1], int(end))
    return {chrom: tuple(value) for chrom, value in rows.items()}


def stream_blocks(vcf: Path, cyto: dict, block_mb: int = 5, filter_expression=FILTER):
    """Collapse calls into depth/AD histograms without storing per-variant records."""
    command = ["bcftools", "query", "-i", filter_expression, "-f",
               "%CHROM\t%POS[\t%GT\t%AD\t%DP]\n", str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    blocks, qc = {}, {"queried_calls": 0, "usable_het_snvs": 0,
                      "excluded_nonprimary": 0, "excluded_nonhet": 0,
                      "excluded_ad": 0, "excluded_centromere_guard": 0}
    assert process.stdout is not None
    for line in process.stdout:
        qc["queried_calls"] += 1
        chrom, position, gt, ad_text, dp_text = line.rstrip().split("\t")
        if chrom not in cyto:
            qc["excluded_nonprimary"] += 1
            continue
        if gt.replace("|", "/") not in ("0/1", "1/0"):
            qc["excluded_nonhet"] += 1
            continue
        try:
            ref, alt = map(int, ad_text.split(","))
            depth, dp = ref + alt, int(dp_text)
        except ValueError:
            qc["excluded_ad"] += 1
            continue
        if not (15 <= depth <= 100 and depth >= .8 * dp and ref > 0 and alt > 0):
            qc["excluded_ad"] += 1
            continue
        position = int(position) - 1
        cent_start, cent_end, _ = cyto[chrom]
        if cent_start - 1_000_000 <= position < cent_end + 1_000_000:
            qc["excluded_centromere_guard"] += 1
            continue
        arm = "p" if position < cent_start else "q"
        key = (chrom, arm, position // (block_mb * 1_000_000))
        if key not in blocks:
            blocks[key] = np.zeros((101, 101), dtype=np.int32)
        blocks[key][depth, alt] += 1
        qc["usable_het_snvs"] += 1
    stderr = process.stderr.read() if process.stderr else ""
    if process.wait():
        raise RuntimeError(f"bcftools failed: {stderr[:1000]}")
    keys = sorted(blocks, key=lambda key: (23 if key[0] == "X" else int(key[0]), key[2]))
    return keys, np.stack([blocks[key] for key in keys]), qc, command


def histogram_statistics(histograms: np.ndarray):
    n, a = np.indices((101, 101))
    moment = np.zeros((101, 101))
    valid = n > 1
    moment[valid] = excess_moment(a[valid], n[valid])
    af = np.divide(a, n, out=np.zeros_like(a, dtype=float), where=n > 0)
    count = histograms.sum(axis=(1, 2))
    qsum = (histograms * moment).sum(axis=(1, 2))
    afsum = (histograms * af).sum(axis=(1, 2))
    depth_sum = (histograms * n).sum(axis=(1, 2))
    return count, qsum, afsum, depth_sum, moment


def null_distribution(block_means, chromosomes, target_chrom, n_blocks, rng, iterations):
    """Resample 5-Mb residual blocks AND empirical between-chromosome offsets.

    The chromosome offset is a conservative floor against long-range technical
    heterogeneity. Its empirical support is only 21 comparator chromosomes.
    """
    pool = chromosomes != target_chrom
    ids = np.unique(chromosomes[pool])
    offsets = np.array([block_means[chromosomes == chrom].mean() for chrom in ids])
    center = block_means[pool].mean()
    centered_blocks = np.concatenate([
        block_means[chromosomes == chrom] - block_means[chromosomes == chrom].mean()
        for chrom in ids
    ])
    draws = rng.choice(offsets - center, iterations)
    # Region means weight spatial blocks equally, so dense SNP regions cannot dominate.
    for start in range(0, iterations, 1000):
        stop = min(iterations, start + 1000)
        draws[start:stop] += rng.choice(centered_blocks, (stop - start, n_blocks)).mean(axis=1)
    return draws, center


def analyze(keys, histograms, iterations=30000, seed=20261024):
    rng = np.random.default_rng(seed)
    count, qsum, afsum, depth_sum, moment = histogram_statistics(histograms)
    chroms = np.array([key[0] for key in keys])
    arms = np.array([key[1] for key in keys])
    autosomal = np.isin(chroms, AUTOSOMES)
    # Exclude tiny sparse blocks to make block exchangeability less implausible.
    usable = (count >= 500) & autosomal
    if usable.sum() < 100:
        raise ValueError("Too few autosomal blocks for empirical calibration")
    depth_hist = histograms.sum(axis=2)
    # Use a leave-one-chromosome depth baseline to avoid attenuating a target signal.
    regions, control = [], []
    for chrom in (*AUTOSOMES, "X"):
        for arm in ("", "p", "q"):
            selected_all = (chroms == chrom) & ((arms == arm) if arm else True)
            selected = selected_all & (count >= 500)
            row = {"region": f"chr{chrom}{arm}", "chromosome": chrom, "arm": arm or "whole",
                   "heterozygous_snvs": int(count[selected_all].sum()),
                   "blocks_used": int(selected.sum()),
                   "mean_alt_fraction": float(afsum[selected_all].sum() / count[selected_all].sum()) if count[selected_all].sum() else None,
                   "mean_ad_depth": float(depth_sum[selected_all].sum() / count[selected_all].sum()) if count[selected_all].sum() else None}
            if chrom == "X" or selected.sum() < 5 or row["heterozygous_snvs"] < 1000:
                row["status"] = "descriptive_only_sex_baseline_unknown" if chrom == "X" else "insufficient_blocks"
                regions.append(row)
                continue
            baseline_hist = histograms[usable & (chroms != chrom)].sum(axis=0)
            baseline_n = baseline_hist.sum(axis=1)
            baseline_q = np.divide((baseline_hist * moment).sum(axis=1), baseline_n,
                                   out=np.zeros(101), where=baseline_n > 0)
            residual_q = (qsum - (depth_hist * baseline_q).sum(axis=1)) / count
            null, center = null_distribution(residual_q[usable], chroms[usable], chrom,
                                            selected.sum(), rng, iterations)
            observed = float(residual_q[selected].mean() - center)
            row.update(status="screened", excess_baf_variance=observed,
                       null_sd=float(np.std(null, ddof=1)),
                       p_empirical=float((1 + np.sum(null >= observed)) / (1 + iterations)),
                       q_effect_ci95_low=float(observed - np.quantile(null, .975)),
                       q_effect_ci95_high=float(observed - np.quantile(null, .025)),
                       _null=null)
            regions.append(row)
    tested = sum(row["status"] == "screened" for row in regions)
    for row in regions:
        if row["status"] != "screened":
            continue
        null = row.pop("_null")
        threshold = float(np.quantile(null, 1 - .05 / tested))
        # Inject expected delta² into an empirical null realization. This is an
        # aggregate signal model; it does not re-call genotypes or model dropout.
        injected_q90 = max(0.0, threshold - float(np.quantile(null, .1)))
        row.update(p_bonferroni=min(1.0, row["p_empirical"] * tested),
                   exploratory_outlier=row["p_empirical"] * tested <= .05,
                   q_familywise_threshold=threshold, q_injection_power90=injected_q90,
                   gain_fraction_power90=fraction_from_broadening(injected_q90, "gain"),
                   loss_fraction_power90=fraction_from_broadening(injected_q90, "loss"))
        for event in ("gain", "loss"):
            for fraction in (.05, .10, .20):
                injected = null + broadening_from_fraction(fraction, event)
                row[f"{event}_power_at_{int(fraction * 100)}pct"] = float(np.mean(injected > threshold))
        if row["arm"] == "whole":
            control.append(row["p_empirical"])
    total_auto = int(count[autosomal].sum())
    baseline_depth = float(depth_sum[autosomal].sum() / total_auto)
    for row in regions:
        row["variant_depth_ratio_to_autosomes"] = row["mean_ad_depth"] / baseline_depth if row["mean_ad_depth"] else None
    summary = {"autosomal_het_snvs": total_auto, "autosomal_mean_ad_depth": baseline_depth,
               "autosomal_mean_alt_fraction": float(afsum[autosomal].sum() / total_auto),
               "autosomal_usable_blocks": int(usable.sum()), "tests": tested,
               "bootstrap_iterations": iterations, "seed": seed,
               "exploratory_outlier_regions": [r["region"] for r in regions if r.get("exploratory_outlier")],
               "chromosome_pseudo_controls_p_below_0_05": sum(p < .05 for p in control),
               "pseudo_controls_total": len(control),
               "caution": "Same-genome pseudo-controls are not known negatives; error rates and detection limits are conditional, not clinically validated."}
    return regions, summary


def plot_regions(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chrom = [r for r in rows if r["arm"] == "whole" and r["status"] == "screened"]
    x = np.arange(len(chrom))
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, layout="constrained")
    values = np.array([r["excess_baf_variance"] for r in chrom])
    low = np.array([r["q_effect_ci95_low"] for r in chrom])
    high = np.array([r["q_effect_ci95_high"] for r in chrom])
    axes[0].errorbar(x, values * 1e4, yerr=np.vstack([values-low, high-values]) * 1e4,
                     fmt="o", color="#153e75", capsize=3)
    axes[0].axhline(0, color="#555555", lw=.8)
    axes[0].set_ylabel("Excess allelic variance × 10⁴\n(95% empirical-null interval)")
    axes[1].bar(x, [r["variant_depth_ratio_to_autosomes"] for r in chrom], color="#277d8e")
    axes[1].axhline(1, color="#555555", lw=.8)
    axes[1].set_ylabel("Depth at called het SNPs\n/ autosomal mean")
    axes[1].set_ylim(.8, 1.2)
    axes[1].set_xticks(x, [r["chromosome"] for r in chrom])
    axes[1].set_xlabel("Chromosome (GRCh38)")
    fig.suptitle("VCF-only mosaicism screen: no independent copy-number measurement", fontsize=13)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_report(rows, summary, output, block_mb):
    screened = [r for r in rows if r["status"] == "screened"]
    whole = [r for r in screened if r["arm"] == "whole"]
    gains = [r["gain_fraction_power90"] for r in whole if r["gain_fraction_power90"] is not None]
    losses = [r["loss_fraction_power90"] for r in whole if r["loss_fraction_power90"] is not None]
    outliers = summary["exploratory_outlier_regions"]
    conclusion = ("Exploratory outliers: " + ", ".join(outliers) + ". These require read-depth and/or cytogenetic confirmation.") if outliers else "No chromosome or arm exceeded the conservative familywise empirical-null threshold. This does not exclude MVA or mosaic aneuploidy."
    text = f"""# VCF-only mosaic aneuploidy screen

{conclusion}

## Observed data

- {summary['autosomal_het_snvs']:,} autosomal high-quality heterozygous SNVs; mean allele-depth {summary['autosomal_mean_ad_depth']:.2f}; mean alternate allele fraction {summary['autosomal_mean_alt_fraction']:.5f}.
- {summary['autosomal_usable_blocks']} eligible {block_mb}-Mb spatial blocks; {summary['tests']} chromosome/arm tests; {summary['bootstrap_iterations']:,} bootstrap replicates; seed {summary['seed']}.
- Same-genome pseudo-controls: {summary['chromosome_pseudo_controls_p_below_0_05']} / {summary['pseudo_controls_total']} whole chromosomes have unadjusted one-sided p<0.05. They are not known negative controls.
- Whole-chromosome 90%-power conditional gain fraction range: {min(gains)*100:.1f}–{max(gains)*100:.1f}%; loss: {min(losses)*100:.1f}–{max(losses)*100:.1f}%. Each region has its own limit in the CSV. These are model-based planning estimates, not patient-specific upper confidence bounds.

## Estimator and null

We use biallelic PASS SNVs with QUAL>=30, MQ>=40, GQ>=30, DP15–100, diploid heterozygous GT, and usable ref+alt depth15–100 and >=80% of DP. Both alleles must have at least one read. We do not filter on a narrow allele-balance range. Centromeres and 1-Mb flanks are excluded; blocks require >=500 SNVs; regions require >=5 blocks. X is descriptive only because sex/ploidy/PAR and tissue baselines are not established. Acrocentric short arms without sufficient markers are untestable.

Exact filter for this run: `{summary['filter']}`. The optional strict sensitivity also requires a dbSNP flag, MQ60, QD5, FS<=20, SOR<=2, QUAL50, GQ50 and DP20–80. Genotype and allele-depth safeguards still apply.

For alternate count a of n reads, q=((a−n/2)²−n/4)/(n(n−1)) is unbiased for (p−0.5)² under independent binomial reads. We subtract a leave-one-chromosome, exact-depth-stratified autosomal baseline. The region statistic weights spatial blocks equally. The null resamples within-chromosome centered {block_mb}-Mb block residuals and adds a randomly sampled comparator chromosome offset. This includes a conservative floor for between-chromosome technical variation. One-sided p-values test excess broadening; Bonferroni controls the nominal 0.05 family across eligible chromosomes and arms. Empirical support is only 21 comparator chromosomes. Confidence intervals are observed effect minus the 2.5/97.5% empirical-null quantiles; calibration relies on exchangeability.

## Conditional detection model

For a simple clonal single-homolog gain in fraction f, the true unphased bands are 0.5 ± f/[2(2+f)]; for loss, 0.5 ± f/[2(2−f)]. Expected broadening is delta². We inject this expected broadening into each empirical-null statistic. The reported 90%-power threshold is the familywise upper null quantile minus the lower 10% null quantile. This is an aggregate signal injection, not a genotype re-calling experiment. It assumes the same measurement noise, retained heterozygotes, and one altered homolog throughout the region. Loss/gain allele identity cannot be resolved without phasing or depth. Opposite-homolog events in different cells, balanced gains/losses, heterogeneous cell clones, and tissue-restricted events can cancel. The genome baseline could itself contain diffuse MVA signal; the screen detects regional contrasts and cannot measure total aneuploid cell burden.

## Depth and limitations

Depth is available only at selected variant loci. The reported depth ratios describe the VCF and cannot serve as GC-corrected read-depth CN calls: missing invariant sites, heterozygote ascertainment, mapping quality, GC/mappability and depth truncation all bias them. Ref/alt bias and local overdispersion can mimic BAF broadening. No matched public normal was processed for this VCF-only scan. The conditional detection estimates are not validated clinical sensitivity and must not be used to rule out low-level, nonclonal, or tissue-specific MVA. Read-depth alignment, a matched processed normal, and/or single-cell karyotyping are the next measurements.

## Reproduction and public sources

`UV_CACHE_DIR=/private/tmp/mva-uv-cache uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz`

Public GRCh38 arm boundaries: [UCSC cytoBand]({CYTO_URL}); cached under work/, SHA-256 in JSON. The source was downloaded and parsed; no patient information was sent.
Methods context, opened: Jakubek et al., *Mosaic chromosomal alterations in blood across ancestries using whole-genome sequencing*. PMID [37904051](https://pubmed.ncbi.nlm.nih.gov/37904051/), DOI [10.1038/s41588-023-01553-1](https://doi.org/10.1038/s41588-023-01553-1). This exploratory implementation is not MoChA and does not claim its validated performance.
"""
    output.write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--cytobands", type=Path, default=Path("work/mosaic_cytoBand.txt.gz"))
    parser.add_argument("--block-mb", type=int, default=5)
    parser.add_argument("--iterations", type=int, default=30000)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--qc", choices=("standard", "strict"), default="standard")
    parser.add_argument("--prefix", default="results/mosaic")
    args = parser.parse_args()
    if not 1 <= args.block_mb <= 20 or args.iterations < 5000:
        parser.error("block-mb must be 1..20 and iterations >=5000")
    args.cytobands.parent.mkdir(parents=True, exist_ok=True)
    if not args.cytobands.exists():
        import requests
        response = requests.get(CYTO_URL, timeout=45)
        response.raise_for_status()
        args.cytobands.write_bytes(response.content)
    filter_expression = STRICT_FILTER if args.qc == "strict" else FILTER
    cache = Path(f"work/mosaic_blocks_{args.block_mb}mb_{args.qc}.npz")
    metadata = cache.with_suffix(".json")
    source_identity = {"vcf_path": str(args.vcf.resolve()), "vcf_bytes": args.vcf.stat().st_size,
                       "vcf_mtime_ns": args.vcf.stat().st_mtime_ns, "filter": filter_expression,
                       "cytobands_sha256": hashlib.sha256(args.cytobands.read_bytes()).hexdigest()}
    if args.reuse_cache and cache.exists() and metadata.exists():
        saved = json.loads(metadata.read_text())
        if saved["source_identity"] != source_identity:
            raise ValueError("Cache input identity mismatch; rerun without --reuse-cache")
        data = np.load(cache)
        keys = [(str(c), str(a), int(b)) for c, a, b in data["keys"]]
        histograms = data["histograms"]
        qc, command = saved["qc"], saved["command"]
    else:
        keys, histograms, qc, command = stream_blocks(args.vcf, boundaries(args.cytobands), args.block_mb, filter_expression)
        np.savez_compressed(cache, keys=np.array(keys), histograms=histograms)
        metadata.write_text(json.dumps({"source_identity": source_identity, "qc": qc, "command": command}, indent=2))
    rows, summary = analyze(keys, histograms, args.iterations, args.seed)
    summary.update(qc=qc, block_mb=args.block_mb, filter=filter_expression,
                   cytobands_sha256=source_identity["cytobands_sha256"],
                   public_cytobands_url=CYTO_URL,
                   method_version="unphased-moment-spatial-and-chromosome-null-v1")
    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(f"{prefix}_regions.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    Path(f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(rows, summary, Path(f"{prefix}_report.md"), args.block_mb)
    plot_regions(rows, Path(f"{prefix}_chromosomes.png"))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
