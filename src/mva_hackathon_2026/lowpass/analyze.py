"""GC-corrected chromosome/arm depth contrasts and conditional detection limits.

This is a prespecified exploratory estimator, not a validated clinical CN caller.
The public normal is a technical comparator, not a certified euploid blood control.
No per-bin input or output is published; only chromosome/arm aggregates can leave
the ignored work directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
CHROMOSOMES = (*AUTOSOMES, "chrX", "chrY")
ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def gc_residual(gc: np.ndarray, depth: np.ndarray):
    """Fit robust cubic log-depth on GC; never use the test region's label."""
    use = np.isfinite(depth) & (depth > 0)
    response = np.log2(np.maximum(depth, 1e-12))
    if use.sum() < 100:
        raise ValueError("Insufficient positive-depth bins for GC fit")
    for _ in range(4):
        fit = np.polynomial.polynomial.polyfit(gc[use], response[use], deg=3)
        residual = response - np.polynomial.polynomial.polyval(gc, fit)
        center = np.median(residual[use])
        mad = max(1.4826 * np.median(abs(residual[use] - center)), 0.01)
        use &= abs(residual - center) < 4 * mad
    corrected = response - np.polynomial.polynomial.polyval(gc, fit)
    return corrected, {"cubic_coefficients": fit.tolist(), "fit_bins": int(use.sum()),
                       "pearson_gc_logdepth_before": float(np.corrcoef(gc, response)[0, 1]),
                       "pearson_gc_logdepth_after": float(np.corrcoef(gc, corrected)[0, 1])}


def prepare(metadata: pd.DataFrame, coverage: pd.DataFrame):
    coordinates = ["chrom", "start", "end"]
    if not set(AUTOSOMES).issubset(set(metadata.chrom)):
        raise ValueError("Reference metadata must include all22 autosomes")
    if metadata.duplicated(coordinates).any() or (metadata.end <= metadata.start).any():
        raise ValueError("Duplicate or invalid reference-bin coordinates")
    for _, group in metadata.groupby("chrom", sort=False):
        if not group.start.is_monotonic_increasing or (group.start.to_numpy()[1:] < group.end.to_numpy()[:-1]).any():
            raise ValueError("Reference bins must be ordered and nonoverlapping")
    if not metadata[coordinates].equals(coverage[coordinates]):
        raise ValueError("BED coverage and reference metadata bins do not match exactly")
    frame = metadata.copy()
    frame["sample_depth"] = coverage["sample"] / (frame.end - frame.start)
    frame["normal_depth"] = coverage["normal"] / (frame.end - frame.start)
    descriptive = frame.groupby("chrom")[["sample_depth", "normal_depth"]].median()
    auto = frame.chrom.isin(AUTOSOMES)
    lo, hi = frame.loc[auto, "normal_depth"].quantile([.01, .99])
    keep = auto & frame.normal_depth.between(max(lo, .1), hi) & (frame.sample_depth > 0)
    frame = frame.loc[keep].copy()
    if len(frame) < 1000:
        raise ValueError("Insufficient usable autosomal bins")
    sample, sample_qc = gc_residual(frame.gc.to_numpy(), frame.sample_depth.to_numpy())
    normal, normal_qc = gc_residual(frame.gc.to_numpy(), frame.normal_depth.to_numpy())
    frame["sample_gc_relative_log2"] = sample - np.median(sample)
    frame["normal_gc_relative_log2"] = normal - np.median(normal)
    frame["log2_ratio"] = sample - normal
    frame["log2_ratio"] -= frame.log2_ratio.median()
    qc = {"input_bins": len(metadata), "retained_autosome_bins": len(frame),
          "normal_depth_mask_min": float(max(lo, .1)), "normal_depth_mask_max": float(hi),
          "sample_gc_fit": sample_qc, "normal_gc_fit": normal_qc,
          "sample_usable_mean_depth": float(frame.sample_depth.mean()),
          "normal_usable_mean_depth": float(frame.normal_depth.mean()),
          "raw_median_autosomal_log2_ratio": float(np.median(np.log2(frame.sample_depth / frame.normal_depth))),
          "descriptive_gc_normalized_depth_by_chromosome": {
              chrom: {"sample": float(2**group.sample_gc_relative_log2.median()),
                      "normal": float(2**group.normal_gc_relative_log2.median()), "bins": len(group)}
              for chrom, group in frame.groupby("chrom", sort=False)},
          "individual_profile_caveat": "Descriptive within-sample normalized depth only; no absolute ploidy or tissue/library correction"}
    return frame, descriptive, qc


def resampled_medians(values, count, rng, iterations):
    out = np.empty(iterations)
    for start in range(0, iterations, 500):
        end = min(iterations, start + 500)
        out[start:end] = np.median(rng.choice(values, (end - start, count)), axis=1)
    return out


def empirical_null(blocks, chromosome, count, rng, iterations):
    """Other-chromosome offsets plus centered spatial blocks make a batch floor.

    The resampling assumptions are conditional. At most 21 comparator chromosome
    offsets exist; 30,000 draws do not create 30,000 independent control genomes.
    """
    other = blocks.loc[blocks.chrom != chromosome]
    center = float(other.log2_ratio.median())
    offsets, residuals = [], []
    for _, group in other.groupby("chrom", sort=False):
        values = group.log2_ratio.to_numpy()
        offset = float(np.median(values))
        offsets.append(offset - center)
        residuals.extend(values - offset)
    null = rng.choice(offsets, iterations) + resampled_medians(np.asarray(residuals), count, rng, iterations)
    return null, center


def detection_fraction(shift, direction):
    # A single one-copy gain/loss in fraction f gives R = 1 +/- f/2.
    fraction = 2 * (2**shift - 1) if direction == "gain" else 2 * (1 - 2**(-shift))
    return float(fraction) if 0 <= fraction <= 1 else None


def analyze(frame: pd.DataFrame, *, block_mb=1, iterations=30000, seed=20261024):
    frame = frame.copy()
    frame["block"] = frame.start // (block_mb * 1_000_000)
    grouped = frame.groupby(["chrom", "arm", "block"], sort=False).log2_ratio
    blocks = grouped.agg(log2_ratio="median", retained_bins="size").reset_index()
    bin_width = float(np.median(frame.end - frame.start))
    minimum_bins = max(1, math.ceil(.5 * block_mb * 1_000_000 / bin_width))
    available_blocks = len(blocks)
    # A nearly empty masked block must not get the same weight as a full block.
    blocks = blocks.loc[blocks.retained_bins >= minimum_bins].copy()
    if len(blocks) < 100:
        raise ValueError("Insufficient spatial blocks for calibration")
    regions = []
    for chrom in CHROMOSOMES:
        for arm in ("whole", "p", "q"):
            group = blocks.loc[(blocks.chrom == chrom) & ((blocks.arm == arm) if arm != "whole" else True)]
            group_bins = frame.loc[(frame.chrom == chrom) & ((frame.arm == arm) if arm != "whole" else True)]
            row = {"region": chrom + (arm if arm != "whole" else ""), "chromosome": chrom,
                   "arm": arm, "bins": len(group_bins), "blocks": len(group), "_values": group.log2_ratio.to_numpy()}
            if chrom not in AUTOSOMES:
                row["status"] = "descriptive_only_sex_baseline_unknown"
            elif chrom in ACROCENTRIC and arm == "p":
                row["status"] = "acrocentric_p_arm_excluded"
            elif len(group) < 5:
                row["status"] = "insufficient_blocks"
            else:
                row["status"] = "screened"
            regions.append(row)
    tests = sum(r["status"] == "screened" for r in regions)
    rng = np.random.default_rng(seed)
    for row in regions:
        values = row.pop("_values")
        if row["status"] != "screened":
            continue
        null, center = empirical_null(blocks, row["chromosome"], len(values), rng, iterations)
        observed = float(np.median(values) - center)
        threshold = float(np.quantile(abs(null), 1 - .05 / tests))
        p = float((1 + np.count_nonzero(abs(null) >= abs(observed))) / (iterations + 1))
        boot = resampled_medians(values, len(values), rng, iterations) - center
        within = 2 * 2**np.quantile(boot, [.025, .975])
        technical = 2 * 2**(observed - np.quantile(null, [.975, .025]))
        gain_shift = max(0., threshold - float(np.quantile(null, .1)))
        loss_shift = max(0., threshold + float(np.quantile(null, .9)))
        row.update(log2_ratio=observed, relative_cn=float(2 * 2**observed),
                   relative_cn_block_bootstrap_ci95_low=float(within[0]),
                   relative_cn_block_bootstrap_ci95_high=float(within[1]),
                   relative_cn_technical_null_interval95_low=float(technical[0]),
                   relative_cn_technical_null_interval95_high=float(technical[1]),
                   technical_null_p=p, bonferroni_p=min(1., p * tests),
                   exploratory_outlier=p * tests <= .05,
                   familywise_abs_log2_threshold=threshold,
                   gain_fraction_power90=detection_fraction(gain_shift, "gain"),
                   loss_fraction_power90=detection_fraction(loss_shift, "loss"),
                   null_sd=float(np.std(null, ddof=1)))
        for fraction in (.025, .05, .1, .2):
            for direction in ("gain", "loss"):
                shift = np.log2(1 + fraction / 2) if direction == "gain" else np.log2(1 - fraction / 2)
                row[f"{direction}_power_at_{fraction:g}_fraction"] = float(np.mean(abs(null + shift) > threshold))
    summary = {"block_mb": block_mb, "iterations": iterations, "seed": seed,
               "tests": tests, "autosomal_blocks": len(blocks),
               "minimum_retained_bins_per_spatial_block": minimum_bins,
               "sparse_spatial_blocks_excluded": available_blocks - len(blocks),
               "exploratory_outlier_regions": [r["region"] for r in regions if r.get("exploratory_outlier")],
               "conditional_null": "Other-chromosome median offsets plus centered spatial block medians",
               "absolute_ploidy_or_total_mva_burden_estimated": False}
    return regions, summary


def write_plot(rows, output: Path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chosen = [r for r in rows if r["arm"] == "whole" and r["status"] == "screened"]
    x = np.arange(len(chosen))
    y = np.array([r["relative_cn"] for r in chosen])
    lower = np.array([r["relative_cn_technical_null_interval95_low"] for r in chosen])
    upper = np.array([r["relative_cn_technical_null_interval95_high"] for r in chosen])
    fig, ax = plt.subplots(figsize=(11, 4.8), layout="constrained")
    ax.errorbar(x, y, yerr=[y - lower, upper - y], fmt="o", capsize=3, color="#143d6b")
    ax.axhline(2, lw=1, color="#666666")
    ax.set_xticks(x, [r["chromosome"].removeprefix("chr") for r in chosen])
    ax.set_xlabel("Autosome (GRCh38)")
    ax.set_ylabel("Relative CN (diploid baseline assumed)\n95% empirical-null interval")
    ax.set_title(title)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--bedcov", type=Path, required=True)
    p.add_argument("--work-prefix", type=Path, required=True)
    p.add_argument("--public-prefix", type=Path, required=True)
    p.add_argument("--comparison", default="sample_vs_public_normal")
    p.add_argument("--block-mb", type=int, default=1, choices=(1, 5))
    p.add_argument("--iterations", type=int, default=30000)
    p.add_argument("--seed", type=int, default=20261024)
    args = p.parse_args()
    private_root = (Path.cwd() / "work").resolve()
    if not args.work_prefix.resolve().is_relative_to(private_root):
        p.error("Per-bin intermediates must stay beneath the ignored work/ directory")
    if args.iterations < 5000:
        p.error("At least5000 resampling iterations are required")
    args.work_prefix.parent.mkdir(parents=True, exist_ok=True)
    args.public_prefix.parent.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.metadata)
    manifest_path = args.metadata.with_suffix(".manifest.json")
    bin_manifest = json.loads(manifest_path.read_text())
    bed_path = args.metadata.with_suffix(".bed")
    if (bin_manifest.get("complete") is not True or bin_manifest["rows"] != len(metadata)
            or bin_manifest["csv_sha256"] != hashlib.sha256(args.metadata.read_bytes()).hexdigest()
            or bin_manifest["bed_sha256"] != hashlib.sha256(bed_path.read_bytes()).hexdigest()
            or any(int(bin_manifest["bins_per_chromosome"].get(chrom, 0)) < 10 for chrom in AUTOSOMES)):
        raise ValueError("Reference-bin completeness/provenance manifest mismatch")
    actual_counts = metadata.groupby("chrom").size().to_dict()
    if actual_counts != bin_manifest["bins_per_chromosome"]:
        raise ValueError("Reference-bin chromosome counts differ from complete manifest")
    coverage = pd.read_csv(args.bedcov, sep="\t", names=["chrom", "start", "end", "sample", "normal"])
    frame, descriptive, qc = prepare(metadata, coverage)
    rows, summary = analyze(frame, block_mb=args.block_mb, iterations=args.iterations, seed=args.seed)
    for data, destination, include_index in [
        (frame, Path(str(args.work_prefix) + "_bins.csv"), False),
        (descriptive, Path(str(args.work_prefix) + "_raw_depth.csv"), True),
        (pd.DataFrame(rows), Path(str(args.public_prefix) + "_regions.csv"), False),
    ]:
        temporary = Path(str(destination) + ".partial")
        data.to_csv(temporary, index=include_index)
        temporary.replace(destination)
    summary.update(comparison=args.comparison, qc=qc,
                   metadata_sha256=hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
                   bin_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                   source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   private_bedcov_sha256=hashlib.sha256(args.bedcov.read_bytes()).hexdigest(),
                   method_version="gc-cubic-spatial-and-chromosome-null-v1",
                   limitations=["Single public normal; different tissue/library and no clinical validation",
                                "Technical null assumes exchangeable comparator chromosomes and spatial blocks",
                                "At most21 chromosome offsets; resampling draws do not create independent controls",
                                "Conditional detection limit, not a patient mosaicism upper confidence bound",
                                "Relative CN assumes most autosomes are diploid; absolute ploidy is unresolved",
                                "No total MVA cell burden; opposing clones and tissue-restricted changes can be missed",
                                "GC is modeled; no independently measured genome-wide mappability track",
                                "Sex chromosomes descriptive only; acrocentric p arms excluded"])
    plot = Path(str(args.public_prefix) + "_chromosomes.png")
    plot_temporary = plot.with_name(plot.stem + ".partial.png")
    write_plot(rows, plot_temporary, args.comparison)
    plot_temporary.replace(plot)
    public_table = Path(str(args.public_prefix) + "_regions.csv")
    summary["public_output_sha256"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (public_table, plot)
    }
    destination = Path(str(args.public_prefix) + "_summary.json")
    temporary = Path(str(destination) + ".partial")
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary.replace(destination)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
