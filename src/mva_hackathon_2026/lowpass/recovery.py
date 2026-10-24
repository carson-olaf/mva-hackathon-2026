"""Measure raw-depth signal recovery through GC correction; no clinical power claim."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .analyze import AUTOSOMES, prepare
from .provenance import atomic_json


def block_statistics(frame: pd.DataFrame) -> dict:
    """Match analyze()'s 1-Mb whole-chromosome statistic without resampling.

    The arm key and half-expected-bin rule are deliberate: centromere-adjacent
    blocks must receive the same treatment as the main chromosome analysis.
    """
    frame = frame.copy()
    frame["block"] = frame.start // 1_000_000
    blocks = frame.groupby(["chrom", "arm", "block"], sort=False).log2_ratio.agg(
        log2_ratio="median", retained_bins="size").reset_index()
    minimum = max(1, math.ceil(.5 * 1_000_000 / float(np.median(frame.end - frame.start))))
    blocks = blocks.loc[blocks.retained_bins >= minimum]
    if len(blocks) < 100:
        raise ValueError("Insufficient spatial blocks for recovery comparison")
    output = {}
    for chrom in AUTOSOMES:
        values = blocks.loc[blocks.chrom == chrom, "log2_ratio"]
        other = blocks.loc[blocks.chrom != chrom, "log2_ratio"]
        output[chrom] = {"statistic": float(values.median() - other.median()) if len(values) >= 5 else None,
                         "retained_bins": int((frame.chrom == chrom).sum()), "blocks": len(values)}
    return output


def recover(metadata: pd.DataFrame, coverage: pd.DataFrame, *, fractions=(.10, .20, .30),
            chromosomes=AUTOSOMES):
    """Inject into covered bases before repeating all prepare() operations."""
    if any(not 0 < fraction <= 1 for fraction in fractions):
        raise ValueError("Injected single-copy mosaic fractions must be in (0, 1]")
    if not set(chromosomes).issubset(AUTOSOMES):
        raise ValueError("Recovery is limited to autosomes")
    baseline_frame, _, _ = prepare(metadata, coverage)
    baseline = block_statistics(baseline_frame)
    rows = []
    for chrom in chromosomes:
        for fraction in fractions:
            for direction, sign in (("gain", 1), ("loss", -1)):
                factor = 1 + sign * fraction / 2
                injected = coverage.copy()
                # Fractional covered bases are intentional, deterministic depth
                # perturbations, not simulated or resampled read alignments.
                injected["sample"] = injected["sample"].astype(float)
                injected.loc[injected.chrom == chrom, "sample"] *= factor
                spiked_frame, _, _ = prepare(metadata, injected)
                spiked = block_statistics(spiked_frame)[chrom]
                base = baseline[chrom]
                expected = float(np.log2(factor))
                detectable = base["statistic"] is not None and spiked["statistic"] is not None
                delta = spiked["statistic"] - base["statistic"] if detectable else None
                rows.append({"chromosome": chrom, "direction": direction, "injected_mosaic_fraction": fraction,
                             "expected_log2_shift": expected, "baseline_statistic": base["statistic"],
                             "spiked_statistic": spiked["statistic"], "observed_statistic_delta": delta,
                             "recovery_ratio": delta / expected if detectable else None,
                             "baseline_retained_target_bins": base["retained_bins"],
                             "spiked_retained_target_bins": spiked["retained_bins"],
                             "retained_target_bins_change": spiked["retained_bins"] - base["retained_bins"],
                             "baseline_retained_genome_bins": len(baseline_frame),
                             "spiked_retained_genome_bins": len(spiked_frame),
                             "retained_genome_bins_change": len(spiked_frame) - len(baseline_frame),
                             "baseline_target_blocks": base["blocks"], "spiked_target_blocks": spiked["blocks"],
                             "status": "evaluated" if detectable else "insufficient_blocks"})
    return rows


def _validated_inputs(metadata_path, coverage_path):
    metadata = pd.read_csv(metadata_path)
    manifest_path = metadata_path.with_suffix(".manifest.json")
    manifest = json.loads(manifest_path.read_text())
    bed = metadata_path.with_suffix(".bed")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    if (manifest.get("complete") is not True or manifest["rows"] != len(metadata)
            or manifest["csv_sha256"] != digest(metadata_path) or manifest["bed_sha256"] != digest(bed)
            or any(int(manifest["bins_per_chromosome"].get(chrom, 0)) < 10 for chrom in AUTOSOMES)
            or metadata.groupby("chrom").size().to_dict() != manifest["bins_per_chromosome"]
            or int(manifest["bin_size"]) != 100_000
            or not ((metadata.end - metadata.start) == 100_000).all()):
        raise ValueError("Complete 100-kb reference-bin manifest does not match metadata")
    coverage = pd.read_csv(coverage_path, sep="\t", names=["chrom", "start", "end", "sample", "normal"])
    # prepare() enforces complete chromosome presence, ordered nonoverlapping
    # metadata and exact coordinate equality with all coverage rows.
    return metadata, coverage, manifest_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--bedcov", type=Path, required=True,
                        help="Mapped-base-matched normal self-split coverage, not patient coverage")
    parser.add_argument("--public-prefix", type=Path, required=True)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    private = root / "work"
    public = root / "results"
    input_paths = (args.metadata, args.bedcov, args.metadata.with_suffix(".bed"),
                   args.metadata.with_suffix(".manifest.json"))
    if any(not path.resolve().is_relative_to(private) for path in input_paths):
        parser.error("Private metadata and coverage inputs must stay beneath work/")
    if not args.public_prefix.resolve().is_relative_to(public):
        parser.error("Aggregate outputs must stay beneath results/")
    # Summary is the completion marker. An interrupted rerun must not retain an
    # old summary next to a newly replaced or partial current CSV.
    Path(str(args.public_prefix) + "_summary.json").unlink(missing_ok=True)
    metadata, coverage, manifest = _validated_inputs(args.metadata, args.bedcov)
    rows = recover(metadata, coverage)
    table = pd.DataFrame(rows)
    csv_path = Path(str(args.public_prefix) + "_recovery.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = csv_path.with_name(csv_path.name + ".partial")
    table.to_csv(temporary, index=False)
    temporary.replace(csv_path)
    ratios = table.loc[table.status == "evaluated", "recovery_ratio"]
    summary = {"generated_utc": datetime.now(timezone.utc).isoformat(),
               "method": "pre-GC deterministic covered-base injection recovery",
               "baseline_null": "Unspiked mapped-base-matched normal self-split; supplied by caller",
               "block_mb": 1, "bin_size": 100_000, "injected_fractions": [.1, .2, .3],
               "rows": len(rows), "evaluated_rows": len(ratios),
               "recovery_ratio_min": float(ratios.min()) if len(ratios) else None,
               "recovery_ratio_median": float(ratios.median()) if len(ratios) else None,
               "recovery_ratio_max": float(ratios.max()) if len(ratios) else None,
               "nominal_lods_adjusted": False,
               "source_hashes": {"metadata": hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
                                 "private_bedcov": hashlib.sha256(args.bedcov.read_bytes()).hexdigest(),
                                 "bin_manifest": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                                 "recovery_code": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                 "prepare_code": hashlib.sha256(Path(__file__).with_name("analyze.py").read_bytes()).hexdigest(),
                                 "aggregate_csv": hashlib.sha256(csv_path.read_bytes()).hexdigest()},
               "caveats": ["Deterministic depth perturbations do not validate alignment, read sampling or mapping bias",
                           "Tests GC correction and refiltering against this control only; not clinical sensitivity or power",
                           "Single copy gain/loss in fraction f assumes ratio1+f/2 or1-f/2 and diploid baseline",
                           "Baseline self-split is not an independent biological or tissue/library-matched control",
                           "Recovery ratios can be above1 or negative; no nominal LOD is automatically corrected"]}
    atomic_json(Path(str(args.public_prefix) + "_summary.json"), summary)
    print(json.dumps({"rows": len(rows), "evaluated_rows": len(ratios), "nominal_lods_adjusted": False}))


if __name__ == "__main__":
    main()
