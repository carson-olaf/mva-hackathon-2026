#!/usr/bin/env python3
"""Reproduce the synthetic low-pass alignment checks; no patient input is used.

Run native tools through the host build lock, after other alignment work ends:

    ~/.agents/bin/build-lock uv run python scripts/validate_lowpass_alignment.py

Requires minimap2 and samtools on PATH plus an installed mva_hackathon_2026.
Use --work-dir to select a fresh output directory. Existing directories are never
removed or reused. The script prints aggregate JSON and saves summary.json.
Intermediate reads, SAM/BAM files and logs contain synthetic sequence only.

The original validated fixture used seed 20261024, minimap2 2.31-r1302 and
samtools 1.24. Its one-batch and 25-pair runs retained 120 of 130 pairs after
global removal of 10 duplicate pairs. Their proper-pair record counts differed
(178 versus 240) among ambiguous repeats; all repeats had MAPQ 0. Unique Q30
alignments matched exactly. This is NOT full alignment identity or a validation
of patient mapping accuracy, production memory, runtime or mosaic sensitivity.
The current run measures these quantities again rather than copying old results.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess

from mva_hackathon_2026.lowpass.align import align

SEED = 20261024


def fixture(directory: Path) -> dict[str, str]:
    """Write the original 120-pair fixture, then append 10 duplicate templates."""
    rng = random.Random(SEED)

    def dna(length):
        return "".join(rng.choice("ACGT") for _ in range(length))

    a, b, repeat = dna(50_000), dna(50_000), dna(5_000)
    a = a[:10_000] + repeat + a[15_000:]
    b = b[:10_000] + repeat + b[15_000:]
    (directory / "reference.fa").write_text(f">synthA\n{a}\n>synthB\n{b}\n")
    left, right, kinds = [], [], {}
    for i in range(120):
        if i < 80:
            sequence = a if i < 40 else b
            start, kind = 20_000 + (i % 40) * 450, "unique"
        else:
            sequence, start, kind = a, 10_500 + (i - 80) * 80, "repeat"
        name = f"synthetic_{kind}_{i:03}"
        kinds[name] = kind
        r1 = sequence[start:start + 150]
        r2 = sequence[start + 350:start + 500].translate(str.maketrans("ACGT", "TGCA"))[::-1]
        left.append(f"@{name}/1\n{r1}\n+\n" + "I" * 150 + "\n")
        right.append(f"@{name}/2\n{r2}\n+\n" + "I" * 150 + "\n")
    for mate, reads in enumerate((left, right), 1):
        (directory / f"R{mate}.fq").write_text("".join(reads))
        duplicate_reads = [record.replace("@synthetic_", "@duplicate_", 1)
                           for record in reads[:10]]
        (directory / f"R{mate}.with_duplicates.fq").write_text("".join(reads + duplicate_reads))
    for i in range(10):
        kinds[f"duplicate_unique_{i:03}"] = "unique"
    return kinds


def parse_sam(text: str) -> tuple[list[dict], int]:
    rows, contigs = [], 0
    for line in text.splitlines():
        if line.startswith("@SQ"):
            contigs += 1
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & 0x900:
            continue
        rows.append({"name": fields[0], "flag": flag, "chrom": fields[2],
                     "pos": int(fields[3]), "mapq": int(fields[4]), "cigar": fields[5]})
    return rows, contigs


def require(condition: bool, message: str):
    if not condition:
        raise RuntimeError(message)


def alignment_tuple(row):
    # Ignore representative name: identical duplicate templates may tie.
    return row["flag"], row["chrom"], row["pos"], row["mapq"], row["cigar"]


def split_check(directory: Path, kinds: dict[str, str]) -> dict:
    parsed, metrics = {}, {}
    for name, options in (("single_index", ["-I", "1M"]),
                          ("split_index", ["-I", "40k", "--split-prefix", str(directory / "split-merge")])):
        command = ["minimap2", "-ax", "sr", "-K", "25M", "-t", "2", "--secondary=no",
                   *options, str(directory / "reference.fa"), str(directory / "R1.fq"),
                   str(directory / "R2.fq")]
        with (directory / f"{name}.sam").open("w") as out, (directory / f"{name}.log").open("w") as log:
            subprocess.run(command, stdout=out, stderr=log, check=True)
        rows, contigs = parse_sam((directory / f"{name}.sam").read_text())
        unique = [row for row in rows if kinds[row["name"]] == "unique"]
        repeats = [row for row in rows if kinds[row["name"]] == "repeat"]
        require(contigs == 2 and len(rows) == 240, "Split-index header or primary count changed")
        require(len(unique) == 160 and all(row["mapq"] >= 50 for row in unique), "Unique mapping quality changed")
        require(len(repeats) == 80 and all(row["mapq"] == 0 for row in repeats), "Ambiguous repeat MAPQ changed")
        parsed[name] = {(row["name"], bool(row["flag"] & 64)): row for row in rows}
        metrics[name] = {"header_contigs": contigs, "primary_records": len(rows),
                         "q30_records": sum(row["mapq"] >= 30 for row in rows),
                         "proper_pair_records": sum(bool(row["flag"] & 2) for row in rows),
                         "repeat_proper_pair_records": sum(bool(row["flag"] & 2) for row in repeats)}
    x, y = parsed["single_index"], parsed["split_index"]
    require(x.keys() == y.keys(), "Single and split read identities changed")
    unique_match = all(alignment_tuple(row) == alignment_tuple(y[key]) for key, row in x.items()
                       if kinds[row["name"]] == "unique")
    require(unique_match, "Unique single/split alignment, flag, CIGAR or MAPQ changed")
    return {"input_pairs": 120, "unique_pairs": 80, "ambiguous_repeat_pairs": 40,
            "runs": metrics, "unique_q30_alignments_flags_cigar_mapq_exact": unique_match,
            "all_ambiguous_repeats_mapq_zero": True,
            "ambiguous_repeat_proper_pair_flag_difference_observed": any(
                bool(row["flag"] & 2) != bool(y[key]["flag"] & 2)
                for key, row in x.items() if kinds[row["name"]] == "repeat"),
            "full_alignment_identity_claimed": False}


def batched_check(directory: Path, kinds: dict[str, str]) -> dict:
    metrics, unique_positions = {}, {}
    for name, pairs, index_size in (("one_batch", 1000, "1M"), ("batches_of_25", 25, "40k")):
        output = directory / f"{name}.bam"
        align(directory / "reference.fa", directory / "R1.with_duplicates.fq",
              directory / "R2.with_duplicates.fq", directory / f"{name}-work", output,
              pairs=pairs, index_size=index_size, threads=2)
        rows, _ = parse_sam(subprocess.check_output(["samtools", "view", str(output)], text=True))
        unique = [row for row in rows if kinds[row["name"]] == "unique"]
        repeats = [row for row in rows if kinds[row["name"]] == "repeat"]
        require(len(rows) == 240, "Global duplicate removal did not retain 120 pairs")
        require(len(unique) == 160 and all(row["mapq"] >= 50 for row in unique), "Unique mapping quality changed")
        require(len(repeats) == 80 and all(row["mapq"] == 0 for row in repeats), "Ambiguous repeat MAPQ changed")
        require(set(Counter(row["name"] for row in rows).values()) == {2}, "Retained primary records are not paired")
        unique_positions[name] = sorted(alignment_tuple(row) for row in rows if row["mapq"] >= 30)
        metrics[name] = {"primary_records": len(rows), "retained_pairs": len(rows) // 2,
                         "q30_unique_records": len(unique), "ambiguous_repeat_records": len(repeats),
                         "proper_pair_records": sum(bool(row["flag"] & 2) for row in rows),
                         "repeat_proper_pair_records": sum(bool(row["flag"] & 2) for row in repeats),
                         "pairs_per_batch": pairs, "index_batch": index_size}
    exact = unique_positions["one_batch"] == unique_positions["batches_of_25"]
    require(exact, "Unique Q30 alignments changed across batches")
    return {"input_pairs": 130, "cross_batch_duplicate_pairs": 10,
            "duplicate_originals_in_first_batch_duplicates_after_original_120": True,
            "global_duplicate_pairs_removed": 10, "global_duplicate_removal_verified": True,
            "runs": metrics, "unique_q30_alignments_flags_cigar_mapq_exact_ignoring_duplicate_name": exact,
            "all_ambiguous_repeats_mapq_zero": True,
            "ambiguous_repeat_proper_pair_flag_counts_differ":
                metrics["one_batch"]["repeat_proper_pair_records"] != metrics["batches_of_25"]["repeat_proper_pair_records"],
            "full_alignment_identity_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work-dir", type=Path, default=Path("work/validation/public_lowpass_alignment"),
                        help="New directory for synthetic intermediates and aggregate summary.json")
    args = parser.parse_args()
    for program in ("minimap2", "samtools"):
        if not shutil.which(program):
            parser.error(f"Required program is not on PATH: {program}")
    directory = args.work_dir
    if directory.exists():
        parser.error("--work-dir already exists; choose a fresh directory (nothing was removed)")
    directory.mkdir(parents=True)
    kinds = fixture(directory)
    summary = {"seed": SEED, "patient_data_read": False,
               "minimap2_version": subprocess.check_output(["minimap2", "--version"], text=True).strip(),
               "samtools_version": subprocess.check_output(["samtools", "--version"], text=True).splitlines()[0],
               "reference_sha256": hashlib.sha256((directory / "reference.fa").read_bytes()).hexdigest(),
               "single_versus_split_index": split_check(directory, kinds),
               "one_batch_versus_25_pair_batches": batched_check(directory, kinds),
               "caveat": "Ambiguous repeat placements and proper-pair flags can differ while all have MAPQ 0; Q30 depth excludes them. Unique Q30 records are compared exactly, ignoring the identity of interchangeable duplicate representatives. No full alignment identity is claimed.",
               "scope": "Synthetic mapping and global duplicate-handling contract only; not patient mapping accuracy, memory/runtime, or aneuploidy sensitivity validation."}
    encoded = json.dumps(summary, indent=2) + "\n"
    (directory / "summary.json").write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
