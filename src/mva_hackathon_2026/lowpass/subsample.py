"""Validate and hash-subsample paired FASTQ without keeping names or genotypes.

Executing this script reads paired FASTQ and writes ignored intermediates.
Downloads are a separate, explicitly approved operation.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path


def records(path: Path):
    with gzip.open(path, "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not header.startswith(b"@") or not plus.startswith(b"+"):
                raise ValueError("Invalid FASTQ record; read content suppressed")
            if len(sequence.rstrip()) != len(quality.rstrip()) or not quality:
                raise ValueError("Incomplete FASTQ record; read content suppressed")
            yield header, sequence, plus, quality


def identity(header: bytes) -> bytes:
    return header.split()[0].removesuffix(b"/1").removesuffix(b"/2")


def pairs(left: Path, right: Path):
    from itertools import zip_longest

    for a, b in zip_longest(records(left), records(right)):
        if a is None or b is None or identity(a[0]) != identity(b[0]):
            raise ValueError("Mates are not synchronized; read content suppressed")
        yield a, b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, action="append", required=True)
    parser.add_argument("--r2", type=Path, action="append", required=True)
    parser.add_argument("--target-bases", type=int, default=3_400_000_000)
    parser.add_argument("--seed", default="mva-lowpass-v1")
    parser.add_argument("--prefix", type=Path, required=True)
    args = parser.parse_args()
    if len(args.r1) != len(args.r2):
        parser.error("Provide equal numbers of --r1 and --r2")
    work_root = (Path.cwd() / "work").resolve()
    data_root = (Path.cwd() / "data").resolve()
    if not args.prefix.resolve().is_relative_to(work_root):
        parser.error("Paired FASTQ intermediates must stay beneath the ignored work/ directory")
    if any(not (p.resolve().is_relative_to(data_root) or p.resolve().is_relative_to(work_root))
           for p in args.r1 + args.r2):
        parser.error("FASTQ inputs must stay under data/ (originals) or work/ (derived subsamples)")
    if args.target_bases <= 0:
        parser.error("target-bases must be positive")
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    out1, out2 = (Path(str(args.prefix) + suffix) for suffix in (".R1.fastq.gz", ".R2.fastq.gz"))
    if out1.exists() or out2.exists():
        raise FileExistsError("Refusing to replace an existing subsample")
    total_pairs = total_bases = 0
    for left, right in zip(args.r1, args.r2, strict=True):
        for a, b in pairs(left, right):
            total_pairs += 1
            total_bases += len(a[1].rstrip()) + len(b[1].rstrip())
    fraction = min(1.0, args.target_bases / total_bases)
    kept_pairs = kept_bases = 0
    with gzip.open(out1, "wb", compresslevel=1) as f1, gzip.open(out2, "wb", compresslevel=1) as f2:
        for lane, (left, right) in enumerate(zip(args.r1, args.r2, strict=True)):
            for i, (a, b) in enumerate(pairs(left, right)):
                # Index and lane avoid selecting correlated subsets across files
                # with repeated read IDs. Both mates share one random decision.
                digest = hashlib.blake2b(f"{args.seed}:{lane}:{i}".encode(), digest_size=8).digest()
                if int.from_bytes(digest, "big") / 2**64 < fraction:
                    f1.writelines(a)
                    f2.writelines(b)
                    kept_pairs += 1
                    kept_bases += len(a[1].rstrip()) + len(b[1].rstrip())
                if i % 100_000 == 0 and shutil.disk_usage(args.prefix.parent).free < 32_000_000_000:
                    raise RuntimeError("STOP: free disk is approaching the 30-GB floor")
    summary = dict(total_pairs=total_pairs, total_bases=total_bases, selection_probability=fraction,
                   seed=args.seed, kept_pairs=kept_pairs, kept_bases=kept_bases,
                   nominal_coverage_grch38=kept_bases / 3_088_269_832,
                   note="Nominal input coverage; aligned usable coverage will be lower")
    Path(str(args.prefix) + ".sampling.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
