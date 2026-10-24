"""Make 100-kb GRCh38 primary-chromosome bins and GC metadata, locally."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--cytoband", type=Path, required=True)
    p.add_argument("--prefix", type=Path, required=True)
    p.add_argument("--bin-size", type=int, default=100_000)
    a = p.parse_args()
    if a.bin_size <= 0:
        p.error("bin-size must be positive")
    primary = {f"chr{x}" for x in range(1, 23)} | {"chrX", "chrY"}
    centromeres, lengths = {}, {}
    with gzip.open(a.cytoband, "rt") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[0] in primary:
                lengths[row[0]] = max(lengths.get(row[0], 0), int(row[2]))
            if row[0] in primary and row[4] == "acen":
                lo, hi = centromeres.get(row[0], (int(row[1]), int(row[2])))
                centromeres[row[0]] = min(lo, int(row[1])), max(hi, int(row[2]))
    a.prefix.parent.mkdir(parents=True, exist_ok=True)
    bed_path, csv_path = Path(str(a.prefix) + ".bed"), Path(str(a.prefix) + ".csv")
    bed_temp, csv_temp = Path(str(bed_path) + ".partial"), Path(str(csv_path) + ".partial")
    seen, counts = set(), Counter()
    with bed_temp.open("w") as bed, csv_temp.open("w") as metadata:
        writer = csv.writer(metadata)
        writer.writerow(["chrom", "start", "end", "arm", "gc", "n_fraction"])

        def emit(chrom, sequence):
            if chrom not in primary:
                return
            if chrom in seen or len(sequence) != lengths.get(chrom):
                raise ValueError("Duplicate/incomplete primary reference chromosome")
            seen.add(chrom)
            cen_lo, cen_hi = centromeres[chrom]
            for start in range(0, len(sequence), a.bin_size):
                end = min(start + a.bin_size, len(sequence))
                b = sequence[start:end].upper()
                canonical = sum(b.count(x) for x in ("A", "C", "G", "T"))
                if end - start < a.bin_size or canonical < 0.8 * len(b):
                    continue
                # Exclude centromeres plus 1-Mb flanks. No p-arm claims on acrocentrics.
                if end > cen_lo - 1_000_000 and start < cen_hi + 1_000_000:
                    continue
                arm = "p" if end <= cen_lo else "q"
                if chrom in {"chr13", "chr14", "chr15", "chr21", "chr22"} and arm == "p":
                    continue
                gc = (b.count("G") + b.count("C")) / canonical
                if not 0.25 <= gc <= 0.70:
                    continue
                bed.write(f"{chrom}\t{start}\t{end}\n")
                writer.writerow([chrom, start, end, arm, gc, 1 - canonical / len(b)])
                counts[chrom] += 1

        chrom, parts = "", []
        with a.reference.open() as f:
            for line in f:
                if line.startswith(">"):
                    if parts:
                        emit(chrom, "".join(parts))
                    chrom, parts = line[1:].split()[0], []
                elif chrom in primary:
                    parts.append(line.strip())
            if parts:
                emit(chrom, "".join(parts))
    if seen != primary or any(counts[f"chr{i}"] < 10 for i in range(1, 23)):
        raise ValueError("Reference-bin generation did not cover the complete primary genome")
    bed_temp.replace(bed_path)
    csv_temp.replace(csv_path)
    def digest(path):
        h = hashlib.sha256()
        with path.open("rb") as f:
            while chunk := f.read(8 * 1024 * 1024):
                h.update(chunk)
        return h.hexdigest()
    manifest = {"complete": True, "generated_utc": datetime.now(timezone.utc).isoformat(),
                "bin_size": a.bin_size, "reference_sha256": digest(a.reference),
                "cytoband_sha256": digest(a.cytoband), "code_sha256": digest(Path(__file__)),
                "primary_reference_lengths": lengths, "bins_per_chromosome": dict(counts),
                "rows": sum(counts.values()), "bed_sha256": digest(bed_path), "csv_sha256": digest(csv_path)}
    destination = Path(str(a.prefix) + ".manifest.json")
    temporary = Path(str(destination) + ".partial")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.replace(destination)


if __name__ == "__main__":
    main()
