"""Map paired read batches to the entire reference, then remove duplicates globally.

Run through the host build lock and an external disk/swap guard. Read batches
bound split-index scratch space; --split-prefix preserves global MAPQ within
each batch. Neither reference partitions nor read batches are separate samples.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import subprocess
import time

from .provenance import atomic_json, fingerprint


def pipeline(commands, log):
    children = []
    previous = None
    try:
        for command in commands:
            child = subprocess.Popen([str(x) for x in command], stdin=previous,
                                     stdout=subprocess.PIPE, stderr=log)
            if previous is not None:
                previous.close()
            children.append(child)
            previous = child.stdout
        # The last command writes to an explicit file, not stdout.
        previous.close()
        statuses = [child.wait() for child in reversed(children)]
        if any(statuses):
            raise RuntimeError(f"Alignment pipeline failed: {statuses}")
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
        for child in children:
            child.wait()


def record(handle):
    first = handle.readline()
    if not first:
        return None
    rows = [first, handle.readline(), handle.readline(), handle.readline()]
    if (not all(rows) or not rows[0].startswith(b"@") or not rows[2].startswith(b"+")
            or len(rows[1].rstrip()) != len(rows[3].rstrip())):
        raise ValueError("Invalid or truncated FASTQ record; content suppressed")
    return rows


def pair_key(row):
    name = row[0].split()[0]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def batches(r1, r2, directory, pairs):
    opener = lambda p: gzip.open(p, "rb") if str(p).endswith(".gz") else Path(p).open("rb")
    with opener(r1) as left, opener(r2) as right:
        number = 0
        while True:
            a, b = record(left), record(right)
            if a is None and b is None:
                return
            number += 1
            target = directory / f"batch-{number:04d}"
            target.mkdir()
            paths = [target / "R1.fastq.gz", target / "R2.fastq.gz"]
            count = 0
            with gzip.open(paths[0], "wb", compresslevel=1) as x, gzip.open(paths[1], "wb", compresslevel=1) as y:
                while a is not None or b is not None:
                    if a is None or b is None or pair_key(a) != pair_key(b):
                        raise ValueError("FASTQ mates differ; content suppressed")
                    x.writelines(a); y.writelines(b)
                    count += 1
                    if count == pairs:
                        break
                    a, b = record(left), record(right)
            yield target, paths, count


def align(reference, r1, r2, directory, output, *, pairs=750000, index_size="500M", threads=6, alt=None):
    work_root = (Path.cwd() / "work").resolve()
    data_root = (Path.cwd() / "data").resolve()
    if any(not Path(p).resolve().is_relative_to(work_root) for p in (directory, output)):
        raise ValueError("Alignment outputs must stay beneath ignored work/")
    if any(not (Path(p).resolve().is_relative_to(work_root) or Path(p).resolve().is_relative_to(data_root))
           for p in (r1, r2)):
        raise ValueError("Read inputs must stay beneath data/ or work/")
    if pairs < 1 or threads < 1 or Path(r1).resolve() == Path(r2).resolve():
        raise ValueError("Use positive batch/thread counts and distinct mate files")
    if Path(output).exists() or Path(str(output) + ".bai").exists():
        raise FileExistsError("Refusing to replace an existing BAM or index")
    source_fingerprint = fingerprint(Path(__file__))
    directory.mkdir(parents=True, exist_ok=False)
    versions = {name: subprocess.check_output([name, "--version"], text=True).splitlines()[0]
                for name in ("minimap2", "samtools")}
    inputs = {str(p): fingerprint(Path(p)) for p in (reference, r1, r2, *([alt] if alt else []))}
    chunks, counts = [], []
    began = time.monotonic()
    for target, reads, count in batches(r1, r2, directory, pairs):
        bam = target / "coordinate.bam"
        mapping = ["minimap2", "-ax", "sr", "-I", index_size, "-K", "25M",
                   "--split-prefix", target / "merge", "--secondary=no", "-t", threads]
        if alt:
            mapping += ["--alt", alt]
        mapping += [reference, *reads]
        with (target / "alignment.log").open("w") as log:
            pipeline([mapping,
                      ["samtools", "sort", "-n", "-@", 1, "-m", "256M", "-T", target / "name", "-O", "BAM", "-"],
                      ["samtools", "fixmate", "-m", "-u", "-", "-"],
                      ["samtools", "sort", "-@", 1, "-m", "256M", "-T", target / "coord", "-o", bam, "-"]], log)
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
        for path in reads:
            path.unlink()
        # Minimap2 normally unlinks its split files at successful merge.
        for path in target.glob("merge*.tmp"):
            path.unlink()
        chunks.append(bam); counts.append(count)
        atomic_json(directory / "progress.json", {"completed_batches": len(chunks), "pairs": sum(counts),
                                                   "elapsed_seconds": time.monotonic() - began})
    if not chunks:
        raise ValueError("No read pairs")
    listing = directory / "chunks.txt"
    listing.write_text("".join(str(p) + "\n" for p in chunks))
    temporary = directory / "aligned.bam"
    with (directory / "merge.log").open("w") as log:
        pipeline([["samtools", "merge", "-u", "-b", listing, "-"],
                  ["samtools", "markdup", "-r", "-@", 1, "-T", directory / "duplicates", "-", temporary]], log)
    subprocess.run(["samtools", "quickcheck", "-v", str(temporary)], check=True)
    subprocess.run(["samtools", "index", "-@", "2", str(temporary)], check=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary.replace(output)
    Path(str(temporary) + ".bai").replace(Path(str(output) + ".bai"))
    result = {"complete": True, "versions": versions, "input_fingerprints": inputs,
              "pairs_per_batch": counts, "total_pairs": sum(counts), "index_batch": index_size,
              "query_batch": "25M", "threads": threads, "global_duplicate_removal": True,
              "split_prefix_global_merge_per_batch": True, "elapsed_seconds": time.monotonic() - began,
              "bam_fingerprint": fingerprint(output), "code_fingerprint": source_fingerprint}
    atomic_json(directory / "complete.json", result)
    for bam in chunks:
        bam.unlink()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("reference", "r1", "r2", "directory", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--alt", type=Path)
    parser.add_argument("--pairs", type=int, default=750000)
    parser.add_argument("--index-size", default="500M")
    parser.add_argument("--threads", type=int, default=6)
    args = vars(parser.parse_args())
    if args["pairs"] < 1 or args["threads"] < 1:
        parser.error("pairs and threads must be positive")
    result = align(**args)
    print(json.dumps({k: result[k] for k in ("complete", "total_pairs", "elapsed_seconds")}))


if __name__ == "__main__":
    main()
