# Low-pass chromosome and arm depth: executable methods

This is an exploratory relative-depth estimator. It does not measure a total MVA-cell fraction or absolute ploidy. The biological comparator is one public normal, not a certified euploid blood control. No accepted copy-number result existed when this methods document was prepared. Completed outcomes belong in `results/lowpass_methods.md`; they must not be inferred from intermediate files.

The prospective local [analysis plan](../results/lowpass_analysis_plan.md) was recorded at **2026-10-24 01:59:15.671492 UTC**, before the first chromosome-depth inference. VCF findings, downloads, sampling and an initial resource stop were already known. The primary analysis is full patient versus public normal, **100-kb bins and 1-Mb spatial blocks**. The other 17 combinations are sensitivity analyses or controls. We do not select the most favorable combination or claim familywise error control across all 18 runs.

## Inputs, provenance and data boundaries

Only Carson-approved reads were downloaded: the patient L001 paired lane, the GRCh38 analysis-set reference and four public normal FASTQs. Additional patient downloads or any new object above 1 GB require separate approval. Patient reads stay under `data/`; derived FASTQ, BAM and per-bin tables stay under ignored `work/`. No read identifiers, raw sequence, BAM or per-bin patient table is published. The original VCF, index and phenotype document are immutable.

| Input | Verified preparation evidence |
|---|---|
| Patient L001 pair | 21,502,947,327 bytes; both authenticated upstream LFS SHA256 checks passed |
| Public normal | HG00096, ENA project PRJNA686136, runs SRR13339854/SRR13339855; four upstream MD5 checks passed; 2,499,907,009 bytes |
| GRCh38 FASTA | 3,263,683,042 bytes; HTTP range boundaries and assembled length checked; local SHA256 recorded |

The reference is the [1000 Genomes GRCh38 full analysis set plus decoy and HLA](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa). A local hash establishes repeatability of these bytes. It is not an upstream cryptographic identity check: the opened upstream listing did not supply a checksum. The mapping reference contains alternate/decoy/HLA contigs and `chr` prefixes. It differs from the original VCF reference. Before inference, all 23 original primary VCF contig lengths (autosomes plus X) must agree after prefix normalization. This checks coordinate lengths; it does not prove full-reference identity.

The normal comes from [ENA PRJNA686136](https://www.ebi.ac.uk/ena/browser/view/PRJNA686136), low-pass HiSeq X Ten B-lymphocyte material. It can differ from patient blood through tissue, cell culture, library preparation and biological copy number. Matched processing reduces some technical differences; it cannot remove those confounders.

| Normal file | Bytes | Upstream ENA MD5 |
|---|---:|---|
| SRR13339854_1.fastq.gz | 630,590,234 | `73077d0317d45bd4a42b5bd953b72218` |
| SRR13339854_2.fastq.gz | 640,217,557 | `9819b776a65d400ae46e6bf4e7391d6a` |
| SRR13339855_1.fastq.gz | 615,450,703 | `41a728a8b907d0f96b4273954214d03d` |
| SRR13339855_2.fastq.gz | 613,648,515 | `42d5905549a944cfd3a2c15adc1fe3c1` |

Two-pass validation checks FASTQ syntax, mate names, equal sequence/quality lengths and complete pairing. One BLAKE2b decision per lane and pair index keeps both mates together. The target is 3.4 billion input bases, rather than an exact read count. Observed preparation results were:

| Sample | Input pairs | Input bases | Selected pairs | Selected bases | Nominal input coverage |
|---|---:|---:|---:|---:|---:|
| Patient | 136,586,642 | 40,702,819,316 | 11,412,959 | 3,401,061,782 | 1.101284x |
| Normal | 15,263,526 | 4,609,584,852 | 11,256,953 | 3,399,599,806 | 1.100810x |

Nominal coverage divides input bases by 3,088,269,832 bp (primary autosomes plus X/Y, excluding mitochondrial sequence). Mapping, duplicate removal and analysis masks reduce usable coverage. These input values do not establish usable autosomal depth. The final results must report both mapped CIGAR bases and actual MAPQ30 usable-bin depths.

## Software, resource stops and alignment design

The run uses Python/uv with the repository lockfile, minimap2 **2.31-r1302**, and samtools/htslib **1.24**. Final per-sample `complete.json` files record tool versions, input/code hashes, batch sizes, elapsed seconds and final BAM hashes. The execution summary carries source and analysis fingerprints. Cached stages require matching inputs, parameters, code and completed output hashes, rather than file existence.

Two alignment attempts were stopped before an accepted BAM existed. The full-reference index attempt crossed the **8-GB swap stop**. A full-read partitioned attempt was stopped when its scratch projection was approximately **15 GB per reference partition**. That projection is not a measured total disk peak. Neither aborted attempt supplies the analysis BAM.

The replacement maps **750,000 paired templates per read batch**. Each batch sees the complete reference in **500-million-target-base partitions**, with a **25-million-query-base batch setting**. Minimap2 `--split-prefix` merges all reference-partition hits before global MAPQ assignment and pairing. Ordinary multipart mapping without this merge is unsuitable: its MAPQ values can be wrong. Target-base size is not a hard RAM limit. Published long-read RAM benchmarks are not applied to these short reads. The [official minimap2 manual](https://lh3.github.io/minimap2/minimap2.html), the [v2.31 merge implementation](https://raw.githubusercontent.com/lh3/minimap2/v2.31/map.c), and the partition-merging paper support this design ([DOI:10.1038/s41598-019-40739-8](https://doi.org/10.1038/s41598-019-40739-8)).

Each read batch passes through name sort, `fixmate -m`, and coordinate sort. All coordinate BAMs are then merged and **duplicate removal runs globally**, including duplicates shared across batches. The implementation is [align.py](../src/mva_hackathon_2026/lowpass/align.py). Batches and reference partitions are computational units, not independent biological samples.

The host runs one native job through `~/.agents/bin/build-lock`. Supervision requires at least **42 GB free before alignment**, stops at **32 GB free** to protect the required 30-GB floor, and stops above **8 GB swap**. No unrelated user apps are closed. The alignment module requires this external supervisor; it does not impose an OS memory limit itself. Resource samples can miss instantaneous peaks. In particular, a final sampled maximum below 8 GB would not erase the earlier swap-trigger stop. Final elapsed times and sampled extrema must come from completed manifests and telemetry, not the initial runtime estimate.

Synthetic mapping checks are in [lowpass_alignment_validation.json](../results/lowpass_alignment_validation.json). A 120-pair fixture contained 80 unique pairs and 40 ambiguous-repeat pairs. Unique alignment flags/MAPQ matched full-index mapping; ambiguous reads had MAPQ0. A second 130-pair fixture included 10 duplicate pairs across read batches: global duplicate removal retained 120 pairs and the high-MAPQ alignments matched the whole-read pipeline. Ambiguous placements and proper-pair flags differed (178 versus 240 proper-pair records); those MAPQ0 reads are excluded from depth. These checks test mapping and duplicate-handling contracts, not patient mapping accuracy or resource use.

## Reproduce with authorized inputs

Commands run from the repository root. They do not download data. Supply the approved, checksum-verified files in the paths below and the GRCh38 UCSC `cytoBand.txt.gz` used by the VCF analysis. Keep a fresh attempt directory for each alignment. Do not reuse an interrupted output on existence alone.

```bash
set -euo pipefail
umask 077
uv sync --locked
mkdir -p work/lowpass data/reference
samtools faidx data/reference/GRCh38.fa
awk '$1 ~ /_alt$/ || $1 ~ /^HLA-/ {print $1}' \
  data/reference/GRCh38.fa.fai > data/reference/alt-contigs.txt
uv run python -m mva_hackathon_2026.lowpass.subsample \
  --r1 data/raw/patient/WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz \
  --r2 data/raw/patient/WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz \
  --target-bases 3400000000 --seed mva-lowpass-v1-patient --prefix work/lowpass/patient
uv run python -m mva_hackathon_2026.lowpass.subsample \
  --r1 data/raw/normal/SRR13339854_1.fastq.gz --r2 data/raw/normal/SRR13339854_2.fastq.gz \
  --r1 data/raw/normal/SRR13339855_1.fastq.gz --r2 data/raw/normal/SRR13339855_2.fastq.gz \
  --target-bases 3400000000 --seed mva-lowpass-v1-normal --prefix work/lowpass/normal
```

The exact alignment entrypoint, once the external disk/swap supervisor is active, is:

```bash
for label in normal patient; do
  ~/.agents/bin/build-lock uv run python -m mva_hackathon_2026.lowpass.align \
    --reference data/reference/GRCh38.fa \
    --r1 work/lowpass/${label}.R1.fastq.gz --r2 work/lowpass/${label}.R2.fastq.gz \
    --directory work/lowpass/${label}-fresh-attempt \
    --output work/lowpass/${label}.bam --alt data/reference/alt-contigs.txt \
    --pairs 750000 --index-size 500M --threads 6
  samtools quickcheck -v work/lowpass/${label}.bam
  samtools flagstat work/lowpass/${label}.bam > work/lowpass/${label}.flagstat.txt
  samtools stats work/lowpass/${label}.bam > work/lowpass/${label}.stats.txt
done
```

The host build lock is a local scheduling utility, not a package dependency. On another host, use an equivalent exclusive native-job scheduler and the same disk/swap stops. Do not launch this command unsupervised on a constrained host. Raw-read cleanup requires a complete manifest bound to the current BAM hash and full selected-pair count, valid BAM/index, at least 5 million post-duplicate-removal reads, mapped fraction at least 0.75 and properly paired fraction at least 0.60. These are gross QC gates, not a clinical accuracy guarantee. The complete reference sequence dictionary and independently recomputed saved QC are additional requirements before CN inference. Delete only the explicitly approved downloaded FASTQs after the cleanup checks; an output filename is insufficient deletion evidence.

### Matched comparisons and complete depth tables

Run the following Python block after BAM QC. It creates a disjoint-template normal split (`-U` complement), then matches each contrast on `samtools stats` mapped CIGAR bases. Both mates share the samtools QNAME decision. Supplementary/secondary records follow their template; the source normal is not independently sequenced twice. Fractions and actual retained bases are recorded. Mapped CIGAR bases include nonprimary contigs; equality here does not guarantee equal filtered autosomal depth.

```bash
uv run python - <<'PY'
import json
from pathlib import Path
import subprocess as sp
w = Path('work/lowpass')
def run(*args):
    sp.run([str(a) for a in args], check=True)
def bases(bam):
    text = sp.check_output(['samtools', 'stats', str(bam)], text=True)
    return next(int(r.split('\t')[2]) for r in text.splitlines()
                if r.startswith('SN\tbases mapped (cigar):\t'))
def finish(tmp, final):
    run('samtools', 'quickcheck', '-v', tmp)
    tmp.replace(final)
    run('samtools', 'index', final)
a, b = w/'normal_halfA.bam', w/'normal_halfB.bam'
ta, tb = w/'normal_halfA.bam.partial', w/'normal_halfB.bam.partial'
run('samtools', 'view', '--subsample', '.5', '--subsample-seed', '20261023',
    '-U', tb, '-b', '-o', ta, w/'normal.bam')
finish(ta, a); finish(tb, b)
count = lambda p: int(sp.check_output(['samtools', 'view', '-c', str(p)]))
assert count(w/'normal.bam') == count(a) + count(b)
records = {}
def match(left, right, label):
    totals = [bases(left), bases(right)]
    target = min(totals)
    chosen, evidence = [], []
    for role, bam, total in zip(('left', 'right'), (left, right), totals):
        fraction = target / total
        out = bam
        if fraction < .999999:
            out = w/f'{label}_{role}_matched.bam'
            tmp = Path(str(out)+'.partial')
            run('samtools', 'view', '--subsample', f'{fraction:.12g}',
                '--subsample-seed', '20261024', '-b', '-o', tmp, bam)
            finish(tmp, out)
        actual = bases(out)
        evidence.append(dict(role=role, fraction=fraction, mapped_cigar_bases=actual))
        chosen.append(str(out))
    assert .99 < evidence[0]['mapped_cigar_bases']/evidence[1]['mapped_cigar_bases'] < 1.01
    records[label] = dict(bams=chosen, sampling=evidence)
    return chosen
match(w/'patient.bam', w/'normal.bam', 'case')
control = match(a, b, 'normal_selfsplit')
match(w/'patient.bam', Path(control[0]), 'case_half')
tmp = w/'comparison_inputs.json.partial'
tmp.write_text(json.dumps(records, indent=2)+'\n')
tmp.replace(w/'comparison_inputs.json')
PY
```

Generate reference-only bins, obtain depth with samtools 1.24 `bedcov`, and analyze the primary before sensitivities. Flag mask `0xF04` excludes unmapped, secondary, QC-failed, duplicate and supplementary reads. `-Q 30` requires MAPQ30; `-j` excludes CIGAR D/N spans from depth. The [samtools bedcov manual](https://www.htslib.org/doc/samtools-bedcov.html) defines these options. Per-bin files stay in ignored `work/`.

```bash
for width in 100000 50000 250000; do
  uv run python -m mva_hackathon_2026.lowpass.bins \
    --reference data/reference/GRCh38.fa --cytoband data/reference/cytoBand.txt.gz \
    --prefix work/lowpass/bins_${width} --bin-size ${width}
done
uv run python - <<'PY'
import json
from pathlib import Path
import subprocess as sp
w = Path('work/lowpass')
contrasts = json.loads((w/'comparison_inputs.json').read_text())
for width in (100000, 50000, 250000):
    meta = w/f'bins_{width}'
    for label in ('case', 'normal_selfsplit', 'case_half'):
        name = f'{label}_{width//1000}kb'
        depth = w/f'depth_{name}.tsv'
        tmp = Path(str(depth)+'.partial')
        with tmp.open('w') as out:
            sp.run(['samtools','bedcov','-Q','30','-G','0xF04','-j',str(meta)+'.bed',
                    *contrasts[label]['bams']], stdout=out, check=True)
        with Path(str(meta)+'.bed').open() as expected, tmp.open() as actual:
            for row in expected:
                fields = actual.readline().rstrip('\n').split('\t')
                assert len(fields) == 5 and fields[:3] == row.rstrip('\n').split('\t')
                assert all(int(value) >= 0 for value in fields[3:])
            assert not actual.readline()
        tmp.replace(depth)
        for block in (1, 5):
            key = f'lowpass_{name}_{block}mb'
            sp.run(['uv','run','python','-m','mva_hackathon_2026.lowpass.analyze',
                    '--metadata',str(meta)+'.csv','--bedcov',str(depth),
                    '--work-prefix',str(w/key),'--public-prefix',f'results/{key}',
                    '--comparison',label,'--block-mb',str(block),
                    '--iterations','30000','--seed','20261024'], check=True)
PY
```

Use fresh run paths if input bytes, code or parameters change. The private controller additionally verifies content-hash stage provenance, exact reference dictionaries and source-file integrity. A completed execution is publishable only when `results/lowpass_run_status.json` has `valid_for_current_run: true` and `execution_summary_sha256` matches the exact current execution-summary bytes. The aggregate publisher verifies all 18 expected analysis summaries against that record and checks 72 unique chromosome/arm rows per table. A partial or stale table must not become a result because its filename exists.

## Estimator, confidence intervals and conditional power

Reference bins are full-width, nonoverlapping primary-chromosome intervals. Exclude bins with over 20% non-ACGT sequence, GC outside 0.25–0.70, centromeres plus 1-Mb flanks, and acrocentric p arms (13/14/15/21/22). All 24 primary reference chromosome lengths must match the cytoband annotation; every autosome must have at least 10 retained reference bins. There is no independently measured mappability track or validated CN blacklist.

Depth is covered bases divided by bin width. As implemented in the opened [samtools 1.24 bedcov source](https://raw.githubusercontent.com/samtools/samtools/1.24/bedcov.c), overlapping mate bases are counted under the same rule in both samples; the analysis does not deduplicate overlapping coverage within a fragment. Different insert-size distributions and read-pair structure can therefore affect effective sampling depth and residual library bias. The public normal insert-size histogram has median 192 bp and a long tail; equal CIGAR-base totals are not equal numbers of independent DNA molecules. Inference uses autosomes only. Normal depth must lie within its autosomal 1st–99th percentiles, with a minimum of 0.1x; patient depth must be positive. At least 1,000 bins must remain. Each sample is fitted separately with cubic `log2(depth) ~ GC`. Four robust iterations remove fit residuals beyond four scaled MADs (MAD multiplier 1.4826, floor 0.01). The patient-minus-normal residual contrast is median-centered across retained autosomal bins. Individual centered profiles are descriptive; they are not absolute CN.

Blocks group bins by chromosome, arm and floor(start/block length). Each qualifying block retains at least half its expected bins, and contributes one median regardless of bin count. At least 100 blocks overall and five blocks per screened region are required. Sex chromosomes remain descriptive; excluded arms and insufficient regions have explicit statuses. Whole-chromosome and arm tests overlap.

For each target region, the observed statistic is its median block contrast minus the pooled median on other chromosomes. Its conditional null draws a comparator chromosome median offset and adds a median of resampled within-chromosome-centered block residuals. There are at most **21 chromosome offsets**; 30,000 draws do not create 30,000 control genomes. The null assumes exchangeable chromosome offsets and blocks. It does not estimate patient/public-normal tissue or library bias.

Let J be the actual eligible test count within one run. A two-sided empirical p-value is `(1 + null exceedances)/(30000 + 1)`; the reported adjusted value is `min(1, J*p)`. The nominal threshold is the `(1 - 0.05/J)` quantile of absolute null statistics. This is conditional per-run Bonferroni calibration, not proven clinical error control or an across-18-run claim. Relative CN is `2 * 2**observed`, assuming a mostly diploid baseline. Both a within-region block-bootstrap 95% interval and an empirical-null 95% interval are reported. These are pointwise intervals, not simultaneous familywise bands; an interval can exclude 2 without passing the Bonferroni call threshold. Neither is a validated biological clone-fraction interval.

Conditional 90% one-copy detection limits use a shift of the already corrected statistic. With threshold T and null N, required gain shift is `max(0, T - q10(N))`; required loss magnitude is `max(0, T + q90(N))`. For one-copy gain in cell fraction f, ratio is `1 + f/2`; for loss it is `1 - f/2`. The inverse converts each shift to f; limits above 100% are reported as not reached. These are directional threshold-crossing approximations. Power at 2.5%, 5%, 10% and 20% modeled fractions is also evaluated by shifting null draws.

Signal injection occurs **after GC correction and filtering**. GC fitting may absorb part of a chromosome-specific change. Thus these limits assume signal preservation; they are not end-to-end validated sensitivity, a patient mosaicism upper bound or proof of no aneuploidy. Whole-genome changes, opposing gain/loss clones, tissue restriction and uncertain baseline ploidy can escape the assay. The normal self-split estimates stochastic same-library noise only. Full-normal and half-normal contrasts reuse one biological comparator and are not independent rescue evidence.

The broader shallow-WGS rationale is consistent with GC-aware depth methods, but this implementation is not QDNAseq and borrows none of its performance claims ([Scheinin et al., PMID:25236618; DOI:10.1101/gr.175141.114](https://doi.org/10.1101/gr.175141.114)). The paper identity and primary abstract were verified through Europe PMC; this does not imply full-text review. A small full-estimator synthetic test exercises null behavior, a known CN gain and 1/5-Mb block sensitivity. It cannot establish clinical sensitivity on this patient.

## Pairing-QC sensitivity

Before patient depth inference, the public normal had 15,975,720 usable Q30 primary records, of which 12,978,429 (81.2385%) had the proper-pair flag. The original analysis remains primary. One additional 100-kb/1-Mb comparison requires proper pairs and re-matches mapped CIGAR bases. It is outside the 18-run grid and supplies no independent biological replication or across-run familywise guarantee.

To reproduce it, filter each accepted BAM with `samtools view -f 2 -F 0xF04 -b`, write a fresh BAM beneath `work/`, run `samtools quickcheck` and index it. Apply the same mapped-CIGAR-base matching procedure above to these two BAMs, retaining seed 20261024 and the 1% realized matching tolerance. Run the same `bedcov -Q 30 -G 0xF04 -j` command on the existing 100-kb BED, then the same analysis CLI with `--block-mb 1 --iterations 30000 --seed 20261024`, using a fresh public prefix. Final records include retained read counts, selection fractions and usable coverage. Requiring proper pairs can alter library/mappability selection, so agreement is a robustness check, not proof that the comparator is biologically normal.

## Raw-depth injection recovery

The prospective addendum tests attenuation before GC correction. After the normal self-split completes, run:

```bash
uv run python -m mva_hackathon_2026.lowpass.recovery \
  --metadata work/lowpass/bins_100000.csv \
  --bedcov work/lowpass/depth_normal_selfsplit_100kb.tsv \
  --public-prefix results/lowpass_normal_selfsplit_pre_gc
```

For each of 22 autosomes, deterministic raw covered-base multipliers represent a single-copy gain/loss in 10%, 20% or 30% of cells. The same filtering and GC fit run again. The change in the chromosome's 1-Mb block statistic is divided by the known log2 shift. The unspiked self-split is the null baseline. The 132 aggregate rows report recovery and retained-bin changes; no positions are published. This checks correction/filter attenuation in this control. It does not validate read alignment, library confounding, clinical sensitivity or biological clone fractions. Conditional detection limits are not automatically adjusted.

## Disclosure and completion record

AI assistance: **OpenAI; ChatGPT Pro; model training/data sharing disabled**, as confirmed by Carson. Only OpenAI models were used for this repository's patient-data work. No private genome was bulk-uploaded to external annotation services. Public references cited above were opened; analysis claims must be tied to current manifests and actual outputs. No submission or push is authorized by running these commands.

Final reporting must include actual filtered coverage, matched-base fractions and counts, primary and all sensitivity outcomes, confidence intervals, conditional detection limits, completed tool/checksum provenance, sampled resource extrema with their sampling caveat, and proof that source files are unchanged and approved raw reads were removed only after accepted BAM QC. Until those checks pass, label the analysis incomplete.
