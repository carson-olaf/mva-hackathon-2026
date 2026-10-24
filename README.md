# MVA Hackathon 2026: test uridine triacetate for chromosome fidelity

Computational Track 2 entry and reproducible Track 1 discovery for Sage Bionetworks' **Rare Disease, Real Kid: MVA Hackathon 2026**. The aim is a testable BUB1B-MVA treatment strategy and a public-data engine for other rare diseases. These are research measurements and hypotheses, not a treatment recommendation.

**Entry proposal:** test approved oral **uridine triacetate (Xuriden)** for a nucleotide-responsive component of BUB1B-MVA. The verified label mean plasma uridine peak is 88.7 µM, versus a 30 µM **nominal** methods comparator; uridine-alone efficacy is unproved. Published DNA-damage findings at 10 or 30 µM make the paired rescue/harm experiment essential. The retrieved Burrell supplement does not establish single-component sufficiency or four-component necessity. This is a laboratory hypothesis, not a clinical treatment selection. [Decision and sources](report/nucleoside_decision.md).

## Start here (judges)

| What | Where |
|---|---|
| Track 2 report | [report/track2_report.md](report/track2_report.md) (PDF: `report/carcango_track2_report.pdf`; this public copy omits individual PGx marker genotypes that appear in the submitted version) |
| 3-minute pitch | [YouTube](https://youtu.be/DIucVnptRfM); narration in [report/pitch_script_v2.md](report/pitch_script_v2.md) |
| One-page summary | [results/codex-summary.md](results/codex-summary.md) |
| Track 1 predictions and methods | [results/track1_carcango_hpo_recessive.csv](results/track1_carcango_hpo_recessive.csv), [report/carcango_track1_methods.pdf](report/carcango_track1_methods.pdf) |
| Innovation: target-identity collision census | [report/collision_audit.md](report/collision_audit.md), [results/collision_audit/](results/collision_audit/) |
| Scalability: reusable public-data engine | [examples/README.md](examples/README.md), `src/mva_repurpose/` |

**Privacy.** No patient reads, VCFs or genotype-scale tables are published. Beyond the diagnosed BUB1B pair, individual genotype results (pharmacogenomic markers and additional gene checks) are withheld to protect the child's privacy; only aggregate counts appear. Running the patient pipeline requires the controlled-access hackathon dataset.

## Results

- Gene-blind Track 1 discovery recovers BUB1B Leu737Ter/Asn1002Lys first. Phase is unproven. Leu737Ter ends 746 nt before the final exon junction: NMD is predicted, not measured.
- 1,779,707 stringent autosomal heterozygous SNVs yield no familywise chromosome/arm outlier. Conditional 90%-power whole-chromosome limits: 7.8–10.3% gain; 7.3–9.3% loss. This cannot exclude heterogeneous MVA clones or estimate total aneuploid-cell burden.
- No rare canonical coding modifier in 49 genes. Seven observed PGx marker genotypes; no full star diplotypes. Zero exact qualifying ClinVar matches in 84 ACMG SF v3.3 genes.
- Two stability methods give no convincing support for a large Asn1002Lys stability effect; benchmark RMSE is not a variant confidence interval, and neutrality or failure of stabilization-based rescue is not established. No nearby pocket passed the exploratory gate; docking stopped. A [sequence-based assay design](results/allele_assay_design.md) verifies UGAA stop context and proposes allele-discriminating peptides with explicit cleavage/phase limits.
- LINCS source QC fails: all eight BUB1B CRISPR vectors contradict provider hashes; drug-consensus reference hashes are unavailable. The pipeline stops before inference. Verified shRNA signatures agree poorly. No drug-evidence bonus or negative efficacy claim is made.
- The earlier Pro shortlist and its scores remain an auditable historical hypothesis set; no drug has a positive BUB1B rescue result. The historical [frontier decision](report/frontier_decision.md) preserves the prior findings; inconclusive routes are deprioritized under the available evidence. The current entry advances uridine triacetate as a falsifiable laboratory hypothesis; the earlier drugs are not co-leads. [Curated replay and score limits](examples/bub1b_curated_readme.md) retain the earlier judgments and missing bridges.
- The [additional functional review](report/renewed_frontier_review.md) checks MEK bypass and nucleoside rescue without transferring results between compounds, assays or diseases. The [collision census](report/collision_audit.md) confirms 13,146 historical wrong gene–ligand pairs. The exact numeric-import bug is fixed upstream; 11 other current identity conflicts remain unresolved. Three additional gene checks found no new coding explanation; their individual calls are withheld from this public repository to protect the child's privacy.
- The public engine abstains from a BUB1B rescue claim. Its positive retrieval control ranks everolimus first for TSC2/TSC without supplying the drug name. This validates retrieval and indication matching, not novel treatment prediction.
- Direct low-pass sequencing is complete: usable autosomal-bin depth is 0.917× patient / 0.801× public normal. The prespecified 100-kb/1-Mb analysis finds no outlier among 61 chromosome/arm tests; all 17 grid sensitivity/control runs also have no outlier. These 18 runs are robustness checks, not biological replications. An additional proper-pair check flags **chr8p** (relative CN 2.0875; within-run adjusted p=0.03457), versus p=0.3192 in the primary analysis. This is an unconfirmed, method-sensitive contrast, not a clone call. Conditional primary whole-chromosome 90%-power limits are 7.07–12.40% gain / 7.77–12.48% loss; these are not a clinical exclusion or total MVA burden. [Results and calibration](results/lowpass_methods.md) retain comparator/library limitations.

Read [the evidence summary](results/codex-summary.md), [draft report](report/track2_report.md), and [Track 1 CSV](results/track1_carcango_hpo_recessive.csv). The official-ready CSV filename includes the account/team and HPO/recessive approach; `track1_carson.csv` is retained as the identical original output. The Pro handoff is integrated, with all 51 reference identities checked and access/correction gaps retained.

## Install and public replay

Python 3.12+ and `uv` are required. The public engine requires no patient genome.

```sh
uv sync
uv run pytest
uv run mva-repurpose --gene BUB1B --mechanism loss_of_function \
  --variant 'ENST00000287598.11:c.2210T>G' \
  --variant 'ENST00000287598.11:c.3006T>G' \
  --hpo HP:0002859 --hpo HP:0000121 --hpo HP:0004322 \
  --hpo HP:0001508 --hpo HP:0003202 --hpo HP:0001622 \
  --hpo HP:0001518 \
  --cache examples/public_api_cache --offline --output work/replay/bub1b
uv run mva-repurpose --gene TSC2 --mechanism loss_of_function \
  --disease 'tuberous sclerosis' --hpo HP:0001250 \
  --cache examples/public_api_cache --offline --output work/replay/positive-control
```

Omit `--offline` and use a fresh cache directory for live APIs. The committed example cache contains public references only. HPO lookup is local; phenotype combinations are not sent to APIs. See [engine usage and scoring](examples/README.md). PGx and literature bridges are optional inputs. Scores are review priorities, not efficacy probabilities. HPO terms currently provide context, not a validated drug score. ChEMBL phase 4 is an approval proxy; current jurisdiction, route, indication and pediatric suitability require separate verification.

## Authorized private-genome reruns

Install `bcftools`; source data must already exist in ignored `data/`. These commands fetch public references, not patient reads. External patient queries are limited to small named-variant shortlists.

```sh
uv run python -m mva_hackathon_2026.track1 --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --team carcango_hpo_recessive
uv run python -m mva_hackathon_2026.mosaic --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --qc strict --prefix results/mosaic_strict
uv run python -m mva_hackathon_2026.clinical_scans --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --step all
uv run python -m mva_hackathon_2026.signature --offline --shrna-qc --control-genes 100 --workers 4
```

Signature replay requires the ignored public cache from a first live run. The full drug workflow intentionally exits with code 2 at its source-quality gate on the audited inputs; historical ranking numbers are retained as unvalidated diagnostics. `--shrna-qc` runs only the checksum-verified knockdown reproducibility analysis. External structure tools and prediction-service constraints are documented in [structure methods](results/structure_methods.md). The public [low-pass methods](report/lowpass_methods.md) contain executable commands and guards; the original authorization plan remains in ignored `work/lowpass-plan.md`. Tests and installation never download patient data.

## Data handling

Patient source files remain in `data/`, a link to a mode-700 directory. Local `work/` and `private/` directories are also mode 700. Candidate tables, reads, alignments, per-bin depth and other patient-scale intermediates remain in ignored `work/` or `data/`. No genotype-scale files are committed. `results/` contains aggregates, a small set of named findings and public-drug analyses. The public engine and tests run without the genome.

No genome is uploaded. APIs receive public gene/drug identifiers, public reference structures or bounded named variants. All analysis agents use OpenAI models. Carson confirmed **Pro; model training and data sharing disabled**. No push, publication or submission is authorized. Official rules require deleting underlying data, intermediates and derived datasets within 30 days of close; subsequent deletion confirmation and public release need human review.

Reference manifests identify hashes and retrieval times. Live references may drift; the frozen public cache supports repeatable engine tests. No-hit research screens are not negative clinical tests. Phase, RNA/NMD, protein function, disease-cell aneuploidy and drug rescue require experimental confirmation.

Genome-analysis cutoff: **2026-10-24 18:00 UTC**; later public-literature follow-up and packaging were separately authorized. Submission deadline: **2026-10-24 23:59 UTC**. The report remains a draft. The three-minute pitch video is a separate deliverable.

## License

Original submission materials use [CC BY 4.0](LICENSE.md), implementing Official Rules item 10. Restricted patient data are excluded. Third-party code, references and publications retain their source terms. The required acknowledgement appears in the draft report. Public release remains a separate authorized step.
