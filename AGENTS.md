# MVA Hackathon 2026: agent rules

Track 2 (drug repurposing) entry for the Sage Bionetworks "Rare Disease, Real Kid" hackathon, plus a Track 1 CSV for authorship. Deadline: 2026-10-24 23:59 UTC. Start with `private/findings.md`, then the task prompt you were given.

## Data rules (from the official rules and organizer guidance; breaking them can disqualify Carson)

- Patient data lives only in `data/` (a link to `~/mva-hackathon-2026-data`, mode 700, excluded from Time Machine). It holds the VCF (GRCh38, contigs named `1`..`22`, `X`, no `chr`), its index and the clinical phenotype .docx. The 85 GB of FASTQ files were not downloaded. Ask Carson before downloading any of them; if he approves, they go into `data/` only, and at least 30 GB of disk must stay free.
- Never commit, copy into tracked files, upload or paste: the VCF, any slice or reformat of it, BAM/CRAM/FASTQ, or any table of per-variant genotypes at genotype scale. Write intermediates to `work/` (git-ignored) only.
- Everything you print to your own context goes to the model provider. Print summaries and counts, not raw VCF blocks. A handful of named variants is a finding and is fine.
- External APIs (Ensembl VEP REST, ClinVar, gnomAD, PharmGKB/CPIC, Open Targets, PubMed) are allowed as tools, but send only a small named-variant list or gene/drug names. Never bulk-upload the genome.
- Do not try to contact the child, the family or MVA Society contacts. Do not write details that could re-identify them.
- Use only OpenAI models in this folder (`gpt-6-astra` or `gpt-6.1-sol`, both logged as provider `openai` by the local codex-router). Never switch to a router model from another provider: their terms are unverified for this data.
- Required AI disclosure line for the methods: provider, plan or tier, and training opt-out setting.

## Track 2 rules

- Drugs must be market-approved somewhere. Combinations of approved drugs are allowed if each drug and the combination have a rationale. Investigational compounds may appear only in a clearly separate "beyond the rules" note.
- Panel weights: scientific rigor 35%, potential impact 25%, innovation 25%, scalability 15%.
- Deliverables: report (Markdown or PDF) that characterizes the variant mechanism (LoF/GoF, pathway, downstream consequence) as the basis for the drug rationale; this repo (public after close); 3-minute pitch video.
- Every citation must be a source you opened, with PMID or DOI. Never invent references, numbers or approvals.

## Repo conventions

- Python 3.12 with uv (`uv run ...`, `uv add ...`). bcftools and tabix are on PATH (Homebrew).
- Code under `src/`, tests under `tests/` (`uv run pytest`). Keep tests few and behavioral.
- Outputs that are safe to publish go to `results/`. Check each one against the data rules before writing it.
