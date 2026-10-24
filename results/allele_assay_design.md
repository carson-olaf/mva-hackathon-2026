# BUB1B allele assay design

Assay planning only. No rescue, readthrough or peptide detectability was measured.

Public reference: **ENST00000287598.11**, UniProt O60566; 1050 aa. Public protein comparisons: uniprot.json: exact match, no mismatches; alphafold_metadata.json: exact match, no mismatches. No patient file was read and no network request was made by the design command.

## Exact stop context

c.2210T>G changes **TTA (Leu) → TGA (stop)**. The local mutant RNA is `CUACCAGAG[UGA]AGUGCCUCU`. The +4 base is **A**, CDS position 2212; the stop-plus-four context is **UGAA**. This is not UGAC. No readthrough percentage or inserted amino acid is predicted.

## Digestion candidates

Coordinates are one-based, inclusive. Trypsin is modelled as cleavage after K/R except before P, with zero or one missed cleavage. The 7–25-aa window is a planning filter only.

| Construct | Site | Missed cuts | Residues | Peptide | Length |
|---|---:|---:|---|---|---:|
| WT | 737 | 0 | 733–752 | `SLPELSASAELCIEDRPMPK` | 20 |
| WT | 737 | 1 | 729–752 | `QLLKSLPELSASAELCIEDRPMPK` | 24 |
| WT | 737 | 1 | 733–757 | `SLPELSASAELCIEDRPMPKLEIEK` | 25 |
| WT | 1002 | 0 | 1000–1033 | `ILNANDEATVSVLGELAAEMNGVFDTTFQSHLNK` | 34 |
| WT | 1002 | 1 | 996–1033 | `FFVRILNANDEATVSVLGELAAEMNGVFDTTFQSHLNK` | 38 |
| WT | 1002 | 1 | 1000–1037 | `ILNANDEATVSVLGELAAEMNGVFDTTFQSHLNKALWK` | 38 |
| N1002K | 737 | 0 | 733–752 | `SLPELSASAELCIEDRPMPK` | 20 |
| N1002K | 737 | 1 | 729–752 | `QLLKSLPELSASAELCIEDRPMPK` | 24 |
| N1002K | 737 | 1 | 733–757 | `SLPELSASAELCIEDRPMPKLEIEK` | 25 |
| N1002K | 1002 | 0 | 1000–1002 | `ILK` | 3 |
| N1002K | 1002 | 1 | 996–1002 | `FFVRILK` | 7 |
| N1002K | 1002 | 1 | 1000–1033 | `ILKANDEATVSVLGELAAEMNGVFDTTFQSHLNK` | 34 |

N1002K creates a cut after residue 1002. Its zero-missed site-containing peptide is `ILK` (3 aa). One missed cleavage gives **`FFVRILK` (7 aa, residues 996–1002)** when the upstream R999 cut is missed, or a 34-aa product when the new K1002 cut is missed. FFVRILK is a short mutant candidate, conditional on that missed cleavage. WT site-containing peptides are 34–38 aa; none meet the 7–25-aa planning filter. Missed-cleavage abundance must be calibrated, not assumed.

An alternative, idealized E-only digestion yields a short pair:

| Construct | Residues | Candidate |
|---|---|---|
| WT | 992–1006 | `LWNKFFVRILNANDE` |
| N1002K | 992–1006 | `LWNKFFVRILKANDE` |

These 15-aa candidates need Glu-C specificity optimization, peptide standards, proteome-uniqueness searches, modification/charge-state selection and measured detection limits. They are sequence-discriminating within the two constructs only. Routine mass-only analysis does not reliably distinguish Leu from Ile. The downstream daughter peptide after K1002 is shared sequence and alone does not establish allele origin.

## Unknown residue after readthrough

The WT L737-spanning tryptic peptide is a useful sequence reference, not proof that readthrough restores Leu. If readthrough inserts Leu, this peptide is shared with the N1002K allele and cannot by itself assign allele origin. The JSON enumerates all 20 hypothetical inserted amino acids and their 0/1-missed digestion products without assigning probabilities. Insertion of Lys or Arg adds a cleavage site. Identify the product with MS/MS and standards before interpreting it as functional rescue.

## Phase and attribution

If the two variants are in trans, the full-length missense allele carries K1002. An N1002-bearing downstream product may support translation from the stop allele, but does not by itself prove traversal of codon 737, full-length integrity, or function. The 737-inserted amino acid remains unknown.

If both changes are in cis, the opposite allele can supply WT full-length protein. A K1002 product can mark downstream translation from the double-variant allele only with isoform/reinitiation controls; a peptide crossing residue 737 plus full-length and RNA controls is needed to establish readthrough.

## Proposed experiment and null calibration

1. Sequence-confirm constructs and phase where possible: WT, L737Ter alone, N1002K alone, and the cis double mutant; do not infer the patient's phase from constructs.
2. Use the exact stop context in a normalized dual reporter with sense/no-stop and vehicle controls; pair signal with RNA abundance, viability and reporter-interference controls. A short reporter does not reproduce endogenous NMD.
3. In endogenous or full-length cDNA models, quantify upstream and downstream BubR1 peptides/protein together with BUB1B transcript abundance. Keep NMD inhibition only as a mechanistic control, not an approved-drug candidate.
4. Search MS/MS for all 20 amino-acid insertion alternatives around 737, then confirm identified products with synthetic standards. K/R insertion can change tryptic cleavage. If Leu is inserted, the 737 peptide is shared by the N1002K allele and is not allele-specific. No inserted amino acid is assumed here.
5. Optimize the one-missed FFVRILK candidate and E-only 15-aa N1002/K1002 peptide pair with standards and independent proteome-uniqueness searches. Calibrate N/K response with known-ratio mixtures, blanks, carryover and replicate dilution series; report lower limit of quantification and false calls.
6. Require concordant full-length protein, checkpoint/kinetochore function and reduced chromosome-missegregation/aneuploidy without growth toxicity before calling a compound a rescue. Measure uncertainty against vehicle and sense/WT controls.

No biological signal measured. Proposed controls: vehicle versus sense/stop reporters, known WT:N1002K mixtures, blank standards and orthogonal function. Peptide absence is uninterpretable until a measured limit of quantification exists.

## Reproduce

```sh
uv run python -m mva_hackathon_2026.allele_assay --cache work/track1-cache --structure-cache work/structure
```

To prepare the two public Ensembl inputs without a private genome, run `uv run python -c 'from mva_hackathon_2026.public_data import PublicCache; from mva_hackathon_2026.track1 import gene_model; gene_model(PublicCache("work/track1-cache"), "BUB1B")'`. This bounded public-reference request does not read a VCF. The design command itself performs no network or patient-data access. Existing UniProt/AlphaFold caches add optional sequence cross-checks; missing comparison files are reported as unavailable. SHA-256 hashes are recorded; a comparison mismatch stops the command.

## Opened primary sources

- DOI [10.7554/eLife.52611](https://elifesciences.org/articles/52611). Stop identity, +4 and wider sequence context affect readthrough in tested human-cell aminoglycoside experiments; not a prediction for this allele or azithromycin.
- DOI [10.1093/hmg/ddx196](https://pmc.ncbi.nlm.nih.gov/articles/PMC5886253/); PMID 28575328. CFTR readthrough can insert multiple context-dependent amino acids with different functional effects; does not establish BUB1B insertion identities.
