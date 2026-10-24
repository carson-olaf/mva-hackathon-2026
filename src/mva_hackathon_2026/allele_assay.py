"""Offline assay planning from public BUB1B reference sequence and named variants.

No patient file is read. Digestion is an explicit sequence rule, not a prediction
of cleavage efficiency, peptide uniqueness, ionization, or biological rescue.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from .public_data import PublicCache
from .track1 import gene_model, translate

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class OfflinePublicCache(PublicCache):
    """Reuse the Track 1 public reference cache; never fetch implicitly."""

    def get(self, url, *, params=None, text=False):
        key = hashlib.sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        path = self.root / f"{key}.{'txt' if text else 'json'}"
        if not path.exists():
            raise FileNotFoundError(f"Missing public reference cache: {path}")
        return super().get(url, params=params, text=text)


def digest(sequence: str, *, cleavage_after="KR", block_before_proline=True,
           max_missed=1) -> list[dict]:
    """Enumerate idealized products, with one-based inclusive coordinates."""
    if max_missed < 0 or any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("Use a protein sequence without stop/ambiguous residues")
    cuts = [0]
    for i, aa in enumerate(sequence):
        if aa in cleavage_after and (not block_before_proline or
                i + 1 == len(sequence) or sequence[i + 1] != "P"):
            cuts.append(i + 1)
    if cuts[-1] != len(sequence):
        cuts.append(len(sequence))
    products = []
    for missed in range(max_missed + 1):
        for i in range(len(cuts) - 1 - missed):
            start, end = cuts[i], cuts[i + missed + 1]
            products.append({"start": start + 1, "end": end,
                             "sequence": sequence[start:end], "length": end - start,
                             "missed_cleavages": missed})
    return products


def site_peptides(sequence: str, position: int, **kwargs) -> list[dict]:
    return [{**p, "within_planning_length_7_to_25": 7 <= p["length"] <= 25,
             "spans_both_sides_of_site": p["start"] < position < p["end"]}
            for p in digest(sequence, **kwargs) if p["start"] <= position <= p["end"]]


def _sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def build_design(cache_root: str | Path, structure_root: str | Path) -> dict:
    model = gene_model(OfflinePublicCache(cache_root), "BUB1B")
    cds = model["sequence"]
    translated = translate(cds)
    if len(cds) % 3 or not translated.endswith("*") or "*" in translated[:-1]:
        raise ValueError("Reference CDS lacks the expected single terminal stop")
    protein = translated[:-1]
    if protein[736] != "L" or protein[1001] != "N" or cds[2209] != "T" or cds[3005] != "T":
        raise ValueError("Public reference does not match the named coding variants")
    comparisons = []
    for filename, extractor in (
        ("uniprot.json", lambda obj: obj["sequence"]["value"]),
        ("alphafold_metadata.json", lambda obj: obj[0]["uniprotSequence"]),
    ):
        path = Path(structure_root) / filename
        if path.exists():
            other = extractor(json.loads(path.read_text()))
            mismatches = [i + 1 for i, (a, b) in enumerate(zip(protein, other)) if a != b]
            comparisons.append({"source_file": str(path), "file_sha256": _sha(path.read_text()),
                                "length": len(other), "protein_sha256": _sha(other),
                                "exact_match": other == protein,
                                "substitution_mismatch_positions": mismatches,
                                "length_difference": len(other) - len(protein)})
        else:
            comparisons.append({"source_file": str(path), "status": "not_available"})
    if any(c.get("exact_match") is False for c in comparisons):
        raise ValueError("Public reference protein sources disagree; inspect before design")
    stop_cds = cds[:2209] + "G" + cds[2210:]
    lys_cds = cds[:3005] + "G" + cds[3006:]
    lys_protein = translate(lys_cds)[:-1]
    if stop_cds[2208:2211] != "TGA" or lys_protein[1001] != "K":
        raise ValueError("Named variant translation check failed")
    context = lambda seq: seq[2199:2208] + "[" + seq[2208:2211] + "]" + seq[2211:2220]
    sequences = {"WT": protein, "N1002K": lys_protein}
    tryptic = {name: {str(pos): site_peptides(seq, pos) for pos in (737, 1002)}
               for name, seq in sequences.items()}
    alternatives = {name: site_peptides(seq, 1002, cleavage_after="E",
                                       block_before_proline=False, max_missed=0)
                    for name, seq in sequences.items()}
    # Exhaustive planning alternatives, with no claim about inserted-AA frequency.
    readthrough = []
    for aa in AMINO_ACIDS:
        hypothetical = protein[:736] + aa + protein[737:]
        readthrough.append({"hypothetical_inserted_amino_acid": aa,
                            "is_prediction": False,
                            "tryptic_site_peptides": site_peptides(hypothetical, 737)})
    return {
        "status": "assay_planning_only_no_experimental_rescue_evidence",
        "reference": {"gene": "BUB1B", "uniprot": "O60566", "transcript": model["transcript"],
                      "cds_length_nt_including_terminal_stop": len(cds), "protein_length_aa": len(protein),
                      "cds_sha256": _sha(cds), "protein_sha256": _sha(protein),
                      "public_reference_comparisons": comparisons,
                      "cache_root": str(cache_root), "network_used": False,
                      "patient_files_read": False},
        "L737Ter_context": {"coding_variant": "c.2210T>G", "wild_type_dna_codon": cds[2208:2211],
                            "mutant_dna_codon": stop_cds[2208:2211],
                            "wild_type_dna_context_9nt_each_side": context(cds),
                            "mutant_rna_context_9nt_each_side": context(stop_cds).replace("T", "U"),
                            "stop_rna": "UGA", "plus4_base": stop_cds[2211],
                            "plus4_cds_position_1based": 2212,
                            "stop_plus4_rna": stop_cds[2208:2212].replace("T", "U"),
                            "readthrough_rate": None, "inserted_amino_acid": None,
                            "interpretation": "Exact context only; neither readthrough rate nor inserted amino acid is inferred."},
        "N1002K_context": {"coding_variant": "c.3006T>G", "wild_type_dna_codon": cds[3003:3006],
                             "mutant_dna_codon": lys_cds[3003:3006],
                             "new_trypsin_cut_after_residue": 1002,
                             "tryptic_site_marker_7_to_25aa_available_by_construct": {name: any(7 <= p["length"] <= 25
                                for p in site_peptides(seq, 1002)) for name, seq in sequences.items()}},
        "digestion_rules": {"trypsin": "After K or R unless the next residue is P; 0 and 1 missed cleavages.",
                            "alternative": "Idealized Glu-C E-only rule: after E, no proline exclusion, 0 missed cleavages; specificity is condition-dependent and needs optimization.",
                            "length_window": "7–25 aa is a planning filter, not a measured detection threshold.",
                            "modifications": "None modelled; methionine oxidation/cysteine alkylation and charge-state transitions require laboratory method development."},
        "tryptic_site_peptides": tryptic,
        "alternative_E_only_site1002_peptides": alternatives,
        "readthrough_737_all20_hypothetical_amino_acid_products": readthrough,
        "phase_interpretation": {
            "phase_proven": False,
            "trans": "If the two variants are in trans, the full-length missense allele carries K1002. An N1002-bearing downstream product may support translation from the stop allele, but does not by itself prove traversal of codon 737, full-length integrity, or function. The 737-inserted amino acid remains unknown.",
            "cis": "If both changes are in cis, the opposite allele can supply WT full-length protein. A K1002 product can mark downstream translation from the double-variant allele only with isoform/reinitiation controls; a peptide crossing residue 737 plus full-length and RNA controls is needed to establish readthrough.",
            "sequence_discrimination": "Candidate sequences distinguish the two specified constructs. Proteome uniqueness, detection, allele expression and phase are not established. Leu/Ile cannot generally be resolved by routine mass-only peptide assignment."},
        "experiment": [
            "Sequence-confirm constructs and phase where possible: WT, L737Ter alone, N1002K alone, and the cis double mutant; do not infer the patient's phase from constructs.",
            "Use the exact stop context in a normalized dual reporter with sense/no-stop and vehicle controls; pair signal with RNA abundance, viability and reporter-interference controls. A short reporter does not reproduce endogenous NMD.",
            "In endogenous or full-length cDNA models, quantify upstream and downstream BubR1 peptides/protein together with BUB1B transcript abundance. Keep NMD inhibition only as a mechanistic control, not an approved-drug candidate.",
            "Search MS/MS for all 20 amino-acid insertion alternatives around 737, then confirm identified products with synthetic standards. K/R insertion can change tryptic cleavage. If Leu is inserted, the 737 peptide is shared by the N1002K allele and is not allele-specific. No inserted amino acid is assumed here.",
            "Optimize the one-missed FFVRILK candidate and E-only 15-aa N1002/K1002 peptide pair with standards and independent proteome-uniqueness searches. Calibrate N/K response with known-ratio mixtures, blanks, carryover and replicate dilution series; report lower limit of quantification and false calls.",
            "Require concordant full-length protein, checkpoint/kinetochore function and reduced chromosome-missegregation/aneuploidy without growth toxicity before calling a compound a rescue. Measure uncertainty against vehicle and sense/WT controls."],
        "null_calibration": "No biological signal measured. Proposed controls: vehicle versus sense/stop reporters, known WT:N1002K mixtures, blank standards and orthogonal function. Peptide absence is uninterpretable until a measured limit of quantification exists.",
        "references": [
            {"doi": "10.7554/eLife.52611", "url": "https://elifesciences.org/articles/52611", "opened_scope": "primary full HTML", "supports": "Stop identity, +4 and wider sequence context affect readthrough in tested human-cell aminoglycoside experiments; not a prediction for this allele or azithromycin."},
            {"pmid": "28575328", "doi": "10.1093/hmg/ddx196", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5886253/", "opened_scope": "primary full HTML", "supports": "CFTR readthrough can insert multiple context-dependent amino acids with different functional effects; does not establish BUB1B insertion identities."}]
    }


def markdown(design: dict) -> str:
    ref, stop = design["reference"], design["L737Ter_context"]
    lines = ["# BUB1B allele assay design", "", "Assay planning only. No rescue, readthrough or peptide detectability was measured.", "",
             f"Public reference: **{ref['transcript']}**, UniProt O60566; {ref['protein_length_aa']} aa. Public protein comparisons: " + "; ".join(f"{Path(c['source_file']).name}: " + ("exact match, no mismatches" if c.get("exact_match") else "unavailable") for c in ref["public_reference_comparisons"]) + ". No patient file was read and no network request was made by the design command.", "",
             "## Exact stop context", "", f"c.2210T>G changes **TTA (Leu) → TGA (stop)**. The local mutant RNA is `{stop['mutant_rna_context_9nt_each_side']}`. The +4 base is **A**, CDS position 2212; the stop-plus-four context is **UGAA**. This is not UGAC. No readthrough percentage or inserted amino acid is predicted.", "",
             "## Digestion candidates", "", "Coordinates are one-based, inclusive. Trypsin is modelled as cleavage after K/R except before P, with zero or one missed cleavage. The 7–25-aa window is a planning filter only.", "", "| Construct | Site | Missed cuts | Residues | Peptide | Length |", "|---|---:|---:|---|---|---:|"]
    for name, sites in design["tryptic_site_peptides"].items():
        for site, peptides in sites.items():
            for p in peptides:
                lines.append(f"| {name} | {site} | {p['missed_cleavages']} | {p['start']}–{p['end']} | `{p['sequence']}` | {p['length']} |")
    lines += ["", "N1002K creates a cut after residue 1002. Its zero-missed site-containing peptide is `ILK` (3 aa). One missed cleavage gives **`FFVRILK` (7 aa, residues 996–1002)** when the upstream R999 cut is missed, or a 34-aa product when the new K1002 cut is missed. FFVRILK is a short mutant candidate, conditional on that missed cleavage. WT site-containing peptides are 34–38 aa; none meet the 7–25-aa planning filter. Missed-cleavage abundance must be calibrated, not assumed.", "", "An alternative, idealized E-only digestion yields a short pair:", "", "| Construct | Residues | Candidate |", "|---|---|---|"]
    for name, peptides in design["alternative_E_only_site1002_peptides"].items():
        for p in peptides:
            lines.append(f"| {name} | {p['start']}–{p['end']} | `{p['sequence']}` |")
    lines += ["", "These 15-aa candidates need Glu-C specificity optimization, peptide standards, proteome-uniqueness searches, modification/charge-state selection and measured detection limits. They are sequence-discriminating within the two constructs only. Routine mass-only analysis does not reliably distinguish Leu from Ile. The downstream daughter peptide after K1002 is shared sequence and alone does not establish allele origin.", "", "## Unknown residue after readthrough", "", "The WT L737-spanning tryptic peptide is a useful sequence reference, not proof that readthrough restores Leu. If readthrough inserts Leu, this peptide is shared with the N1002K allele and cannot by itself assign allele origin. The JSON enumerates all 20 hypothetical inserted amino acids and their 0/1-missed digestion products without assigning probabilities. Insertion of Lys or Arg adds a cleavage site. Identify the product with MS/MS and standards before interpreting it as functional rescue.", "", "## Phase and attribution", "", design["phase_interpretation"]["trans"], "", design["phase_interpretation"]["cis"], "", "## Proposed experiment and null calibration", ""]
    lines += [f"{i}. {step}" for i, step in enumerate(design["experiment"], 1)]
    lines += ["", design["null_calibration"], "", "## Reproduce", "", "```sh", "uv run python -m mva_hackathon_2026.allele_assay --cache work/track1-cache --structure-cache work/structure", "```", "", "To prepare the two public Ensembl inputs without a private genome, run `uv run python -c 'from mva_hackathon_2026.public_data import PublicCache; from mva_hackathon_2026.track1 import gene_model; gene_model(PublicCache(\"work/track1-cache\"), \"BUB1B\")'`. This bounded public-reference request does not read a VCF. The design command itself performs no network or patient-data access. Existing UniProt/AlphaFold caches add optional sequence cross-checks; missing comparison files are reported as unavailable. SHA-256 hashes are recorded; a comparison mismatch stops the command.", "", "## Opened primary sources", ""]
    for r in design["references"]:
        lines += [f"- DOI [{r['doi']}]({r['url']})" + (f"; PMID {r['pmid']}" if "pmid" in r else "") + f". {r['supports']}"]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="work/track1-cache")
    parser.add_argument("--structure-cache", default="work/structure")
    parser.add_argument("--output", default="results/allele_assay_design")
    args = parser.parse_args()
    design = build_design(args.cache, args.structure_cache)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(json.dumps(design, indent=2) + "\n")
    output.with_suffix(".md").write_text(markdown(design))
    print(json.dumps({"transcript": design["reference"]["transcript"],
                      "stop_plus4": design["L737Ter_context"]["stop_plus4_rna"],
                      "status": design["status"], "outputs": [str(output.with_suffix(s)) for s in (".json", ".md")]}))


if __name__ == "__main__":
    main()
