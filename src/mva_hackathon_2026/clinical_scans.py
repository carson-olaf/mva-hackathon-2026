"""Bounded research modifier, key-marker PGx and ACMG/ClinVar scans.

Patient-scale intermediates remain in work/. Missing VCF records are unknown.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import urllib.parse
import xml.etree.ElementTree as ET
import zipfile

import requests

from .public_data import PublicCache
from .track1 import gene_model, extract_gene, IMPACT

MODIFIER_PANELS = {
    "checkpoint_CIN": "BUB1 BUB3 MAD2L1 MAD1L1 CDC20 TRIP13 CEP57 TTK KNL1 ZWINT PLK1 AURKB CENPE TP53 BUB1B".split(),
    "rhabdomyosarcoma_predisposition_or_pathway": "DICER1 NF1 HRAS KRAS NRAS FGFR4 PTCH1 SUFU RB1 CDKN2A PAX3 PAX7 FOXO1".split(),
    "nephrocalcinosis": "SLC34A1 SLC34A3 SLC9A3R1 CASR CLDN16 CLDN19 CYP24A1 CLCN5 OCRL AGXT HOGA1 GRHPR ATP6V1B1 ATP6V0A4 SLC4A1 WDR72 FAM20A SLC12A1 KCNJ1 BSND CLCNKB".split(),
}
PGX_MARKERS = {
    "CYP2C9": ["rs1799853", "rs1057910", "rs28371686", "rs9332131", "rs7900194", "rs28371685"],
    "CYP2C19": ["rs4244285", "rs4986893", "rs12248560"],
    "CYP3A5": ["rs776746", "rs10264272", "rs41303343"],
    "CYP3A4": ["rs35599367"],
    "TPMT": ["rs1800462", "rs1800460", "rs1142345"],
    "NUDT15": ["rs116855232", "rs186364861", "rs746071566"],
    "DPYD": ["rs3918290", "rs55886062", "rs67376798", "rs75017182", "rs56038477"],
    "SLCO1B1": ["rs4149056", "rs2306283"],
    "UGT1A1": ["rs3064744", "rs4148323", "rs887829"],
    "G6PD": ["rs1050828", "rs1050829", "rs5030868"],
    "ABCB1": ["rs1045642", "rs2032582", "rs1128503"],
    "CYP2D6": ["rs3892097"],
}
DRUG_RELEVANCE = {
    "CYP2C9": "warfarin; selected NSAIDs; fluvastatin (drug-specific CPIC guidance)",
    "CYP2C19": "clopidogrel; PPIs; selected antidepressants (supportive-care relevance)",
    "CYP3A5": "tacrolimus (CPIC); sirolimus/everolimus exposure is not a validated extension of that CPIC dose rule",
    "CYP3A4": "sirolimus/everolimus and other CYP3A substrates: exposure hypothesis; no CPIC genotype dose rule assigned",
    "TPMT": "azathioprine/mercaptopurine/thioguanine toxicity if considered; not an MVA efficacy rationale",
    "NUDT15": "azathioprine/mercaptopurine/thioguanine toxicity if considered; not an MVA efficacy rationale",
    "DPYD": "fluorouracil/capecitabine toxicity if considered; not the usual rhabdomyosarcoma regimen",
    "SLCO1B1": "statins including pravastatin/simvastatin: drug-specific myopathy/exposure guidance",
    "UGT1A1": "atazanavir (CPIC); irinotecan toxicity relevance is a separate oncology labeling/guideline question",
    "G6PD": "oxidant drugs including rasburicase; confirm enzyme activity where clinically indicated",
    "ABCB1": "vincristine and CYP3A-substrate transport hypotheses; no validated CPIC ABCB1 dose rule assigned",
    "CYP2D6": "codeine/tramadol and selected supportive-care drugs; SNV-only call cannot resolve CYP2D6 CNV/hybrids",
}
SF_URL = "https://search.clinicalgenome.org/kb/reports/acmg-activity-summary-report"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
SF_EXTRA = {"SMAD3", "TSC1", "TSC2"}  # Present in the paper, absent from current curation CSV.
AR_SF = {"ATP7B", "BTD", "CASQ2", "CYP27A1", "GAA", "MUTYH", "RPE65", "TRDN"}
SF_STARS = {"criteria_provided,_multiple_submitters,_no_conflicts": 2,
            "reviewed_by_expert_panel": 3, "practice_guideline": 4}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
GENE_ALIASES = {"SLC9A3R1": "NHERF1"}  # Current HGNC name, verified NCBI Gene9368.
PGX_PMIDS = ("41618934", "29152729", "36049896", "26417955", "25801146", "35152405")


def download_public(url, path):
    path = Path(path)
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, timeout=120, stream=True) as response:
        response.raise_for_status()
        if int(response.headers.get("content-length", 0)) > 1_000_000_000:
            raise ValueError("Public download exceeds 1GB approval boundary")
        temp, size = path.with_suffix(path.suffix + ".part"), 0
        try:
            with temp.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    size += len(chunk)
                    if size > 1_000_000_000:
                        raise ValueError("Public download exceeds 1GB approval boundary")
                    handle.write(chunk)
            temp.replace(path)
        finally:
            temp.unlink(missing_ok=True)
    return path


class ReferenceCache(PublicCache):
    """Read prior Track1 public cache; write new responses only to clinical cache."""
    def get(self, url, *, params=None, text=False):
        key = hashlib.sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        previous = Path("work/track1-cache") / f"{key}.{'txt' if text else 'json'}"
        if previous.exists():
            return previous.read_text() if text else json.loads(previous.read_text())
        return super().get(url, params=params, text=text)


def workbook_rows(path):
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in root.findall("m:si", NS)]
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in root.findall(".//m:row", NS):
            values = {}
            for cell in row.findall("m:c", NS):
                column = re.match("[A-Z]+", cell.attrib["r"])[0]
                value = cell.find("m:v", NS)
                if cell.attrib.get("t") == "s" and value is not None:
                    value = shared[int(value.text)]
                elif cell.attrib.get("t") == "inlineStr":
                    value = "".join(cell.find("m:is", NS).itertext())
                else:
                    value = value.text if value is not None else ""
                values[column] = value
            rows.append(values)
    return rows


def parse_cpic_markers(path, requested):
    rows = workbook_rows(path)
    genomic = next((r for r in rows if "GRCh38" in r.get("A", "") and "NC_" in r.get("A", "")), None)
    rsids = next((r for r in rows if r.get("A", "").lower() == "rsid"), None)
    if not genomic or not rsids:
        return {}
    accession = re.search(r"NC_0*(\d+)\.", genomic["A"])
    if not accession:
        return {}
    chromosome = "X" if accession[1] == "23" else accession[1]
    markers = {}
    for column, rsid in rsids.items():
        matching = [rid for rid in requested if rid in re.split(r"[,;\s]+", rsid)]
        if not matching:
            continue
        hgvs = genomic.get(column, "")
        snv = re.fullmatch(r"g\.(\d+)([ACGT])>([ACGT])", hgvs)
        position = re.search(r"g\.(\d+)", hgvs)
        if not position:
            continue
        aliases = [r["A"] for r in rows if r.get("A", "").startswith("*") and r.get(column, "")]
        for rid in matching:
            markers[rid] = {"rsid": rid, "chrom": chromosome, "pos": int(position[1]),
                            "ref": snv[2] if snv else None, "alt": snv[3] if snv else None,
                            "definition_hgvs": hgvs, "allele_names_with_marker": aliases,
                            "source": "CPIC_current_allele_definition_workbook",
                            "complex": snv is None}
    return markers


def public_marker(cache, gene, rsid):
    data = cache.get(f"https://rest.ensembl.org/variation/human/{rsid}")
    matches = [m for m in data.get("mappings", []) if m.get("assembly_name") == "GRCh38"
               and m.get("seq_region_name") in [str(i) for i in range(1, 23)] + ["X"]]
    if not matches:
        return None
    mapping = matches[0]
    alleles = mapping.get("allele_string", "").split("/")
    return {"gene": gene, "rsid": rsid, "chrom": mapping["seq_region_name"],
            "pos": mapping["start"], "ref": alleles[0], "alt": alleles[1] if len(alleles) == 2 else None,
            "alleles": alleles, "definition_hgvs": "", "source": "Ensembl_dbsnp_GRCh38_public_marker",
            "complex": mapping["end"] != mapping["start"] or any(len(a) != 1 for a in alleles)}


def query_marker(vcf, marker):
    """Return only explicitly called, QC-passing marker genotypes; absence is unknown."""
    if marker["complex"]:
        return {"status": "complex_marker_not_called", "genotype": "unknown", "reason": "repeat/indel requires normalized validation; absence is not reference"}
    region = f"{marker['chrom']}:{marker['pos']}-{marker['pos']}"
    text = subprocess.check_output(["bcftools", "query", "-r", region, "-f",
        "%POS\t%REF\t%ALT\t%FILTER\t%QUAL[\t%GT\t%DP\t%GQ]\n", str(vcf)], text=True, stderr=subprocess.PIPE)
    for line in text.splitlines():
        pos, ref, alt, filt, qual, gt, dp, gq = line.split("\t")
        if int(pos) != marker["pos"] or ref != marker["ref"]:
            continue
        if filt != "PASS" or "." in (qual, dp, gq) or float(qual) < 30 or int(dp) < 10 or int(gq) < 30:
            return {"status": "failed_qc", "genotype": "unknown"}
        values = [ref] + alt.split(",")
        try:
            genotype = "/".join(values[int(a)] for a in re.split(r"[/|]", gt))
        except (ValueError, IndexError):
            return {"status": "missing_genotype", "genotype": "unknown"}
        return {"status": "observed", "genotype": genotype, "gt": gt, "dp": int(dp), "gq": int(gq),
                "chrom": marker["chrom"], "pos": marker["pos"], "ref": ref, "alt": alt}
    return {"status": "absent_variant_only_vcf", "genotype": "unknown",
            "reason": "no record; no gVCF or reference coverage available"}


def marker_phenotype(gene, observations):
    observed = {rid: result for rid, result in observations.items() if result["status"] == "observed"}
    if gene == "CYP2D6":
        return "Indeterminate: CNV, hybrid genes, phase and paralog mapping not resolved"
    if gene == "ABCB1":
        return "No validated CPIC metabolizer phenotype/dose assignment"
    if gene == "CYP3A5" and "rs776746" in observed:
        genotype = observed["rs776746"]["genotype"]
        if genotype == "C/C":
            return "Marker compatible with two CYP3A5*3 splice alleles/reduced expression; full diplotype unconfirmed"
        if "C" in genotype:
            return "CYP3A5*3 splice marker present; other haplotype function and full diplotype unresolved"
    if gene == "SLCO1B1" and "rs4149056" in observed:
        genotype = observed["rs4149056"]["genotype"]
        if "C" in genotype:
            return "Decreased-function c.521T>C marker present; full transporter phenotype/diplotype unresolved"
    if gene == "CYP2C19" and "rs4244285" in observed and "A" in observed["rs4244285"]["genotype"]:
        if "rs12248560" in observed and "T" in observed["rs12248560"]["genotype"]:
            return "No-function *2 splice marker AND increased-expression *17 marker present; phase/full diplotype unresolved; no rapid-metabolizer assignment"
        return "No-function *2 splice marker present; full phenotype/diplotype unresolved"
    if gene == "CYP2C19" and "rs12248560" in observed and "T" in observed["rs12248560"]["genotype"]:
        return "Increased-expression *17 promoter marker present; full metabolizer phenotype unresolved"
    if not observed:
        return "Indeterminate: no QC-passing selected marker genotype; absence is not normal function"
    return "Indeterminate full phenotype: selected marker genotype(s) only; phase/CNV/coverage incomplete"


def pgx(vcf, cache, work, output):
    all_markers, observations, failures, tables = {}, {}, {}, {}
    for gene, requested in PGX_MARKERS.items():
        markers = {}
        if gene not in ("CYP3A4", "ABCB1", "CYP2D6"):
            url = f"https://files.cpicpgx.org/data/report/current/allele_definition/{gene}_allele_definition_table.xlsx"
            path = work / f"clinical_{gene}_alleles.xlsx"
            try:
                download_public(url, path)
                markers = parse_cpic_markers(path, requested)
                tables[gene] = {"url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                "requested_markers_found": len(markers)}
            except Exception as exc:
                failures[gene] = f"CPIC_workbook_{type(exc).__name__}"
        for rsid in requested:
            if rsid not in markers:
                try:
                    marker = public_marker(cache, gene, rsid)
                    if marker:
                        markers[rsid] = marker
                except Exception as exc:
                    failures[rsid] = type(exc).__name__
        all_markers[gene] = markers
        observations[gene] = {rid: query_marker(vcf, marker) for rid, marker in markers.items()}
    (work / "clinical_pgx_markers_and_calls.json").write_text(json.dumps({"markers": all_markers, "observations": observations}, indent=2))
    rows = []
    for gene, requested in PGX_MARKERS.items():
        results = observations[gene]
        observed = [f"{rid} {r['genotype']} (DP{r['dp']}; GQ{r['gq']})" for rid, r in results.items() if r["status"] == "observed"]
        rows.append({"gene": gene, "genotype": "; ".join(observed) or "No QC-passing selected genotype; unknown",
                     "phenotype": marker_phenotype(gene, results), "candidate_drugs_affected": DRUG_RELEVANCE[gene],
                     "selected_markers": len(requested), "resolved_public_marker_coordinates": len(results),
                     "observed_markers": len(observed), "full_star_diplotype": "not_assigned",
                     "clinical_action": "Confirm with clinical PGx assay/coverage and drug-specific guidance before dosing"})
    write_csv(output / "pgx_dosing_relevance.csv", rows)
    publications_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=" + ",".join(PGX_PMIDS) + "&retmode=xml"
    publications_path = download_public(publications_url, work / "clinical_pgx_guideline_abstracts.xml")
    publication_root = ET.fromstring(publications_path.read_bytes())
    references = []
    for article in publication_root.findall("PubmedArticle"):
        references.append({"PMID": article.findtext(".//PMID"), "title": article.findtext(".//ArticleTitle"),
                           "DOI": next((item.text for item in article.findall(".//ArticleId") if item.attrib.get("IdType") == "doi"), None)})
    if {r["PMID"] for r in references} != set(PGX_PMIDS):
        raise ValueError("PGx primary-reference metadata incomplete")
    summary = {"genes": len(rows), "selected_markers": sum(map(len, PGX_MARKERS.values())),
               "observed_markers": sum(r["observed_markers"] for r in rows),
               "called_complete_diplotypes": 0, "definition_tables": tables, "reference_failures": failures,
               "primary_guidelines": references, "primary_metadata_url": publications_url,
               "primary_metadata_sha256": hashlib.sha256(publications_path.read_bytes()).hexdigest(),
               "marker_status_counts": dict(Counter(r["status"] for gene in observations.values() for r in gene.values())),
               "null_calibration": "No association/enrichment test. Missing markers are unknown; common marker controls are not pathogenic findings. Only a marker screen, no phenotype-negative claim."}
    (output / "pgx_summary.json").write_text(json.dumps(summary, indent=2))
    table = "\n".join(f"| {r['gene']} | {r['genotype']} | {r['phenotype']} | {r['candidate_drugs_affected']} |" for r in rows)
    (output / "pgx_report.md").write_text("# Key-marker pharmacogenomics\n\n" +
        f"Screened {summary['selected_markers']} prespecified markers in 12 genes; {summary['observed_markers']} explicit QC-passing marker genotypes. No full star diplotypes or normal-function calls are assigned.\n\n" +
        "| Gene | Explicitly observed genotype | Phenotype boundary | Drug relevance if considered |\n|---|---|---|---|\n" + table +
        "\n\nMissing records in this variant-only VCF are **unknown**, not homozygous reference. GRCh38 coordinates come from current CPIC allele-definition workbooks, with public Ensembl/dbSNP fallback. CPIC genomic strand is used; alleles are reported on the forward GRCh38 strand. Complex repeats/indels (including UGT1A1*28) are deliberately not called. CYP2D6 CNV/hybrids and paralog mapping are unresolved. ABCB1 variants do not create a validated CPIC dosing recommendation. Functional marker associations do not establish complete haplotypes; no absence-of-risk inference is made. These genes cannot personalize drugs whose relevant pharmacogenes are not measured.\n\n" +
        "Methods/control: PASS, QUAL>=30, DP>=10 and GQ>=30 at explicitly represented markers; reference and position must match. The behavioral test proves absent VCF records cannot become reference or normal metabolizers. No population association or enrichment is claimed, and no clinical sensitivity or specificity is estimated from this single genome.\n\n" +
        "Sources opened: CPIC [CYP3A5/tacrolimus](https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-CYP3A5/), PMID [25801146](https://pubmed.ncbi.nlm.nih.gov/25801146/), DOI [10.1002/cpt.113](https://doi.org/10.1002/cpt.113); CPIC [statins](https://cpicpgx.org/guidelines/cpic-guideline-for-statins/), PMID [35152405](https://pubmed.ncbi.nlm.nih.gov/35152405/), DOI [10.1002/cpt.2557](https://doi.org/10.1002/cpt.2557). Full workbook URLs and SHA-256 values are in pgx_summary.json. Drug relevance is conditional on the drug being selected; it is not evidence of MVA efficacy.\n")
    report_path = output / "pgx_report.md"
    report_path.write_text(report_path.read_text() + "\nAdditional primary-guideline metadata and abstracts were opened through PubMed E-utilities: thiopurines TPMT/NUDT15 PMID[41618934](https://pubmed.ncbi.nlm.nih.gov/41618934/), DOI10.1002/cpt.70209; DPYD PMID[29152729](https://pubmed.ncbi.nlm.nih.gov/29152729/), DOI10.1002/cpt.911; G6PD PMID[36049896](https://pubmed.ncbi.nlm.nih.gov/36049896/), DOI10.1002/cpt.2735; UGT1A1/atazanavir PMID[26417955](https://pubmed.ncbi.nlm.nih.gov/26417955/), DOI10.1002/cpt.269. Cached metadata hashes and publication titles are in pgx_summary.json.\n")
    return summary


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def local_clinvar_frequencies(variants, work):
    """Legacy ExAC/1000G fields are useful common-allele exclusions, not gnomAD absence."""
    clinvar = work / "clinical_clinvar.vcf.gz"
    if not clinvar.exists() or not Path(str(clinvar) + ".tbi").exists():
        return {}
    positions = work / "clinical_modifier_positions.tsv"
    positions.write_text("".join(f"{v['chrom']}\t{v['pos']}\n" for v in variants))
    text = subprocess.check_output(["bcftools", "query", "-R", str(positions), "-f",
        "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_EXAC\t%INFO/AF_TGP\n", str(clinvar)], text=True)
    frequencies = {}
    requested = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in variants}
    for line in text.splitlines():
        chrom, position, ref, alt, exac, tgp = line.split("\t")
        values = []
        for field in (exac, tgp):
            for number in field.split(","):
                try:
                    values.append(float(number))
                except ValueError:
                    pass
        key = (chrom, int(position), ref, alt)
        if values and key in requested:
            frequencies[key] = max(values)
    return frequencies


def gnomad_gene_reference(symbol, gene_id, work):
    """Fetch public gene data only; patient variant coordinates are never sent."""
    path = work / f"clinical_gnomad_{symbol}.json"
    query = ('query { gene(gene_id: "' + gene_id + '", reference_genome: GRCh38) '
             '{ gene_id variants(dataset: gnomad_r4) '
             '{ variant_id exome { ac an af } genome { ac an af } } } }')
    if path.exists():
        raw = path.read_bytes()
    else:
        response = requests.post("https://gnomad.broadinstitute.org/api", json={"query": query}, timeout=90, stream=True)
        response.raise_for_status()
        chunks, size = [], 0
        for chunk in response.iter_content(1024 * 1024):
            size += len(chunk)
            if size > 10_000_000:
                raise ValueError("Gene reference exceeded bounded 10MB response")
            chunks.append(chunk)
        raw = b"".join(chunks)
        path.write_bytes(raw)
    data = json.loads(raw)
    gene = data.get("data", {}).get("gene")
    if data.get("errors") or not gene or gene.get("gene_id") != gene_id:
        raise ValueError("gnomAD gene reference query/schema failure")
    frequencies = {}
    for variant in gene.get("variants", []):
        values = [variant[cohort]["af"] for cohort in ("exome", "genome")
                  if variant.get(cohort) and variant[cohort].get("af") is not None]
        if values:
            frequencies[variant["variant_id"]] = max(values)
    metadata = {"url": "https://gnomad.broadinstitute.org/api", "gene": symbol, "gene_id": gene_id,
                "dataset": "gnomad_r4", "reference_genome": "GRCh38", "query": query,
                "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "public_variants": len(gene["variants"]),
                "privacy": "Gene identifier only; no patient variant/genotype sent"}
    path.with_suffix(".meta.json").write_text(json.dumps(metadata, indent=2))
    return frequencies, metadata


def modifiers(vcf, cache, work, output, max_variants=12):
    symbols = sorted({gene for panel in MODIFIER_PANELS.values() for gene in panel})
    models, failures = {}, {}
    def fetch(gene):
        try:
            model = gene_model(cache, GENE_ALIASES.get(gene, gene))
            if model:
                model["symbol"] = gene
            return gene, model, None
        except Exception as exc:
            return gene, None, type(exc).__name__
    with ThreadPoolExecutor(max_workers=4) as executor:
        for gene, model, error in executor.map(fetch, symbols):
            if model:
                models[gene] = model
            else:
                failures[gene] = error or "no_canonical_coding_transcript"
    calls, counts = [], Counter()
    for gene, model in models.items():
        variants, stats = extract_gene(vcf, model)
        calls.extend(variants)
        counts.update(stats)
    # The known causal pair is an internal retrieval control, not a modifier.
    controls = [v for v in calls if v["gene"] == "BUB1B"]
    legacy_af = local_clinvar_frequencies([v for v in calls if v["gene"] != "BUB1B"], work)
    candidates = sorted([v for v in calls if v["gene"] != "BUB1B"],
                        key=lambda v: (-IMPACT[v["consequence"]], v["gene"], v["pos"]))
    annotated = []
    for variant in candidates[:max_variants]:
        v = dict(variant)
        if len(v["ref"]) == len(v["alt"]) == 1:
            hgvs = f"{v['chrom']}:g.{v['pos']}{v['ref']}>{v['alt']}"
            try:
                response = cache.get("https://rest.ensembl.org/vep/human/hgvs/" + urllib.parse.quote(hgvs, safe=""), params={"canonical": 1, "hgvs": 1})
                frequencies = []
                for colocated in response[0].get("colocated_variants", []):
                    for name, af in colocated.get("frequencies", {}).get(v["alt"], {}).items():
                        if isinstance(af, (int, float)):
                            frequencies.append(af)
                v["max_population_af"] = max(frequencies) if frequencies else None
                tx = next((t for t in response[0].get("transcript_consequences", []) if t.get("transcript_id") == v["transcript"].split(".")[0]), {})
                v["vep"] = {k: tx[k] for k in ("hgvsc", "hgvsp", "consequence_terms", "sift_prediction", "polyphen_prediction") if k in tx}
                v["clinical_significance"] = sorted({s for colocated in response[0].get("colocated_variants", []) for s in colocated.get("clin_sig", [])})
                v["rarity"] = "rare_returned_AF<=0.001" if frequencies and max(frequencies) <= .001 else "common_returned_AF>0.001" if frequencies else "unknown_no_frequency_returned"
            except Exception as exc:
                v["rarity"] = "unknown_query_failed_" + type(exc).__name__
        else:
            v["rarity"] = "unknown_indel_not_sent"
        annotated.append(v)
    annotation_by_key = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v for v in annotated}
    all_candidates = []
    for variant in candidates:
        key = (variant["chrom"], variant["pos"], variant["ref"], variant["alt"])
        v = dict(annotation_by_key.get(key, variant))
        frequencies = [af for af in (legacy_af.get(key), v.get("max_population_af")) if af is not None]
        v["max_population_af"] = max(frequencies) if frequencies else None
        v["rarity"] = "rare_returned_AF<=0.001" if frequencies and max(frequencies) <= .001 else "common_returned_AF>0.001" if frequencies else "unknown_no_population_frequency"
        v["population_source"] = "local_ClinVar_legacy_ExAC_1000G" if key in legacy_af else "VEP" if key in annotation_by_key and frequencies else "unavailable"
        all_candidates.append(v)
    gnomad_metadata, gnomad_failures = {}, {}
    unknown_genes = sorted({v["gene"] for v in all_candidates if v["max_population_af"] is None})
    for gene in unknown_genes[:5]:
        try:
            reference = cache.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{GENE_ALIASES.get(gene, gene)}", params={"expand": 1})
            frequencies, metadata = gnomad_gene_reference(gene, reference["id"], work)
            gnomad_metadata[gene] = metadata
            for variant in all_candidates:
                if variant["gene"] != gene or variant["max_population_af"] is not None:
                    continue
                variant_id = f"{variant['chrom']}-{variant['pos']}-{variant['ref']}-{variant['alt']}"
                if variant_id in frequencies:
                    af = frequencies[variant_id]
                    variant.update(max_population_af=af, population_source="gnomAD_r4_public_gene_reference_local_match",
                                   rarity="rare_returned_AF<=0.001" if af <= .001 else "common_returned_AF>0.001")
        except Exception as exc:
            gnomad_failures[gene] = type(exc).__name__
    published = sorted(all_candidates, key=lambda v: (v["rarity"].startswith("common"), -IMPACT[v["consequence"]], v["gene"], v["pos"]))[:12]
    (work / "clinical_modifier_calls.json").write_text(json.dumps({"all_protein_altering_calls": calls, "shortlist": published, "all_candidate_annotations": all_candidates}, indent=2))
    # At most twelve named candidates are publishable findings, not a genotype-scale table.
    rows = [{"gene": v["gene"], "variant": f"chr{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}",
             "protein": v["protein"], "transcript": v["transcript"], "consequence": v["consequence"],
             "genotype": "heterozygous" if v["het"] else "homozygous_alternate", "rarity": v["rarity"],
             "max_population_af": v.get("max_population_af"), "population_source": v["population_source"],
             "interpretation": "unvalidated candidate; no evidence this modifies BUB1B-MVA"} for v in published]
    write_csv(output / "modifier_shortlist.csv", rows, fields=list(rows[0]) if rows else ["gene", "variant", "interpretation"])
    summary = {"panel_genes": MODIFIER_PANELS, "genes_requested": len(symbols), "canonical_genes_scanned": len(models),
               "resolved_gene_aliases": GENE_ALIASES,
               "failures": failures, "counts": dict(counts), "protein_altering_calls": len(calls),
               "BUB1B_retrieval_control_calls": len(controls), "external_named_variant_shortlist": len(annotated),
               "shortlist_rarity_counts": dict(Counter(v["rarity"] for v in annotated)),
               "all_noncausal_coding_rarity_counts": dict(Counter(v["rarity"] for v in all_candidates)),
               "gnomad_gene_references": gnomad_metadata, "gnomad_reference_failures": gnomad_failures,
               "local_ClinVar_AF_available_noncausal_variants": len(legacy_af),
               "local_ClinVar_AF_common_gt_0_001": sum(af > .001 for af in legacy_af.values()),
               "noncausal_calls_without_local_AF": len(candidates) - len(legacy_af),
               "verified_MVA_modifiers": 0,
               "null_calibration": "No cohort/population null: no modifier enrichment or disease-modification claim. Known BUB1B pair is an internal retrieval positive control; common AF-returned alleles are negative prioritization controls."}
    (output / "modifier_summary.json").write_text(json.dumps(summary, indent=2))
    table = "\n".join(f"| {r['gene']} | {r['variant']} / {r['protein']} | {r['rarity']} |" for r in rows)
    rarity_counts = summary["all_noncausal_coding_rarity_counts"]
    conclusion = "No rare coding candidate was found where frequency is known" if not rarity_counts.get("rare_returned_AF<=0.001", 0) else "Rare coding candidates require functional and segregation follow-up"
    (output / "modifier_report.md").write_text(f"# Modifier/second-hit research screen\n\nScanned {len(models)}/{len(symbols)} prespecified checkpoint/CIN, rhabdomyosarcoma and nephrocalcinosis genes locally. Found {len(calls)} protein-altering canonical-CDS calls including {len(controls)} BUB1B retrieval controls. Among the {len(all_candidates)} remaining coding calls, {rarity_counts.get('common_returned_AF>0.001',0)} have returned population AF>0.001, {rarity_counts.get('rare_returned_AF<=0.001',0)} AF<=0.001 and {rarity_counts.get('unknown_no_population_frequency',0)} unknown frequency. Annotated a bounded shortlist of {len(annotated)} named candidates externally. **{conclusion}; no variant is established as a BUB1B-MVA modifier.**\n\n| Gene | Named finding | Rarity evidence |\n|---|---|---|\n" + table +
        "\n\nQC and translation reuse Track1. Canonical coding exons only; gene-level phase, splice regions, fusions (including PAX3/PAX7-FOXO1), deletions/CNV and tumor-specific second hits are not assessed. This VCF is not a paired tumor-normal study. Rarity is confirmed only where a population frequency is returned; unknown is not rare. External candidates are ordered by truncating/in-frame/missense impact, then gene/coordinate. Up to twelve named variants are sent to VEP. Local exact-allele ClinVar AF_EXAC/AF_TGP fields additionally exclude common alleles throughout the panel; these are legacy ExAC/1000G frequencies, not contemporary gnomAD absence evidence. For residual unknowns, at most five bounded public gnomAD gene references can be queried by gene identifier only, then matched locally: no additional patient variant is sent. Absence in the public gene response remains unknown. Unknown/rare candidates appear first in the publication shortlist; common calls are controls. The full aggregate rarity counts and public-reference hashes/queries are in modifier_summary.json; remaining calls stay private. SLC9A3R1 was resolved through its current HGNC name NHERF1 (NCBI Gene9368).\n\nNull/control: one genome cannot establish modifier enrichment or penetrance. No population enrichment inference or null p-value is reported. BUB1B is an internal retrieval control only; common variants in the shortlist are deprioritization controls, not disease-causing discoveries. Functional causality needs isogenic rescue/perturbation and segregation.\n")
    return summary


def high_confidence_clinvar(info):
    return info.get("CLNSIG") in PLP and SF_STARS.get(info.get("CLNREVSTAT"), 0) >= 2 and "CLNSIGCONF" not in info


def secondary(vcf, work, output):
    sf_path = download_public(SF_URL, work / "clinical_acmg.csv")
    sf_rows = list(csv.reader(sf_path.open()))[2:]
    genes = {r[0] for r in sf_rows if r} | SF_EXTRA
    if len(genes) != 84:
        raise ValueError(f"ACMG v3.3 reference changed or incomplete: found {len(genes)} genes")
    clinvar = download_public(CLINVAR_URL, work / "clinical_clinvar.vcf.gz")
    known, counts, clinvar_date = {}, Counter(), None
    with gzip.open(clinvar, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("##fileDate="):
                    clinvar_date = line.strip().split("=", 1)[1]
                continue
            chrom, pos, variant_id, ref, alt, _, _, raw_info = line.rstrip().split("\t")[:8]
            info = dict(field.split("=", 1) if "=" in field else (field, "") for field in raw_info.split(";"))
            matched = genes & {x.split(":")[0] for x in info.get("GENEINFO", "").split("|")}
            if not matched:
                continue
            counts["ClinVar_records_in_SF_genes"] += 1
            if not high_confidence_clinvar(info) or "," in alt:
                continue
            counts["ClinVar_high_confidence_PLP_records"] += 1
            known[(chrom, int(pos), ref, alt)] = {"genes": sorted(matched), "clinvar_id": variant_id,
                 "classification": info["CLNSIG"], "review_stars": SF_STARS[info["CLNREVSTAT"]],
                 "disease": info.get("CLNDN", ""), "hgvs": info.get("CLNHGVS", "")}
    positions = work / "clinical_secondary_positions.tsv"
    positions.write_text("".join(f"{chrom}\t{pos}\n" for chrom, pos in sorted({(k[0], k[1]) for k in known}, key=lambda x: (int(x[0]) if x[0].isdigit() else 23, x[1]))))
    command = ["bcftools", "query", "-R", str(positions), "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ]\n", str(vcf)]
    result = subprocess.check_output(command, text=True, stderr=subprocess.PIPE)
    hits = []
    for line in result.splitlines():
        chrom, pos, ref, alts, qual, filt, gt, dp, gq = line.split("\t")
        counts["patient_records_at_reference_positions"] += 1
        if filt != "PASS" or "." in (qual, dp, gq) or float(qual) < 30 or int(dp) < 10 or int(gq) < 30:
            continue
        for index, alt in enumerate(alts.split(","), 1):
            if str(index) not in re.split(r"[/|]", gt):
                continue
            match = known.get((chrom, int(pos), ref, alt))
            if match:
                hits.append({"variant": f"chr{chrom}:{pos} {ref}>{alt}", "gt": gt, "dp": int(dp), "gq": int(gq), **match})
    (work / "clinical_secondary_matches.json").write_text(json.dumps(hits, indent=2))
    rows = []
    for hit in hits:
        for gene in hit["genes"]:
            inheritance = "AR" if gene in AR_SF or gene == "HFE" else "XL" if gene in {"ABCD1", "GLA", "OTC"} else "AD_or_semidominant"
            state = "candidate_requires_gene_disease_and_clinical_validation"
            alleles = re.split(r"[/|]", hit["gt"])
            if inheritance == "AR" and len(set(alleles)) > 1:
                state = "single_heterozygous_AR_allele_not_reportable_without_second_trans_allele"
            if gene == "HFE":
                state = "requires_p_C282Y_homozygote_specific_rule"
            if gene == "ABCD1" and len(alleles) == 2 and len(set(alleles)) > 1:
                state = "single_heterozygous_ABCD1_not_returnable_under_v3_3_rule"
            rows.append({"gene": gene, "variant": hit["variant"], "genotype": hit["gt"], "classification": hit["classification"],
                         "review_stars": hit["review_stars"], "clinvar_id": hit["clinvar_id"], "condition": hit["disease"],
                         "finding_type": "secondary_candidate", "reportability": state, "inheritance": inheritance})
    write_csv(output / "secondary_candidates.csv", rows, fields=list(rows[0]) if rows else ["gene", "variant", "finding_type", "reportability"])
    summary = {"panel": "ACMG_SF_v3.3", "panel_gene_count": len(genes), "panel_source": SF_URL,
               "panel_reconciled_from_paper": sorted(SF_EXTRA), "paper_PMID": "40568962",
               "public_ClinVar_fileDate": clinvar_date, "public_ClinVar_url": CLINVAR_URL,
               "public_ClinVar_bytes": clinvar.stat().st_size, "public_ClinVar_sha256": hashlib.sha256(clinvar.read_bytes()).hexdigest(),
               "counts": dict(counts), "exact_high_confidence_PLP_variant_matches": len(hits),
               "clinically_confirmed_secondary_findings": 0,
               "null_calibration": "Deterministic reference-match screen, not statistical burden inference. >=2-star P/LP and exact allele/QC safeguards reduce false candidate reporting; absent matches do not establish a negative clinical screen."}
    (output / "secondary_summary.json").write_text(json.dumps(summary, indent=2))
    (output / "secondary_report.md").write_text(f"# Secondary findings research screen\n\nScreened all {len(genes)} ACMG SF v3.3 genes against the public GRCh38 ClinVar file. Found **{len(hits)} exact, QC-passing matches** to >=2-star pathogenic/likely pathogenic alleles without aggregate conflicts. **0 clinically confirmed secondary findings**; clinical confirmation is outside this research pipeline.\n\n" +
        "The current ClinGen curation CSV contains 81 distinct genes; SMAD3, TSC1 and TSC2 were reconciled from the official v3.3 paper's Table1 to recover the full84. The public ClinVar reference stays under ignored work/. Patient calls are matched locally by GRCh38 chromosome, position, reference and alternate allele, with PASS, QUAL>=30, DP>=10 and GQ>=30. Multiallelic called alleles are evaluated separately. Two-star multiple-submitters/no-conflict, expert panel and practice-guideline classifications qualify; conflicts and weaker evidence are excluded.\n\n" +
        "All matches are secondary candidates only. AR single heterozygotes do not establish a reportable finding; phase and second-allele evidence are needed. HFE requires the specific C282Y homozygote rule; ABCD1, GLA and OTC follow distinct X-linked return rules. RYR1/CACNA1S require the malignant-hyperthermia gene-disease relationship; TTN requires truncating variants in relevant exons. Clinical gene-disease/variant-type review is required before any return.\n\n" +
        "Limitations: exact representation matching can miss differently normalized indels/complex alleles; no CNV/SV, novel loss-of-function classification, splice rescue or full coverage assessment was performed. The screen can miss clinically reportable variants absent from ClinVar or carrying less review evidence. A zero-match result is not a negative ACMG clinical analysis. This deterministic screen has no cohort null or validated clinical sensitivity; no enrichment inference is made. Behavioral tests verify rejection of one-star/conflicting/benign variants.\n\n" +
        "Source opened: Lee et al., ACMG SF v3.3, PMID [40568962](https://pubmed.ncbi.nlm.nih.gov/40568962/), DOI [10.1016/j.gim.2025.101454](https://doi.org/10.1016/j.gim.2025.101454); [ClinGen resource](https://search.clinicalgenome.org/kb/genes/acmgsf). ClinVar sourceURL, size and SHA-256 are recorded in secondary_summary.json.\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--step", choices=("modifiers", "pgx", "secondary", "all"), default="all")
    parser.add_argument("--work", type=Path, default=Path("work"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    if not args.work.resolve().is_relative_to((Path.cwd() / "work").resolve()):
        parser.error("Patient intermediates must remain beneath ignored work/")
    args.work.mkdir(parents=True, exist_ok=True)
    if subprocess.run(["git", "check-ignore", "-q", str(args.work / "clinical_secondary_matches.json")]).returncode:
        parser.error("Patient intermediate path must be git-ignored")
    args.output.mkdir(parents=True, exist_ok=True)
    cache = ReferenceCache(args.work / "clinical-cache")
    summary = {}
    if args.step in ("modifiers", "all"):
        summary["modifiers"] = modifiers(args.vcf, cache, args.work, args.output)
    if args.step in ("pgx", "all"):
        summary["pgx"] = pgx(args.vcf, cache, args.work, args.output)
    if args.step in ("secondary", "all"):
        summary["secondary"] = secondary(args.vcf, args.work, args.output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
