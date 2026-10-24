"""A conservative public graph engine: association is not a rescue mechanism.

ChEMBL mechanism records generate candidates. Reactome membership supplies a
bounded graph link, never a causal direction. Phase-4 indication labels and explicit
literature evidence can promote candidates; otherwise the engine abstains from
calling graph proximity a treatment strategy.
"""

from collections import defaultdict
import csv
import json
from pathlib import Path
import re

from .cache import PublicCache
from .curation import curated_hypotheses


OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
# A curated regulator check can override an outdated phase-4 database proxy.
# This says nothing about other jurisdictions; they require separate verification.
REGULATORY_HOLDS = {
    "ATALUREN": {
        "reason": "EU conditional marketing authorisation was not renewed on 2025-03-28; another current eligible jurisdiction has not been verified",
        "source": "https://www.ema.europa.eu/en/medicines/human/EPAR/translarna",
        "effective_date": "2025-03-28", "status": "hold_pending_current_jurisdiction_verification",
    },
}
TARGET_QUERY = """query Target($id: String!) {
  target(ensemblId: $id) {
    id approvedSymbol proteinIds { id source }
    pathways { pathwayId pathway }
    drugAndClinicalCandidates { count rows { id maxClinicalStage } }
  }
}"""


def chunks(items, size=40):
    items = list(items)
    for index in range(0, len(items), size):
        yield items[index:index + size]


def words(text: str) -> set[str]:
    return set(re.findall(r"[a-z0-9]+", text.lower()))


def resolve_gene(cache: PublicCache, gene: str | None, variants: list[str]) -> str:
    """Resolve only a small named HGVS list; the API never sees a private genome."""
    if len(variants) > 10:
        raise ValueError("Supply at most 10 named variants, including when --gene is provided")
    for variant in variants:
        if len(variant) > 200 or not re.fullmatch(r"[A-Za-z0-9_.-]+:[cgmnpr]\.\S+", variant):
            raise ValueError("Each variant must be one named HGVS expression of at most 200 characters")
    if gene:
        if not re.fullmatch(r"[A-Za-z0-9-]{2,30}", gene):
            raise ValueError("Provide one HGNC gene symbol")
        return gene.upper()
    if not variants or len(variants) > 10:
        raise ValueError("Without --gene, supply 1–10 named HGVS variants")
    symbols = set()
    for variant in variants:
        from urllib.parse import quote
        payload = cache.json(
            "https://rest.ensembl.org/vep/human/hgvs/" + quote(variant, safe=""),
            params={"content-type": "application/json"},
        )
        symbols.update(
            t["gene_symbol"] for v in payload for t in v.get("transcript_consequences", [])
            if t.get("gene_symbol") and t.get("canonical") == 1
        )
    if len(symbols) != 1:
        raise ValueError(f"HGVS did not resolve to one canonical gene: {sorted(symbols)}; supply --gene")
    return symbols.pop()


def hpo_labels(cache: PublicCache, hpo: list[str]) -> list[dict]:
    # Fetch one complete public ontology. Patient-derived phenotype identifiers
    # are looked up locally and are never included in an external request.
    for term in hpo:
        if not re.fullmatch(r"HP:\d{7}", term):
            raise ValueError(f"Invalid HPO identifier: {term}")
    ontology = cache.json("https://purl.obolibrary.org/obo/hp.json") if hpo else {}
    labels = {
        node["id"].rsplit("/", 1)[-1].replace("_", ":"): node.get("lbl", "unresolved")
        for graph in ontology.get("graphs", []) for node in graph.get("nodes", [])
    }
    output = []
    for term in sorted(set(hpo)):
        output.append({"id": term, "label": labels.get(term, "unresolved")})
    return output


def graph(cache: PublicCache, gene: str, max_pathways: int, max_proteins: int):
    lookup = cache.json(
        f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}",
        params={"content-type": "application/json"},
    )
    target = cache.json(OT_URL, body={"query": TARGET_QUERY, "variables": {"id": lookup["id"]}})["data"]["target"]
    if target is None:
        raise ValueError(f"Open Targets has no target entry for {gene}")
    direct = {p["id"] for p in target["proteinIds"] if p["source"] == "uniprot_swissprot"}
    proteins: dict[str, dict] = {
        p: {"gene": gene, "pathways": [], "direct": True} for p in direct
    }
    # All leaf pathways reported for this gene are inspected in stable-ID order.
    # The documented cap bounds runtime and prevents genome-scale API use.
    pathways = sorted(target["pathways"], key=lambda p: p["pathwayId"])[:max_pathways]
    for pathway in pathways:
        refs = cache.json(
            f"https://reactome.org/ContentService/data/participants/{pathway['pathwayId']}/referenceEntities"
        )
        for ref in refs:
            if ref.get("databaseName") != "UniProt" or ref.get("schemaClass") != "ReferenceGeneProduct":
                continue
            accession = ref["identifier"].split("-")[0]
            entry = proteins.setdefault(accession, {
                "gene": (ref.get("geneName") or ref.get("name") or [accession])[0],
                "pathways": [], "direct": accession in direct,
            })
            entry["pathways"].append(pathway["pathwayId"])
    # Retain genes in the most shared pathways first; this is discovery priority,
    # not a rescue score. No biological sign is inferred from graph adjacency.
    available_proteins = len(proteins)
    retained = sorted(proteins, key=lambda p: (-int(proteins[p]["direct"]), -len(proteins[p]["pathways"]), p))[:max_proteins]
    proteins = {p: proteins[p] for p in retained}
    target_records = []
    for group in chunks(proteins):
        target_records += cache.pages("target", "targets", {"target_components__accession__in": ",".join(group)})
    # Human single-protein/complex records only. Multi-species or nonhuman
    # records cannot silently become human pharmacology.
    target_records = [t for t in target_records if t.get("organism") == "Homo sapiens"]
    # ChEMBL can annotate a drug to a binding cofactor rather than the inhibited
    # enzyme (rapalogs are attached to FKBP1A, for example). Follow one explicit
    # protein-complex membership edge so this data-model choice does not erase
    # such drugs. Keep the extra edge visible; it does not imply rescue.
    cofactors: dict[str, dict] = {}
    for record in target_records:
        if record.get("target_type") != "PROTEIN COMPLEX":
            continue
        components = [c for c in record.get("target_components", []) if c.get("accession")]
        linked = [c["accession"] for c in components if c["accession"] in proteins]
        for component in components:
            accession = component["accession"]
            if accession not in proteins and linked:
                cofactors[accession] = {
                    "gene": component.get("component_description") or accession,
                    "pathways": sorted({p for link in linked for p in proteins[link]["pathways"]}),
                    "direct": False, "via_complex": record["target_chembl_id"],
                }
    for accession in sorted(cofactors)[:50]:
        proteins[accession] = cofactors[accession]
    for group in chunks(sorted(cofactors)[:50]):
        target_records += [t for t in cache.pages("target", "targets", {
            "target_components__accession__in": ",".join(group)
        }) if t.get("organism") == "Homo sapiens"]
    target_map = {t["target_chembl_id"]: t for t in target_records}
    mechanisms = []
    for group in chunks(target_map):
        mechanisms += cache.pages("mechanism", "mechanisms", {"target_chembl_id__in": ",".join(group)})
    counts = {
        "available_leaf_pathways": len(target["pathways"]), "selected_leaf_pathways": len(pathways),
        "available_pathway_proteins": available_proteins, "retained_pathway_proteins": len(retained),
        "available_complex_cofactors": len(cofactors), "retained_complex_cofactors": min(50, len(cofactors)),
        "pathway_truncated": len(target["pathways"]) > max_pathways,
        "protein_truncated": available_proteins > max_proteins or len(cofactors) > 50,
    }
    return target, pathways, proteins, target_map, mechanisms, counts


def classify_direction(mechanism: str) -> str:
    value = mechanism.lower()
    # Partial agonists can activate or antagonize according to receptor tone;
    # inverse agonists reduce constitutive activity. Neither is a full agonist.
    if "partial agonist" in value:
        return "unknown"
    if "inverse agonist" in value:
        return "inhibit"
    if any(w in value for w in ("inhibitor", "antagonist", "negative modulator", "blocker", "degrader")):
        return "inhibit"
    if any(w in value for w in ("activator", "agonist", "positive modulator")):
        return "activate"
    return "unknown"


def load_pgx(path: Path | None) -> list[dict]:
    if path is None:
        return []
    with path.open() as handle:
        if path.suffix == ".json":
            value = json.load(handle)
            return value if isinstance(value, list) else value.get("rows", [])
        return list(csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))


def normalized_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def molecule_aliases(molecule: dict, molecules: list[dict]) -> set[str]:
    """Use explicit ChEMBL parent/synonym links; do not guess active moieties."""
    hierarchy = molecule.get("molecule_hierarchy") or {}
    parent = hierarchy.get("parent_chembl_id") or molecule["molecule_chembl_id"]
    related = [molecule] + [m for m in molecules if
        ((m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m["molecule_chembl_id"]) == parent
    ]
    aliases = {normalized_name(m.get("pref_name") or "") for m in related}
    aliases.update(normalized_name(s.get("molecule_synonym") or s.get("synonyms") or "")
                   for m in related for s in m.get("molecule_synonyms", []))
    return {a for a in aliases if len(a) >= 4 and not a.isdigit()}


def pgx_drug_text(row: dict) -> str:
    drugs = row.get("candidate_drugs_affected", row.get("candidate_drugs", row.get("drugs", row.get("which_candidate_drugs_it_affects", ""))))
    return "; ".join(drugs) if isinstance(drugs, list) else str(drugs)


def pgx_for(name: str, pgx: list[dict], aliases: set[str] | None = None) -> list[dict]:
    """Match complete ChEMBL names in prose, slash lists and salt/parent aliases."""
    out = []
    for index, row in enumerate(pgx):
        text = " " + normalized_name(pgx_drug_text(row)) + " "
        matched = sorted(a for a in (aliases or {normalized_name(name)}) if f" {a} " in text)
        if matched:
            out.append({
                "input_row": index, "gene": row.get("gene", ""), "phenotype": row.get("phenotype", ""),
                "recommendation": row.get("recommendation", row.get("clinical_action", "")),
                "status": row.get("status", ""), "matched_aliases": matched,
                "candidate_drugs_affected": pgx_drug_text(row),
                "match_caveat": "Name association only; preserve the input's CPIC, indication, phase and assay limitations",
            })
    return out


def pgx_penalty(rows: list[dict]) -> int:
    """Flag review needs without penalizing a normal-metabolizer result."""
    for row in rows:
        recommendation = str(row.get("recommendation", "")).casefold()
        status = str(row.get("status", "")).casefold()
        if status == "actionable" or re.search(r"\b(avoid|reduce|reduced|adjust|alternative)\b", recommendation):
            return -10
    for row in rows:
        if re.search(r"\b(poor|intermediate|ultrarapid)\b", str(row.get("phenotype", "")).casefold()):
            return -3
    return 0


def published_evidence(cache: PublicCache, evidence: dict | None, gene: str) -> dict[str, list[dict]]:
    """Supplied bridges are explicit inputs; PubMed existence is verified.

    Retrieval verifies citation identity, not the curator's causal interpretation.
    These assertions remain marked as supplied and must be inspected by a reviewer.
    """
    out: dict[str, list[dict]] = defaultdict(list)
    if not evidence:
        return out
    for entry in evidence.get("evidence", []):
        if entry.get("gene") != gene:
            continue
        if entry.get("direction") not in {"rescue", "against", "unknown", "hypothesis"}:
            raise ValueError("Evidence bridge direction must be rescue, against, unknown, or hypothesis")
        score = entry.get("evidence_score", 0)
        if not isinstance(score, int) or not 0 <= score <= 25:
            raise ValueError("Supplied evidence_score must be an integer from 0 to 25")
        pmids = entry.get("pmids", [])
        if not pmids or not all(re.fullmatch(r"\d{1,10}", str(p)) for p in pmids):
            raise ValueError("Each supplied evidence bridge requires valid PMID identifiers")
        opened = []
        for pmid in pmids:
            response = cache.json("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params={
                "query": f"EXT_ID:{pmid} AND SRC:MED", "format": "json", "resultType": "lite",
            })
            matches = [r for r in response.get("resultList", {}).get("result", []) if r.get("id") == str(pmid)]
            if not matches:
                raise ValueError(f"Supplied PMID not found in Europe PMC: {pmid}")
            opened.append({"pmid": str(pmid), "title": matches[0].get("title"), "doi": matches[0].get("doi")})
        item = {**entry, "opened_citations": opened, "origin": "supplied_literature_bridge"}
        out[entry["drug"].casefold()].append(item)
    return out


def run(cache: PublicCache, *, gene: str | None, variants: list[str], hpo: list[str],
        mechanism: str = "unknown", disease: str | None = None,
        pgx: list[dict] | None = None, evidence: dict | None = None,
        approval_audit: dict | None = None,
        max_pathways: int = 20, max_proteins: int = 200) -> dict:
    if len(variants) > 10:
        raise ValueError("Only small named variant lists (at most 10) are accepted")
    if mechanism not in {"loss_of_function", "gain_of_function", "unknown"}:
        raise ValueError("Unknown disease mechanism")
    gene = resolve_gene(cache, gene, variants)
    expected_variants = (evidence or {}).get("curation", {}).get("applicable_variants", [])
    if expected_variants and sorted(variants) != sorted(expected_variants):
        raise ValueError("The curated variant-specific assessment does not match the supplied named variants")
    phenotype = hpo_labels(cache, hpo)
    target, pathways, proteins, target_map, mechanisms, counts = graph(cache, gene, max_pathways, max_proteins)
    by_drug: dict[str, list[dict]] = defaultdict(list)
    for item in mechanisms:
        by_drug[item["molecule_chembl_id"]].append(item)
    bridges = published_evidence(cache, evidence, gene)
    # Literature may propose a drug outside the one-hop pathway neighborhood.
    # Its identity and approval still have to pass the same public API checks.
    for name in bridges:
        molecules = cache.pages("molecule", "molecules", {"pref_name__iexact": name})
        for molecule in molecules:
            by_drug.setdefault(molecule["molecule_chembl_id"], [])
    molecules = []
    for group in chunks(by_drug):
        molecules += cache.pages("molecule", "molecules", {"molecule_chembl_id__in": ",".join(group)})
    # Retrieve explicit parent molecules for active-moiety PGx matching. A salt
    # record need not list the unsalted name among its synonyms. These parent
    # records supply names only; they do not generate new drug hypotheses.
    known_ids = {m["molecule_chembl_id"] for m in molecules}
    parent_ids = sorted({(m.get("molecule_hierarchy") or {}).get("parent_chembl_id")
                         for m in molecules if float(m.get("max_phase") or 0) == 4} - known_ids - {None})
    alias_molecules = list(molecules)
    for group in chunks(parent_ids):
        alias_molecules += cache.pages("molecule", "molecules", {"molecule_chembl_id__in": ",".join(group)})
    candidates, excluded = [], []
    for molecule in molecules:
        name = molecule.get("pref_name") or molecule["molecule_chembl_id"]
        # max_phase=4 is the hackathon's requested market-approval proxy. It
        # does not establish current availability in a particular jurisdiction.
        if float(molecule.get("max_phase") or 0) != 4 or molecule.get("withdrawn_flag"):
            excluded.append({"drug": name, "reason": "not max_phase=4 or flagged withdrawn"})
            continue
        hold = REGULATORY_HOLDS.get(name.upper())
        if hold:
            excluded.append({"drug": name, "reason": hold["reason"], "regulatory_check": hold})
            continue
        mid = molecule["molecule_chembl_id"]
        edges = []
        direct_conflict = False
        direct_aligned = False
        for item in by_drug[mid]:
            record = target_map.get(item["target_chembl_id"], {})
            components = [c["accession"] for c in record.get("target_components", []) if c.get("accession")]
            linked = [p for p in components if p in proteins]
            direct = any(proteins[p]["direct"] for p in linked)
            direction = classify_direction(item.get("mechanism_of_action") or "")
            conflict = direct and ((mechanism == "loss_of_function" and direction == "inhibit") or
                                  (mechanism == "gain_of_function" and direction == "activate"))
            aligned = direct and ((mechanism == "loss_of_function" and direction == "activate") or
                                 (mechanism == "gain_of_function" and direction == "inhibit"))
            direct_conflict |= conflict
            direct_aligned |= aligned
            edges.append({
                "target_chembl_id": item["target_chembl_id"], "target": record.get("pref_name"),
                "proteins": linked, "mechanism": item.get("mechanism_of_action"),
                "direction": direction, "direct": direct,
                "shared_pathways": sorted({r for p in linked for r in proteins[p]["pathways"]}),
                "causal_direction_supported": aligned,
            })
        indications = cache.pages("drug_indication", "drug_indications", {"molecule_chembl_id": mid})
        phase4_indications = [i for i in indications if float(i.get("max_phase_for_ind") or 0) == 4]
        # Requiring an exact disease label avoids conflating cancer treatment
        # with rescue of a cancer-predisposition syndrome.
        matched = [i for i in phase4_indications if disease and disease.casefold() in {
            (i.get("efo_term") or "").casefold(), (i.get("mesh_heading") or "").casefold(),
        }]
        supplied = bridges.get(name.casefold(), [])
        # A source identity or an indirect causal bridge cannot establish rescue.
        # Legacy "rescue" assertions need an explicit functional observation and
        # a complete path before they qualify as supported literature evidence.
        positive = [e for e in supplied if e.get("direction") == "rescue"
                    and e.get("functional_rescue_demonstrated") is True and not e.get("missing_bridges")]
        hypotheses = [e for e in supplied if e.get("direction") in {"hypothesis", "rescue"} and e not in positive]
        negative = [e for e in supplied if e.get("direction") == "against"]
        pgx_rows = pgx_for(name, pgx or [], molecule_aliases(molecule, alias_molecules))
        # Small integer components are audit priorities, not clinical efficacy
        # probabilities. Unknown pediatric suitability receives no credit.
        components = {
            "approval": 10,
            "mechanism_link": 20 if direct_aligned else (5 if edges else 0),
            "evidence_level": 30 if matched else max((int(e.get("evidence_score", 0)) for e in positive), default=0),
            "pediatric_use": max((5 if e.get("pediatric_use") == "established" else 0 for e in positive), default=0),
            "pgx": pgx_penalty(pgx_rows),
            "direction_conflict": -100 if direct_conflict else 0,
            "negative_evidence": -min(50, 15 * len(negative)),
        }
        score = sum(components.values())
        if direct_conflict:
            status, grade = "exclude_direction_conflict", "X"
        elif matched:
            status, grade = "known_disease_label_match", "L"
        elif positive:
            status, grade = "supplied_mechanistic_hypothesis", "C"
        elif hypotheses:
            status, grade = "curated_unverified_laboratory_hypothesis", "H"
        elif direct_aligned:
            status, grade = "direction_aligned_target_hypothesis", "D"
        else:
            status, grade = "pathway_context_only", "E"
        candidates.append({
            "drug": name, "chembl_id": mid, "score": score, "components": components,
            "grade": grade, "status": status, "mechanism_edges": edges,
            "approval": {"max_phase": molecule["max_phase"], "withdrawn_flag": molecule.get("withdrawn_flag"),
                         "source": f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mid}.json",
                         "meaning": "ChEMBL phase-4 proxy; current jurisdiction-specific approval not verified",
                         "current_approval_verified": False},
            "phase4_indication_labels": sorted({i.get("efo_term") or i.get("mesh_heading") for i in phase4_indications
                                            if i.get("efo_term") or i.get("mesh_heading")}),
            "disease_matches": matched, "supplied_evidence": supplied, "pgx": pgx_rows,
            "indication_scope": "Disease label only; route, formulation, disease manifestation and age applicability are unknown",
            "pediatric_suitability": "supplied evidence" if components["pediatric_use"] else "unknown; no pediatric safety credit",
            "caveat": "Hypothesis for experimental review. No efficacy, dose, or patient suitability is inferred from a shared pathway.",
        })
    candidates.sort(key=lambda r: (-r["score"], r["drug"], r["chembl_id"]))
    supported = [r for r in candidates if r["grade"] in {"C", "D"}]
    pgx_matches = {row["input_row"] for candidate in candidates for row in candidate["pgx"]}
    unmatched_pgx = [{
        "input_row": index, "gene": row.get("gene", ""),
        "candidate_drugs_affected": pgx_drug_text(row),
        "reason": "No explicit drug name/ChEMBL alias matched this bounded candidate list; drug-class labels are not expanded and this is not a normal PGx result",
    } for index, row in enumerate(pgx or []) if index not in pgx_matches]
    return {
        "schema_version": 2, "gene": gene, "variants": variants, "hpo": phenotype,
        "variant_role": ("Named-variant identity checked against the curated bundle; effects and phase are not inferred from the strings"
                         if expected_variants else "Provenance only when --gene is supplied; variant effects do not drive the score"),
        "hpo_role": "Local label/context lookup only; HPO terms do not drive drug scores",
        "mechanism": mechanism, "disease": disease,
        "graph": {"ensembl_id": target["id"], "pathways": pathways, "proteins": proteins,
                  "max_pathways": max_pathways, "max_proteins": max_proteins, "discovery_counts": counts},
        "ranking": candidates, "supported_hypotheses": supported,
        "known_disease_label_matches": [r["drug"] for r in candidates if r["grade"] == "L"],
        "abstention": None if supported else "No direction-supported rescue hypothesis with a ChEMBL phase-4 proxy was established. Graph proximity or a broad disease-label match alone is insufficient.",
        "excluded_not_eligible": excluded,
        "pgx_integration": {"input_rows": len(pgx or []), "matched_input_rows": len(pgx_matches), "unmatched_rows": unmatched_pgx},
        "curated_hypotheses": curated_hypotheses(bridges, evidence, approval_audit, candidates),
        "limitations": [
            "No patient genome is required or uploaded; input is a small named-variant list and HPO terms.",
            "Phase-4 ChEMBL status is a market-approval proxy; jurisdiction and current availability require regulator verification.",
            "Exact disease-label matches are retrieval benchmarks, not novel predictions; route, formulation, manifestation and age may differ (for example topical sirolimus for TSC facial angiofibromas).",
            "With --gene supplied, named variants are provenance only; HPO labels/context do not drive scores.",
            "Shared Reactome membership does not establish causal direction, tissue efficacy, or genotype rescue.",
            "HPO terms are resolved for provenance, but do not assign unvalidated phenotype-to-drug scores.",
            "Missing PGx matches are unknown, never proof of normal metabolism; any matched PGx needs manual review.",
            "Pediatric use is unknown unless a supplied evidence bridge explicitly documents it.",
            "Literature bridge citations are retrieved; causal assertions remain the curator's responsibility.",
            "Evidence scores are transparent review priorities, not calibrated treatment probabilities.",
        ],
        "api_provenance": list(cache.used.values()),
    }


def write_outputs(result: dict, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    with prefix.with_suffix(".csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", "drug", "chembl_id", "score", "grade", "status", "pgx_review"])
        writer.writeheader()
        for rank, item in enumerate(result["ranking"], 1):
            writer.writerow({"rank": rank, **{k: item[k] for k in ("drug", "chembl_id", "score", "grade", "status")}, "pgx_review": bool(item["pgx"])})
    lines = [f"# Public-data repurposing engine: {result['gene']}", "", result["abstention"] or "Direction-supported hypotheses are listed below; this is research triage, not a treatment recommendation.", ""]
    curated = result.get("curated_hypotheses", {}).get("ranked", [])
    if curated:
        with prefix.with_suffix(".curated.csv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=["drug", "chembl_id", "external_judgment_priority", "evidence_cap", "grade", "approval_gate", "laboratory_eligible", "rejection_reason", "clinical_eligible", "positive_computational_rescue_lines"])
            writer.writeheader()
            for row in curated:
                writer.writerow({"drug": row["display_name"], "chembl_id": row["chembl_id"],
                                 "external_judgment_priority": row["assessment"]["research_priority"],
                                 "evidence_cap": row["assessment"]["evidence_cap"],
                                 "grade": row["grade"], "laboratory_eligible": row["laboratory_eligible"],
                                 "rejection_reason": row["rejection_reason"],
                                 "approval_gate": row["approval"]["gate"], "clinical_eligible": row["clinical_eligible"],
                                 "positive_computational_rescue_lines": row["independent_positive_computational_rescue_lines"]})
        lines += ["## Curated laboratory hypotheses", "", "These priorities are computed from explicit supplied ordinal judgments. Missing causal bridges remain hypotheses; zero independent positive computational rescue lines and no clinical eligibility are inferred. Direction conflicts remain hard rejections regardless of priority. Clinical gate assessments never authorize treatment.", "",
                  "| Drug | Judgment-based priority | E cap | Product-approval gate | Laboratory eligible | Clinical eligible |", "|---|---:|---:|---|---|---|"]
        for row in curated:
            a = row["assessment"]
            lines.append(f"| {row['display_name']} | {a['research_priority']} | {a['evidence_cap']} | {row['approval']['gate']} | {row['laboratory_eligible']} | {row['clinical_eligible']} |")
        lines += ["", result["curated_hypotheses"]["selection_warning"], ""]
        for row in curated:
            lines += [f"### {row['display_name']}", "", "Missing bridges: " + "; ".join(row["missing_bridges"]), "", "Next experiment: " + str(row["next_experiment"]), ""]
            if row["rejection_reason"]:
                lines += ["Hard rejection: " + row["rejection_reason"], ""]
        lines += ["## Automatic graph context (separate from curated priority)", ""]
    lines += [
             "| Rank | Drug with ChEMBL phase-4 proxy | Score | Grade | Evidence status |", "|---|---|---:|---|---|"]
    for rank, item in enumerate(result["ranking"][:10], 1):
        lines.append(f"| {rank} | {item['drug']} | {item['score']} | {item['grade']} | {item['status']} |")
    lines += ["", "Grades: H = unverified curated laboratory hypothesis; L = known disease-label match, with route/manifestation/age unverified; C = supplied complete functional-rescue bridge; D = direct target with compatible direction; E = graph context only; X = incompatible direct action.", "", "Scores are review priorities. They are not efficacy probabilities. Grade E drugs are not rescue candidates. Grade L is a retrieval benchmark, not a novel prediction or automatic rescue recommendation. Equal scores are sorted alphabetically; their displayed order is not biological evidence.", "", "## PGx integration", "", f"Input rows: {result['pgx_integration']['input_rows']}; matched rows: {result['pgx_integration']['matched_input_rows']}. Unmatched rows are explicit in JSON and do not establish normal metabolism.", "", "## Limitations", ""]
    lines += [f"- {item}" for item in result["limitations"]]
    prefix.with_suffix(".md").write_text("\n".join(lines) + "\n")
