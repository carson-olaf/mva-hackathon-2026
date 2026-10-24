"""Audit explicit external judgments without upgrading incomplete causal paths."""

from datetime import datetime, timezone
import re


GATES = (
    "phase_and_genotype", "exact_genotype_functional_rescue", "feasible_unbound_exposure",
    "pediatric_pk_formulation", "renal_hepatic_interactions", "cancer_safety",
    "chronic_safety", "monitoring", "institutional_clinical_authorization",
)


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def score_assessment(entry: dict, weights: dict) -> dict:
    """Compute supplied ordinal priorities; efficacy caps dominate arithmetic."""
    components = entry["assessed_components"]
    if set(components) != {"G", "E", "S", "F", "N"} or set(weights) != set(components):
        raise ValueError("Curated scores/weights must have exactly G, E, S, F and N")
    if any(not isinstance(v, int) or not 0 <= v <= 5 for v in components.values()):
        raise ValueError("Curated components must be integers from 0 to 5")
    if any(not isinstance(v, int) or v <= 0 for v in weights.values()) or sum(weights.values()) != 20:
        raise ValueError("Curated positive integer weights must sum to 20")
    applicability = entry.get("evidence_applicability", "computation_only")
    caps = {"indirect_or_analogous_protein": 1, "bubr1_abundance_without_function": 2,
            "computation_only": 0, "exact_genotype_functional_rescue": 3}
    if applicability not in caps:
        raise ValueError("Unknown curated evidence applicability")
    cap = min(caps[applicability], 1) if entry.get("missing_bridges") else caps[applicability]
    effective = {**components, "E": min(components["E"], cap)}
    total = sum(effective[k] * weights[k] for k in weights)
    scenarios = []
    for component in weights:
        for change in (-1, 1):
            altered = {**effective, component: min(5, max(0, effective[component] + change))}
            altered["E"] = min(altered["E"], cap)
            scenarios.append({"component": component, "change": change,
                              "score": sum(altered[k] * weights[k] for k in weights)})
    return {"requested_components": components, "effective_components": effective,
            "weights": weights, "evidence_cap": cap, "research_priority": total,
            "one_point_sensitivity": scenarios,
            "one_point_score_range": [min(x["score"] for x in scenarios), max(x["score"] for x in scenarios)],
            "meaning": "Computed from externally supplied ordinal judgments; not learned discovery, efficacy probability or a clinical score"}


def approval_for(entry: dict, audit: dict | None) -> dict:
    if not audit:
        return {"gate": "UNKNOWN", "reason": "No independent official-product approval audit supplied"}
    names = {normalized(n) for n in [entry["drug"], *entry.get("drug_aliases", [])]}
    matches = [r for r in audit.get("drugs", []) if normalized(r.get("drug", "")) in names]
    if not matches:
        return {"gate": "UNKNOWN", "reason": "Drug identity did not match the independent approval audit"}
    record = matches[0]
    official = [s for s in record.get("sources", []) if re.match(
        r"https://(?:api\.fda\.gov|www\.accessdata\.fda\.gov|www\.fda\.gov|www\.ema\.europa\.eu)/", s.get("url", ""))
        and re.fullmatch(r"[0-9a-f]{64}", s.get("sha256", ""))]
    verified = record.get("approval_status") == "verified_current_prescription_product" and bool(official)
    checked = audit.get("verified_at_utc", "")
    try:
        reviewed = datetime.fromisoformat(checked.replace("Z", "+00:00").replace(" UTC", "+00:00"))
        if reviewed.tzinfo is None:
            reviewed = reviewed.replace(tzinfo=timezone.utc)
        fresh = 0 <= (datetime.now(timezone.utc) - reviewed).total_seconds() <= 7 * 86400
    except ValueError:
        fresh = False
    return {"gate": "PASS" if verified and fresh else "UNKNOWN", "product_record": record,
            "verified_at_utc": checked, "fresh_within_7_days": fresh,
            "scope": "Approved product and original indication only; not MVA efficacy, off-label access, monotherapy suitability or pediatric MVA safety"}


def curated_hypotheses(bridges: dict, bundle: dict | None, approval_audit: dict | None,
                       candidates: list[dict]) -> dict:
    if not bundle or not bundle.get("curation"):
        return {"ranked": [], "market_approved_laboratory_hypotheses": [], "unresolved_approval": []}
    by_name = {normalized(r["drug"]): r for r in candidates}
    metadata = bundle["curation"]
    weights = metadata["weights"]
    origin = metadata.get("assessment_origin", "Explicit supplied ordinal judgments; independent source review not established")
    ranked = []
    for group in bridges.values():
        for entry in group:
            if "assessed_components" not in entry:
                continue
            candidate = by_name.get(normalized(entry["drug"]))
            assessment = score_assessment(entry, weights)
            approval = approval_for(entry, approval_audit)
            gates = {gate: "UNKNOWN" for gate in GATES}
            gates.update(entry.get("clinical_gates", {}))
            if set(gates) != set(GATES) or any(v not in {"PASS", "FAIL", "UNKNOWN"} for v in gates.values()):
                raise ValueError("Invalid curated clinical gate")
            if entry.get("missing_bridges"):
                gates["exact_genotype_functional_rescue"] = "UNKNOWN"
            if not candidate:
                approval = {**approval, "gate": "UNKNOWN", "identity_caveat": "No eligible public molecule identity resolved"}
            rejected = bool(candidate and (candidate.get("grade") == "X" or
                                           candidate.get("status") == "exclude_direction_conflict"))
            rejection_reason = "Automatic target-action direction conflicts with the supplied gene mechanism" if rejected else None
            ranked.append({
                "drug": entry["drug"], "display_name": entry.get("display_name", entry["drug"]),
                "chembl_id": candidate["chembl_id"] if candidate else None,
                "grade": "X_direction_conflict" if rejected else "H_indirect_laboratory_hypothesis", "assessment": assessment,
                "assessment_origin": origin,
                "proposed_path": entry.get("proposed_path", []), "missing_bridges": entry.get("missing_bridges", []),
                "contrary_evidence": entry.get("contrary_evidence", []),
                "citations": entry["opened_citations"], "next_experiment": entry.get("next_experiment"),
                "independent_positive_computational_rescue_lines": 0,
                "convergence": metadata.get("convergence_note", "No independent computational rescue evidence supplied through this curation branch"),
                "approval": approval, "clinical_gates": gates,
                "rejection_reason": rejection_reason,
                "laboratory_eligible": approval["gate"] == "PASS" and not rejected,
                "clinical_eligible": approval["gate"] == "PASS" and not rejected and all(v == "PASS" for v in gates.values()),
                "clinical_eligibility_meaning": "Supplied curator gate assessment only; never software authorization to treat",
                "expected_effect_size": None, "pgx": candidate["pgx"] if candidate else [],
                "pgx_state": "marker_screen_requires_review" if candidate and candidate["pgx"] else "no_selected_drug_gene_mapping; not a safety pass",
                "pgx_role": "Review flags only; no marker coverage or genotype result demonstrates efficacy or chronic safety",
                "scenario": entry.get("scenario", metadata.get("scenario", "No scenario supplied; allele and phase applicability remain unknown")),
            })
    ranked.sort(key=lambda r: (bool(r["rejection_reason"]), -r["assessment"]["research_priority"], r["drug"]))
    return {"ranked": ranked,
            "market_approved_laboratory_hypotheses": [r["drug"] for r in ranked if r["laboratory_eligible"]],
            "unresolved_approval": [r["drug"] for r in ranked if r["approval"]["gate"] != "PASS"],
            "rejected_direction_conflicts": [r["drug"] for r in ranked if r["rejection_reason"]],
            "scoring_origin": origin + "; this engine computes caps, arithmetic, gates and sensitivity only",
            "selection_warning": metadata.get("selection_warning", "Ordinal priorities can change under sensitivity analysis; they do not establish efficacy or clinical authorization")}
