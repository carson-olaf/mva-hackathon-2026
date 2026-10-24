"""Behavioral tests against frozen real public API responses, with network denied."""

import json
from pathlib import Path

import requests

from mva_repurpose.cache import PublicCache
from mva_repurpose.cli import main
from mva_repurpose.engine import classify_direction, run


CACHE = Path(__file__).resolve().parents[1] / "examples" / "public_api_cache"


def no_network(*args, **kwargs):
    raise AssertionError("An offline reproducibility test attempted a network request")


def test_bub1b_end_to_end_abstains_when_only_pathway_proximity_exists(tmp_path, monkeypatch):
    monkeypatch.setattr(requests.Session, "request", no_network)
    prefix = tmp_path / "bub1b"
    assert main([
        "--gene", "BUB1B", "--variant", "ENST00000287598.11:c.2210T>G",
        "--variant", "ENST00000287598.11:c.3006T>G",
        "--hpo", "HP:0002859", "--hpo", "HP:0000121",
        "--hpo", "HP:0004322", "--hpo", "HP:0001508",
        "--hpo", "HP:0003202", "--hpo", "HP:0001622",
        "--hpo", "HP:0001518",
        "--mechanism", "loss_of_function",
        "--disease", "mosaic variegated aneuploidy syndrome",
        "--cache", str(CACHE), "--offline", "--output", str(prefix),
    ]) == 0
    result = json.loads(prefix.with_suffix(".json").read_text())
    assert len(result["ranking"]) >= 10
    assert {r["drug"] for r in result["ranking"]} >= {"COLCHICINE", "PACLITAXEL"}
    assert all(float(r["approval"]["max_phase"]) == 4 and not r["approval"]["withdrawn_flag"]
               for r in result["ranking"])
    assert all(r["grade"] in {"E", "X"} for r in result["ranking"])
    assert "ATALUREN" not in {r["drug"] for r in result["ranking"]}
    assert any(r["drug"] == "ATALUREN" and r.get("regulatory_check") for r in result["excluded_not_eligible"])
    assert result["supported_hypotheses"] == []
    assert result["abstention"]
    assert len(result["hpo"]) == 7
    assert all(term["label"] != "unresolved" for term in result["hpo"])
    assert prefix.with_suffix(".csv").stat().st_size > 100
    assert prefix.with_suffix(".md").stat().st_size > 100
    # No phenotype-specific public query was needed to resolve the input.
    assert all("HP:" not in json.dumps(item["request"]) for item in result["api_provenance"])


def test_tsc2_discovers_published_everolimus_positive_control_and_rejects_wrong_disease(monkeypatch):
    monkeypatch.setattr(requests.Session, "request", no_network)
    common = dict(gene="TSC2", variants=[], hpo=["HP:0001250"], mechanism="loss_of_function")
    result = run(PublicCache(CACHE, offline=True), disease="tuberous sclerosis", **common)
    names = [r["drug"] for r in result["ranking"]]
    assert names.index("EVEROLIMUS") < 3
    everolimus = next(r for r in result["ranking"] if r["drug"] == "EVEROLIMUS")
    assert everolimus["grade"] == "L"
    assert everolimus["status"] == "known_disease_label_match"
    assert not everolimus["approval"]["current_approval_verified"]
    assert result["supported_hypotheses"] == []
    assert everolimus["disease_matches"]
    assert any(edge["mechanism"] for edge in everolimus["mechanism_edges"])
    assert not everolimus["supplied_evidence"]  # The drug name was not supplied as an answer.
    negative = run(PublicCache(CACHE, offline=True), disease="mosaic variegated aneuploidy syndrome", **common)
    same_drug = next(r for r in negative["ranking"] if r["drug"] == "EVEROLIMUS")
    assert same_drug["grade"] == "E"
    assert same_drug["score"] < everolimus["score"]
    assert negative["supported_hypotheses"] == []
    assert negative["abstention"]
    normal = run(PublicCache(CACHE, offline=True), disease="tuberous sclerosis", pgx=[{
        "gene": "CYP3A5", "phenotype": "normal metabolizer", "candidate_drugs": "EVEROLIMUS",
    }], **common)
    assert next(r for r in normal["ranking"] if r["drug"] == "EVEROLIMUS")["components"]["pgx"] == 0
    review = run(PublicCache(CACHE, offline=True), disease="tuberous sclerosis", pgx=[{
        "gene": "CYP3A5", "phenotype": "intermediate metabolizer", "candidate_drugs": "EVEROLIMUS",
    }], **common)
    assert next(r for r in review["ranking"] if r["drug"] == "EVEROLIMUS")["components"]["pgx"] < 0


def test_direction_edges_and_clinical_pgx_column_respect_their_boundaries(monkeypatch):
    """Synthetic boundary inputs exercise the real ranking/output contract."""
    from mva_repurpose import engine

    assert classify_direction("Receptor inverse agonist") == "inhibit"
    assert classify_direction("Receptor partial agonist") == "unknown"
    proteins = {"PRIMARY": {"direct": True, "pathways": []}, "NEIGHBOR": {"direct": False, "pathways": []}}
    targets = {
        "PRIMARY_TARGET": {"target_components": [{"accession": "PRIMARY"}]},
        "NEIGHBOR_TARGET": {"target_components": [{"accession": "NEIGHBOR"}]},
    }
    mechanisms = [
        {"molecule_chembl_id": "SYNTHETIC", "target_chembl_id": "PRIMARY_TARGET", "mechanism_of_action": "Receptor agonist"},
        {"molecule_chembl_id": "SYNTHETIC", "target_chembl_id": "NEIGHBOR_TARGET", "mechanism_of_action": "Receptor inverse agonist"},
        {"molecule_chembl_id": "SYNTHETIC", "target_chembl_id": "PRIMARY_TARGET", "mechanism_of_action": "Receptor partial agonist"},
    ]
    monkeypatch.setattr(engine, "graph", lambda *args: ({"id": "SYNTHETIC"}, [], proteins, targets, mechanisms, {}))

    class SyntheticCache:
        used = {}

        def pages(self, resource, collection, params):
            if resource == "molecule":
                return [{"molecule_chembl_id": "SYNTHETIC", "pref_name": "SYNTHETIC ACTIVE SULFATE",
                         "max_phase": 4, "withdrawn_flag": False,
                         "molecule_synonyms": [{"molecule_synonym": "SYNTHETIC ACTIVE"}]}]
            assert resource == "drug_indication"
            return []

    result = run(SyntheticCache(), gene="TEST", variants=[], hpo=[], mechanism="loss_of_function", pgx=[
        {"gene": "SYNTHETIC_PGX", "phenotype": "normal metabolizer", "candidate_drugs_affected": "synthetic active/another agent (context)", "clinical_action": "Review assay limitations"},
        {"gene": "UNMATCHED", "phenotype": "unknown", "candidate_drugs_affected": "unexpanded drug class"},
    ])
    candidate = result["ranking"][0]
    assert [e["causal_direction_supported"] for e in candidate["mechanism_edges"]] == [True, False, False]
    assert candidate["pgx"][0]["gene"] == "SYNTHETIC_PGX"
    assert candidate["pgx"][0]["recommendation"] == "Review assay limitations"
    assert candidate["components"]["pgx"] == 0
    assert result["pgx_integration"]["matched_input_rows"] == 1
    assert result["pgx_integration"]["unmatched_rows"][0]["gene"] == "UNMATCHED"


def test_curated_indirect_paths_keep_caps_and_clinical_gates(monkeypatch):
    """Real public replay cannot turn supplied priorities into rescue evidence."""
    from copy import deepcopy
    from mva_repurpose.curation import score_assessment

    monkeypatch.setattr(requests.Session, "request", no_network)
    root = Path(__file__).resolve().parents[1]
    bundle = json.loads((root / "examples/bub1b_curated_evidence.json").read_text())
    # Even an overstated supplied tag cannot close a visibly missing bridge.
    bundle["evidence"][0]["direction"] = "rescue"
    bundle["evidence"][0]["functional_rescue_demonstrated"] = True
    result = run(PublicCache(root / "examples/bub1b_curated_cache", offline=True, fallbacks=(CACHE,)),
                 gene="BUB1B", variants=bundle["curation"]["applicable_variants"], hpo=[], mechanism="loss_of_function",
                 evidence=bundle, disease="mosaic variegated aneuploidy syndrome")
    rows = result["curated_hypotheses"]["ranked"]
    assert {r["drug"] for r in rows} == {"ARIMOCLOMOL", "AZITHROMYCIN", "NIACIN"}
    assert all(r["assessment"]["effective_components"]["E"] <= 1 for r in rows)
    assert all(r["missing_bridges"] and r["clinical_eligible"] is False for r in rows)
    assert all(r["approval"]["gate"] == "UNKNOWN" for r in rows)
    assert all(r["independent_positive_computational_rescue_lines"] == 0 for r in rows)
    assert not result["supported_hypotheses"]
    altered = deepcopy(bundle["evidence"][0])
    altered["assessed_components"]["E"] = 5
    capped = score_assessment(altered, bundle["curation"]["weights"])
    assert capped["effective_components"]["E"] == 1
    before = capped["research_priority"]
    altered["assessed_components"]["G"] += 1
    assert score_assessment(altered, bundle["curation"]["weights"])["research_priority"] == before + bundle["curation"]["weights"]["G"]
    import pytest
    from mva_repurpose.engine import resolve_gene
    with pytest.raises(ValueError, match="at most 10 named variants"):
        resolve_gene(None, "BUB1B", ["BUB1B:p.Arg1Lys"] * 11)
    with pytest.raises(ValueError, match="curated variant-specific assessment"):
        run(PublicCache(root / "examples/bub1b_curated_cache", offline=True, fallbacks=(CACHE,)),
            gene="BUB1B", variants=["BUB1B:p.Arg1Lys"], hpo=[], mechanism="loss_of_function", evidence=bundle)
    from mva_repurpose import curation
    # Product approval and supplied all-PASS clinical gates cannot override a
    # target-action conflict found by the automatic mechanism layer.
    entry = deepcopy(bundle["evidence"][0])
    entry.update(opened_citations=[], missing_bridges=[], clinical_gates={k: "PASS" for k in curation.GATES})
    rejected = deepcopy(next(r for r in result["ranking"] if r["drug"] == entry["drug"]))
    rejected.update(grade="X", status="exclude_direction_conflict")
    monkeypatch.setattr(curation, "approval_for", lambda *args: {"gate": "PASS"})
    generic_bundle = {"curation": {"weights": bundle["curation"]["weights"]}}
    guarded = curation.curated_hypotheses({entry["drug"]: [entry]}, generic_bundle, None, [rejected])
    assert guarded["market_approved_laboratory_hypotheses"] == []
    assert guarded["ranked"][0]["grade"] == "X_direction_conflict"
    assert guarded["ranked"][0]["rejection_reason"]
    assert guarded["ranked"][0]["clinical_eligible"] is False
    assert "arimoclomol" not in guarded["selection_warning"].casefold()
    assert "Pro" not in guarded["scoring_origin"]
    assert guarded["ranked"][0]["scenario"].startswith("No scenario supplied")
