"""Behavioral contracts for the public-only namespace census.

All records below are deliberately tiny synthetic fixtures. No private dataset,
network service, or downloaded public export is required by these tests.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/audit_target_collisions.py"
SPEC = importlib.util.spec_from_file_location("public_collision_audit", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_identifier_namespace_controls_collision_classification():
    targets = {"701": {"hgnc": {"hgnc:2021"}}}
    ncbi = {"701": {"hgnc:1149"}}

    def classify(claim, gene):
        return AUDIT.classify(
            {"gene_claim_name": claim, "gene_concept_id": gene}, targets, ncbi
        )

    assert classify("NCBIGENE:701", "hgnc:1149") == "confirmed_numeric_namespace_collision"
    assert classify("NCBIGENE:701", "hgnc:2021") == "identity_agrees"
    assert classify("IUPHAR.RECEPTOR:701", "hgnc:2021") == "identity_agrees"
    assert classify("IUPHAR.RECEPTOR:701", "hgnc:1149") == "identity_conflict_other_or_unresolved_cause"
    assert classify("NCBIGENE:701", "hgnc:9999") == "identity_conflict_other_or_unresolved_cause"
    assert classify("NCBIGENE:999", "hgnc:9999") == "unresolved_source_target"


def test_foreign_gene_namespace_is_quarantined_and_boolean_vocabularies_work():
    with pytest.raises(ValueError, match="Not a human HGNC identifier"):
        AUDIT.hgnc("MGI:103220")
    target = AUDIT.target_map([{
        "Target id": "3348", "HGNC id": "MGI:103220|7013",
        "Human Entrez Gene": "4222", "HGNC symbol": "MEOX1",
        "Human SwissProt": "P50221",
    }])["3348"]
    assert target["hgnc"] == {"hgnc:7013"}
    assert target["invalid_hgnc_values"] == {"MGI:103220"}
    assert AUDIT.hgnc("HGNC:7013") == "hgnc:7013"
    assert all(AUDIT.truth(value) for value in ["yes", "TRUE", "t", "1"])
    assert not any(AUDIT.truth(value) for value in ["", "false", "no", "0", None])


def test_normalized_drug_counts_do_not_count_salts_or_duplicate_rows_twice():
    rows = [
        {"gene_concept_id": "hgnc:1", "drug_claim_name": "IUPHAR.LIGAND:10",
         "drug_concept_id": "rxcui:shared", "interaction_types": "inhibitor"},
        {"gene_concept_id": "hgnc:1", "drug_claim_name": "IUPHAR.LIGAND:11",
         "drug_concept_id": "rxcui:shared", "interaction_types": "inhibitor"},
        {"gene_concept_id": "hgnc:2", "drug_claim_name": "IUPHAR.LIGAND:10",
         "drug_concept_id": "", "interaction_types": "inhibitor"},
    ]
    counts = AUDIT.count_edges(rows + [rows[0]])
    assert counts["rows"] == 4
    assert counts["distinct_gene_ligand_pairs"] == 3
    assert counts["distinct_gene_normalized_drug_pairs"] == 1
    assert counts["distinct_genes"] == 2


def write_tsv(path, records):
    with path.open("w", newline="") as handle:
        handle.write("# Synthetic public-schema fixture\n")
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


@pytest.fixture
def public_fixture(tmp_path):
    cache = tmp_path / "public_cache"
    cache.mkdir()
    write_tsv(cache / "gtop_targets.tsv", [
        {"Target id": "701", "HGNC id": "2021", "Human Entrez Gene": "1188",
         "HGNC symbol": "CLCNKB", "Human SwissProt": "P51801"},
        {"Target id": "1950", "HGNC id": "1149", "Human Entrez Gene": "701",
         "HGNC symbol": "BUB1B", "Human SwissProt": "O60566"},
    ])
    write_tsv(cache / "gtop_interactions.tsv", [
        {"Target ID": "701", "Target Species": "Human", "Target Ligand": "",
         "Ligand ID": lid, "Type": "inhibitor", "PubMed ID": "synthetic-reference"}
        for lid in ["10", "11"]
    ])
    write_tsv(cache / "gtop_ligands.tsv", [
        {"Ligand ID": "10", "Name": "Synthetic active form", "Approved": "yes", "Withdrawn": ""},
        {"Ligand ID": "11", "Name": "Synthetic withdrawn salt", "Approved": "yes", "Withdrawn": "yes"},
    ])
    with gzip.open(cache / "human_gene_info.gz", "wt") as handle:
        handle.write("#tax_id\tGeneID\tdbXrefs\n9606\t701\tHGNC:HGNC:1149\n")
        handle.write("10090\t701\tHGNC:HGNC:9999\n")

    def claim(prefix, target, gene_id, gene_name, ligand):
        return {
            "gene_claim_name": f"{prefix}:{target}", "gene_concept_id": f"hgnc:{gene_id}",
            "gene_name": gene_name, "drug_claim_name": f"IUPHAR.LIGAND:{ligand}",
            "drug_concept_id": "rxcui:shared", "drug_name": "Synthetic normalized drug",
            "drug_is_approved": "true", "interaction_source_db_name": "GuideToPharmacology",
            "interaction_source_db_version": "synthetic", "interaction_types": "inhibitor",
        }

    write_tsv(cache / "dgidb_old_interactions.tsv", [
        claim("NCBIGENE", "701", "1149", "BUB1B", "10"),
        claim("NCBIGENE", "701", "1149", "BUB1B", "11"),
        claim("NCBIGENE", "1950", "1149", "BUB1B", "10"),
        claim("NCBIGENE", "701", "9999", "SYNTHETIC_OTHER", "10"),
    ])
    write_tsv(cache / "dgidb_2026-10_interactions.tsv", [
        claim("IUPHAR.RECEPTOR", "701", "2021", "CLCNKB", "10"),
        claim("IUPHAR.RECEPTOR", "701", "2021", "CLCNKB", "11"),
    ])
    (cache / "orphanet_genes.xml").write_text("""<JDBOR date="synthetic">
      <Disorder><OrphaCode>123</OrphaCode><Name>Synthetic test disorder</Name>
        <DisorderGroup><Name>Disorder</Name></DisorderGroup>
        <DisorderGeneAssociationList><DisorderGeneAssociation>
          <SourceOfValidation>synthetic-reference</SourceOfValidation>
          <Gene><Symbol>BUB1B</Symbol><ExternalReferenceList><ExternalReference>
            <Source>HGNC</Source><Reference>1149</Reference>
          </ExternalReference></ExternalReferenceList></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation></DisorderGeneAssociationList>
      </Disorder></JDBOR>""")
    records = [{"filename": path.name, "complete": True,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
               for path in cache.iterdir()]
    (cache / "manifest.json").write_text(json.dumps({"requests": records}))
    return cache


def test_offline_cli_measures_observed_errors_and_the_repaired_release(public_fixture, tmp_path):
    output = tmp_path / "audit_results"
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--cache", str(public_fixture), "--output", str(output),
         "--snapshots", "old", "2026-10", "--exposure-snapshot", "old"],
        check=True, capture_output=True, text=True,
    )
    summary = json.loads(completed.stdout)
    old = summary["snapshots"]["old"]
    assert old["gtop_source_rows"] == 4
    assert old["confirmed_numeric_collisions"]["distinct_gene_ligand_pairs"] == 2
    assert old["confirmed_numeric_collisions"]["distinct_gene_normalized_drug_pairs"] == 1
    assert old["classifications"]["identity_agrees"] == 1
    assert old["classifications"]["identity_conflict_other_or_unresolved_cause"] == 1
    assert summary["snapshots"]["2026-10"]["confirmed_numeric_collisions"]["rows"] == 0
    assert summary["snapshots"]["2026-10"]["classifications"]["identity_agrees"] == 2
    assert summary["exposed_orpha_entities"] == 1
    with (output / "rare_disease_exposure.csv").open() as handle:
        disease = next(csv.DictReader(handle))
    assert int(disease["unique_approved_flag_ligands"]) == 1
    with (output / "identity_conflicts.csv").open() as handle:
        conflicts = list(csv.DictReader(handle))
    assert len(conflicts) == 3
    assert {r["correct_pair_present_2026-10"] for r in conflicts} == {"True"}


def test_tampered_public_cache_is_rejected_before_counting(public_fixture, tmp_path):
    path = public_fixture / "gtop_targets.tsv"
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="Hash mismatch: gtop_targets.tsv"):
        AUDIT.audit(public_fixture, tmp_path / "results", ["old", "2026-10"], "old")


def test_empty_findings_replace_a_previous_run(tmp_path):
    path = tmp_path / "findings.csv"
    AUDIT.write_csv(path, [{"gene": "previous finding"}])
    AUDIT.write_csv(path, [])
    with path.open() as handle:
        assert list(csv.DictReader(handle)) == []
