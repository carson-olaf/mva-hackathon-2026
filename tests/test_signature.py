"""A failed source refresh must never leave earlier positive evidence current."""
import hashlib
import json
import sys

import pandas as pd
import pytest

from mva_hackathon_2026 import signature


def test_failed_source_refresh_invalidates_existing_positive_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(signature, "OFFLINE", False)
    cache, output = tmp_path / "cache", tmp_path / "results"
    cache.mkdir()
    output.mkdir()
    (cache / "geneinfo_beta.txt").write_text(
        "gene_symbol\tfeature_space\n"
        + "".join(f"GENE{i}\tlandmark\n" for i in range(25))
    )
    raw = b"symbol\tCD-coefficient\nGENE1\t0.1\nGENE2\t-0.2\n"

    def cached_signature(identifier, library, advertised_sha):
        (cache / f"{identifier}.tsv").write_bytes(raw)
        return {
            "id": identifier,
            "library": library,
            "meta": {
                "pert_name": "synthetic_public_fixture",
                "persistent_id": f"https://example.invalid/{identifier}.tsv",
                "cell_line": "A549",
                "pert_time": "96 h",
                "sha256": advertised_sha,
                "size_in_bytes": len(raw) if advertised_sha else 0,
            },
        }

    valid = cached_signature("valid", signature.SHRNA, hashlib.sha256(raw).hexdigest())
    contradictory = cached_signature("contradictory", signature.CRISPR, "0" * 64)
    unadvertised = cached_signature("unadvertised", signature.MEAN_CP, "")
    (cache / "bub1b_metadata.json").write_text(json.dumps([valid, contradictory]))
    (cache / "controls_20_20261024.json").write_text("[]")
    monkeypatch.setattr(
        signature,
        "approved_compounds",
        lambda _: ({"fixture": {"signature": unadvertised}}, [unadvertised]),
    )

    def inference_must_not_run(*args, **kwargs):
        raise AssertionError("Inference ran after a source-provenance failure")

    monkeypatch.setattr(signature, "cell_consensus", inference_must_not_run)
    previous_diagnostics = {
        "evidence_status": "source_verified",
        "source_gate_passed": True,
        "ranking_bonus": 3,
        "drug_reversal_inference_validated": True,
        "top10": [{"drug": "fixture", "score": 0.9, "ranking_bonus": 3}],
    }
    original_json = json.dumps(previous_diagnostics).encode()
    (output / "signature_diagnostics.json").write_bytes(original_json)
    original_csv = b"drug,score,evidence_status,ranking_bonus\nfixture,0.9,source_verified,3\n"
    (output / "signature_ranked_approved_drugs.csv").write_bytes(original_csv)
    original_figure = b"prior figure evidence bytes"
    (output / "signature_null_calibration.png").write_bytes(original_figure)
    monkeypatch.setattr(sys, "argv", [
        "signature", "--offline", "--cache", str(cache), "--output", str(output),
        "--control-genes", "20", "--permutations", "99",
    ])
    # main() runs real cache reads, byte hashes, source validation and output
    # invalidation. Only the unrelated public-compound discovery adapter is fixed.
    with pytest.raises(SystemExit) as stopped:
        signature.main()
    assert stopped.value.code == 2
    qc = json.loads((output / "signature_source_qc.json").read_text())
    assert {row["id"]: row["status"] for row in qc["records"]} == {
        "valid": "verified", "contradictory": "mismatch", "unadvertised": "not_advertised",
    }
    current = json.loads((output / "signature_diagnostics.json").read_text())
    assert current["source_gate_passed"] is False
    assert current["current_valid"] is False
    assert current["drug_reversal_inference_validated"] is False
    assert current["evidence_status"] == "uninformative_quality_failure"
    assert current["ranking_bonus"] == current["top10"][0]["ranking_bonus"] == 0
    assert current["top10"][0]["score"] == 0.9
    table = pd.read_csv(output / "signature_ranked_approved_drugs.csv")
    assert table.loc[0, "score"] == 0.9
    assert table.loc[0, "evidence_status"] == "uninformative_quality_failure"
    assert not bool(table.loc[0, "current_valid"])
    assert table.loc[0, "ranking_bonus"] == 0
    assert not (output / "signature_null_calibration.png").exists()
    status = json.loads((output / "signature_run_status.json").read_text())
    assert status["exit_code"] == 2 and status["current_valid"] is False
    originals = {
        "signature_diagnostics.json": original_json,
        "signature_ranked_approved_drugs.csv": original_csv,
        "signature_null_calibration.png": original_figure,
    }
    for archived in status["historical_artifacts"]:
        assert (output / archived["archive"]).read_bytes() == originals[archived["original"]]
