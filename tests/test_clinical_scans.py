from mva_hackathon_2026.clinical_scans import (
    high_confidence_clinvar, local_clinvar_frequencies, marker_phenotype, query_marker,
)


def test_absent_vcf_record_is_unknown_not_homozygous_reference(monkeypatch):
    monkeypatch.setattr("subprocess.check_output", lambda *args, **kwargs: "")
    marker = {"chrom": "7", "pos": 1, "ref": "T", "alt": "C", "complex": False}
    result = query_marker("synthetic.vcf.gz", marker)
    assert result["status"] == "absent_variant_only_vcf"
    assert result["genotype"] == "unknown"
    assert "Indeterminate" in marker_phenotype("CYP3A5", {"rs776746": result})


def test_secondary_filter_rejects_conflict_weak_evidence_and_benign():
    good = {"CLNSIG": "Pathogenic", "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts"}
    assert high_confidence_clinvar(good)
    assert not high_confidence_clinvar({**good, "CLNSIGCONF": "Likely_benign(1)"})
    assert not high_confidence_clinvar({**good, "CLNREVSTAT": "criteria_provided,_single_submitter"})
    assert not high_confidence_clinvar({**good, "CLNSIG": "Benign"})


def test_local_population_frequency_requires_exact_allele(monkeypatch, tmp_path):
    (tmp_path / "clinical_clinvar.vcf.gz").touch()
    (tmp_path / "clinical_clinvar.vcf.gz.tbi").touch()
    monkeypatch.setattr("subprocess.check_output", lambda *args, **kwargs:
                        "1\t10\tA\tC\t0.2\t.\n1\t10\tA\tT\t0.9\t.\n")
    variants = [{"chrom": "1", "pos": 10, "ref": "A", "alt": "C"}]
    frequencies = local_clinvar_frequencies(variants, tmp_path)
    assert frequencies == {("1", 10, "A", "C"): .2}
