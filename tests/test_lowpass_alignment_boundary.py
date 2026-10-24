"""Private-output and overwrite boundaries are enforced before native tools run."""
import pytest

from mva_hackathon_2026.lowpass.align import align


def test_alignment_rejects_public_outputs_external_reads_and_existing_bam(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "work"
    work.mkdir()
    r1, r2 = work / "left.fq", work / "right.fq"
    for path in (r1, r2):
        path.write_text("synthetic placeholder; must not be read")
    reference = work / "reference.fa"
    with pytest.raises(ValueError, match="outputs must stay"):
        align(reference, r1, r2, work / "attempt", tmp_path / "results.bam")
    with pytest.raises(ValueError, match="Read inputs"):
        align(reference, tmp_path / "external.fq", r2, work / "attempt", work / "out.bam")
    output = work / "out.bam"
    output.write_bytes(b"preserved prior artifact")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        align(reference, r1, r2, work / "attempt", output)
    assert output.read_bytes() == b"preserved prior artifact"
    assert not (work / "attempt").exists()
    assert not (tmp_path / "results.bam").exists()
