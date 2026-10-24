"""Behavior at the GC-corrected chromosome-screen boundary, without patient data."""
import numpy as np
import pandas as pd
import pytest

from mva_hackathon_2026.lowpass.analyze import analyze, prepare
from mva_hackathon_2026.lowpass.provenance import complete, reusable, signature


def test_completion_marker_rejects_interruption_stale_inputs_parameters_code_and_output(tmp_path):
    source, code, output, marker = [tmp_path / name for name in ("input", "code.py", "depth.tsv", "complete.json")]
    source.write_text("synthetic input")
    code.write_text("method v1")
    output.write_text("complete synthetic output")
    key = signature([source], {"mapq": 30}, [code])
    assert not reusable(marker, key, [output])  # Bare output from interrupted stage.
    complete(marker, key, [output])
    assert reusable(marker, key, [output])
    assert not reusable(marker, signature([source], {"mapq": 20}, [code]), [output])
    source.write_text("changed input")
    assert not reusable(marker, signature([source], {"mapq": 30}, [code]), [output])
    source.write_text("synthetic input")
    code.write_text("method v2")
    assert not reusable(marker, signature([source], {"mapq": 30}, [code]), [output])
    code.write_text("method v1")
    output.write_text("truncated")
    assert not reusable(marker, key, [output])
    output.unlink()
    assert not reusable(marker, key, [output])
    marker.write_text("[]")
    assert not reusable(marker, key, [output])


def test_incomplete_or_duplicate_genome_coordinates_cannot_enter_inference():
    metadata = pd.DataFrame([(f"chr{i}", 0, 100_000, "q", .5, 0.) for i in range(1, 23)],
                            columns=["chrom", "start", "end", "arm", "gc", "n_fraction"])
    coverage = metadata[["chrom", "start", "end"]].copy()
    coverage["sample"] = coverage["normal"] = 100_000
    with pytest.raises(ValueError, match="all22"):
        prepare(metadata.iloc[:-1], coverage.iloc[:-1])
    with pytest.raises(ValueError, match="match exactly"):
        prepare(metadata, coverage.iloc[:-1])
    with pytest.raises(ValueError, match="Duplicate"):
        prepare(pd.concat([metadata, metadata.iloc[:1]]), pd.concat([coverage, coverage.iloc[:1]]))


def test_gc_corrected_null_and_known_copy_gain_survive_spatial_sensitivity():
    rng = np.random.default_rng(29)
    reference, depths = [], []
    for chromosome in range(1, 23):
        for number in range(500):
            start = number * 100_000
            gc = rng.uniform(.3, .7)
            coordinate = (f"chr{chromosome}", start, start + 100_000)
            reference.append((*coordinate, "p" if number < 250 else "q", gc, 0.))
            # Deliberately opposing library GC effects must not become CN calls.
            normal = 1.1 * np.exp(2 * (gc - .5) + rng.normal(0, .04))
            sample = 1.1 * np.exp(-1.5 * (gc - .5) + rng.normal(0, .04))
            depths.append((*coordinate, sample * 100_000, normal * 100_000))
    metadata = pd.DataFrame(reference, columns=["chrom", "start", "end", "arm", "gc", "n_fraction"])
    coverage = pd.DataFrame(depths, columns=["chrom", "start", "end", "sample", "normal"])
    null_frame, _, qc = prepare(metadata, coverage)
    assert abs(qc["sample_gc_fit"]["pearson_gc_logdepth_before"]) > .8
    assert abs(qc["sample_gc_fit"]["pearson_gc_logdepth_after"]) < .05
    _, null = analyze(null_frame, iterations=5000)
    assert null["exploratory_outlier_regions"] == []

    # One extra copy in 30% of cells gives diploid-relative CN 2.3, ratio 1.15.
    coverage.loc[coverage.chrom == "chr7", "sample"] *= 1.15
    signal_frame, _, signal_qc = prepare(metadata, coverage)
    profile = signal_qc["descriptive_gc_normalized_depth_by_chromosome"]["chr7"]
    assert 1.10 < profile["sample"] < 1.20
    assert .98 < profile["normal"] < 1.02
    for block_mb in (1, 5):
        rows, summary = analyze(signal_frame, block_mb=block_mb, iterations=5000)
        target = next(row for row in rows if row["region"] == "chr7")
        assert target["exploratory_outlier"]
        assert target["bonferroni_p"] < .05
        assert abs(target["relative_cn"] - 2.3) < .04
        assert target["relative_cn_block_bootstrap_ci95_low"] < target["relative_cn"]
        assert target["relative_cn_block_bootstrap_ci95_high"] > target["relative_cn"]
        assert target["gain_fraction_power90"] < .30
        assert summary["absolute_ploidy_or_total_mva_burden_estimated"] is False
