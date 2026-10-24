"""Pre-GC injection recovery on synthetic, variable-GC covered-base profiles."""
import numpy as np
import pandas as pd

from mva_hackathon_2026.lowpass.analyze import analyze, prepare
from mva_hackathon_2026.lowpass.recovery import block_statistics, recover


def test_pre_gc_depth_injection_recovers_known_shift_and_matches_main_statistic():
    rng = np.random.default_rng(101)
    metadata, coverage = [], []
    for chromosome in range(1, 23):
        for number in range(100):
            coordinate = (f"chr{chromosome}", number * 100_000, (number + 1) * 100_000)
            gc = rng.uniform(.3, .7)
            metadata.append((*coordinate, "p" if number < 50 else "q", gc, 0.))
            normal = np.exp(2 * (gc - .5) + rng.normal(0, .005))
            sample = np.exp(-1.5 * (gc - .5) + rng.normal(0, .005))
            coverage.append((*coordinate, 100_000 * sample, 100_000 * normal))
    metadata = pd.DataFrame(metadata, columns=["chrom", "start", "end", "arm", "gc", "n_fraction"])
    coverage = pd.DataFrame(coverage, columns=["chrom", "start", "end", "sample", "normal"])
    original = coverage.copy()
    baseline, _, _ = prepare(metadata, coverage)
    statistics = block_statistics(baseline)
    main_rows, _ = analyze(baseline, iterations=5000)
    for chromosome in ("chr7", "chr22"):
        main_statistic = next(row["log2_ratio"] for row in main_rows if row["region"] == chromosome)
        assert abs(main_statistic - statistics[chromosome]["statistic"]) < 1e-12
        assert abs(main_statistic) < .01
    rows = recover(metadata, coverage, fractions=(.1, .3), chromosomes=("chr7", "chr22"))
    assert len(rows) == 8
    assert all(row["status"] == "evaluated" for row in rows)
    assert all(.94 < row["recovery_ratio"] < 1.06 for row in rows)
    assert all(row["retained_genome_bins_change"] == 0 for row in rows)
    assert all(row["retained_target_bins_change"] == 0 for row in rows)
    assert all(np.sign(row["observed_statistic_delta"]) == (1 if row["direction"] == "gain" else -1) for row in rows)
    pd.testing.assert_frame_equal(original, coverage)
