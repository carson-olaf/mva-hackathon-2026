import numpy as np

from mva_hackathon_2026.mosaic import (
    broadening_from_fraction, excess_moment, fraction_from_broadening,
)


def test_binomial_null_and_mosaic_signal_recovery():
    """The boundary statistic removes sampling broadening and recovers a true shift."""
    rng = np.random.default_rng(42)
    n = 35
    null = excess_moment(rng.binomial(n, .5, 200_000), n)
    assert abs(null.mean()) < 6e-5
    for event in ("gain", "loss"):
        expected = broadening_from_fraction(.30, event)
        delta = np.sqrt(expected)
        shifted = excess_moment(rng.binomial(n, .5 + delta, 200_000), n)
        assert abs(shifted.mean() - expected) < 1e-4
        assert abs(fraction_from_broadening(expected, event) - .30) < 1e-10


def test_impossible_gain_is_not_reported_as_a_valid_detection_fraction():
    assert fraction_from_broadening(.10, "gain") is None


def test_full_regional_screen_null_and_injected_chromosome():
    from mva_hackathon_2026.mosaic import analyze
    rng = np.random.default_rng(121)
    keys = [(str(chrom), arm, block) for chrom in range(1,23) for arm in ('p','q') for block in range(6)]
    hist = np.zeros((len(keys),101,101),dtype=np.int32)
    for i,_ in enumerate(keys):
        hist[i,40,:41]=np.bincount(rng.binomial(40,.5,3000),minlength=41)
    regions, summary = analyze(keys,hist,iterations=5000,seed=41)
    assert summary['tests']==66
    assert not any(r.get('exploratory_outlier') for r in regions)
    delta=np.sqrt(broadening_from_fraction(.3,'gain'))
    for i,key in enumerate(keys):
        if key[0]=='15':
            hist[i,40,:41]=np.bincount(rng.binomial(40,.5+delta,3000),minlength=41)
    regions,_=analyze(keys,hist,iterations=5000,seed=41)
    assert all(next(r for r in regions if r['region']==region)['exploratory_outlier'] for region in ('chr15','chr15p','chr15q'))
