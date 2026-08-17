"""
Characterization tests for the existing (pre-refactor) Weibull behavior.

These pin down today's numeric output for mle()/mrr() plus Fisher and
Beta-Binomial bounds, using the exact datasets from docs/classes.md. They
exist to catch accidental behavior changes while predictr is restructured
to support multiple distributions (dist= plugin architecture) - not to
validate correctness against an external reference.
"""
import numpy as np
import pytest

from predictr import Analysis

FAILURES = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
SUSPENSIONS = [1.9, 2.0, 2.0]


def test_mle_uncensored():
    a = Analysis(df=FAILURES, bounds='fb', show=False)
    a.mle()
    assert a.beta == pytest.approx(2.511134261086927)
    assert a.eta == pytest.approx(1.1136702424830867)


def test_mle_uncensored_fisher_bounds():
    a = Analysis(df=FAILURES, bounds='fb', show=False)
    a.mle()
    expected_lower = [0.01486765, 0.02261472, 0.02889973, 0.03935792, 0.04823643,
                       0.05984536, 0.09102804, 0.11639564, 0.15882946, 0.19514885,
                       0.24316787, 0.37625443, 0.49066944, 0.59730177, 0.70098116,
                       0.80498416, 0.83885957, 0.91274685, 1.03012305, 1.17372223,
                       1.27655538, 1.44046786, 1.59394762]
    expected_upper = [0.34047517, 0.388925, 0.42051836, 0.46417928, 0.49553517,
                       0.53126818, 0.60908331, 0.66058682, 0.7331727, 0.78667859,
                       0.84960306, 0.99818206, 1.11206575, 1.21610863, 1.32139261,
                       1.43710071, 1.47812314, 1.57532537, 1.75885648, 2.05322059,
                       2.32801478, 2.90571867, 3.62686569]
    np.testing.assert_allclose(a.bounds_lower, expected_lower, rtol=1e-6)
    np.testing.assert_allclose(a.bounds_upper, expected_upper, rtol=1e-6)


def test_mle_censored():
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, bounds=None, show=False)
    a.mle()
    assert a.beta == pytest.approx(1.6924809036741737)
    assert a.eta == pytest.approx(1.7783001813795736)


def test_mrr_uncensored():
    a = Analysis(df=FAILURES, bounds='bbb', show=False)
    a.mrr()
    assert a.beta == pytest.approx(2.1231598259170834)
    assert a.eta == pytest.approx(1.1303814671874557)


def test_mrr_uncensored_bbb_bounds():
    a = Analysis(df=FAILURES, bounds='bbb', show=False)
    a.mrr()
    expected_lower = [0.008512444610847119, 0.06284989170835438, 0.15316111797522317,
                       0.27133837251975246, 0.41819659074797405, 0.6069622310029172]
    expected_upper = [0.39303776899708265, 0.5818034092520259, 0.7286616274802475,
                       0.8468388820247768, 0.9371501082916456, 0.9914875553891529]
    np.testing.assert_allclose(np.array(a.bounds_lower, dtype=float), expected_lower, rtol=1e-6)
    np.testing.assert_allclose(np.array(a.bounds_upper, dtype=float), expected_upper, rtol=1e-6)


def test_mrr_censored():
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, bounds=None, show=False)
    a.mrr()
    assert a.beta == pytest.approx(1.8019397878368828)
    assert a.eta == pytest.approx(1.6135654263192374)
