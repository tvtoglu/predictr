"""
Tests for the configurable probability-plot y-axis range (y_min=/y_max=).

The tricky part this covers: fill_between()/fill_betweenx() (used to shade
the confidence bounds) makes matplotlib auto-expand the y-view when it's
added, regardless of an earlier plt.ylim() call - so the final visible
range has to be set *after* all plotting, or a bounds plot silently ignores
y_min/y_max. These tests render real figures and read back the actual
Axes ylim rather than just checking that no exception was raised.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.stats import norm

from predictr import Analysis, PROBABILITY_PLOT_TICKS

FAILURES = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
NORMAL_DATA = list(np.random.default_rng(1).normal(100.0, 15.0, 25))


def _weibull_prob_paper(p):
    return np.log(-np.log(1 - p))


@pytest.mark.parametrize('bad_kwargs', [
    dict(y_min=0.0, y_max=0.99),
    dict(y_min=0.01, y_max=1.0),
    dict(y_min=0.5, y_max=0.4),
    dict(y_min=0.5, y_max=0.5),
])
def test_invalid_y_range_raises(bad_kwargs):
    with pytest.raises(ValueError):
        Analysis(df=FAILURES, show=False, **bad_kwargs)


def test_default_y_range_is_001_099():
    a = Analysis(df=FAILURES, show=False)
    assert a.y_min == 0.01
    assert a.y_max == 0.99


def test_weibull_plot_respects_default_range_with_bounds():
    a = Analysis(df=FAILURES, bounds='fb', show=False)
    a.mle()
    a.plot()
    ylim = plt.gca().get_ylim()
    expected = (_weibull_prob_paper(0.01), _weibull_prob_paper(0.99))
    np.testing.assert_allclose(ylim, expected)


def test_weibull_plot_respects_custom_range_with_bounds():
    a = Analysis(df=FAILURES, bounds='fb', show=False, y_min=0.05, y_max=0.95)
    a.mle()
    a.plot()
    ylim = plt.gca().get_ylim()
    expected = (_weibull_prob_paper(0.05), _weibull_prob_paper(0.95))
    np.testing.assert_allclose(ylim, expected)


def test_weibull_plot_mrr_respects_custom_range_with_bounds():
    a = Analysis(df=FAILURES, bounds='bbb', show=False, y_min=0.02, y_max=0.98)
    a.mrr()
    a.plot_mrr()
    ylim = plt.gca().get_ylim()
    expected = (_weibull_prob_paper(0.02), _weibull_prob_paper(0.98))
    np.testing.assert_allclose(ylim, expected)


def test_normal_plot_respects_custom_range_with_bounds():
    a = Analysis(df=NORMAL_DATA, dist='normal', bounds='fb', show=False,
                 y_min=0.001, y_max=0.999)
    a.mle()
    a.plot()
    ylim = plt.gca().get_ylim()
    expected = (norm.ppf(0.001), norm.ppf(0.999))
    np.testing.assert_allclose(ylim, expected)


def test_tick_candidates_are_filtered_not_regenerated():
    """The step spacing (PROBABILITY_PLOT_TICKS) must stay the same set of
    values regardless of range - narrowing y_min/y_max should only drop
    ticks outside the range, not change the remaining ones."""
    a = Analysis(df=FAILURES, show=False, y_min=0.05, y_max=0.95)
    a.mle()
    a.plot()
    yticks = plt.gca().get_yticks()
    visible_probs = 1 - np.exp(-np.exp(yticks))
    expected = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= 0.05)
                                       & (PROBABILITY_PLOT_TICKS <= 0.95)]
    np.testing.assert_allclose(np.sort(visible_probs), np.sort(expected), atol=1e-6)
