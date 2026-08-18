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

from predictr import Analysis, PlotAll, PROBABILITY_PLOT_TICKS

FAILURES = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
FAILURES_B = [1.8506941739639076, 2.2685555679846954, 2.380993183650987, 2.642404955035375,
              2.777082863078587, 2.89527127055147, 2.9099992138728927, 3.1425481097241,
              3.3758727398694406, 3.8274990886889997]
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


def test_normal_mle_line_spans_full_visible_range():
    """The MLE line must be evaluated across a wide, fixed percentile range
    (like Weibull's) rather than padded around the data's x-range, or it
    stops short of the plot edges whenever the data doesn't happen to
    cover the configured [y_min, y_max] window - see the line, not just
    the axis limits themselves."""
    data = list(np.random.default_rng(3).normal(50.0, 5.0, 10))
    a = Analysis(df=data, dist='normal', show=False, y_min=0.001, y_max=0.999)
    a.mle()
    a.plot()
    ax = plt.gca()
    mle_line = ax.get_lines()[0]
    y_of_line = mle_line.get_ydata()
    ylim = ax.get_ylim()
    assert min(y_of_line) <= ylim[0] + 1e-6
    assert max(y_of_line) >= ylim[1] - 1e-6


def test_plotall_mult_weibull_default_and_custom_range():
    a = Analysis(df=FAILURES, bounds='lrb', bounds_type='2s')
    a.mle()
    b = Analysis(df=FAILURES_B, bounds='pbb', bounds_type='2s')
    b.mle()

    fig_default = PlotAll({'a': a, 'b': b}).mult_weibull(show=False)
    ylim_default = fig_default.axes[0].get_ylim()
    expected_default = (_weibull_prob_paper(0.01), _weibull_prob_paper(0.99))
    np.testing.assert_allclose(ylim_default, expected_default)

    fig_custom = PlotAll({'a': a, 'b': b}).mult_weibull(y_min=0.05, y_max=0.95, show=False)
    ylim_custom = fig_custom.axes[0].get_ylim()
    expected_custom = (_weibull_prob_paper(0.05), _weibull_prob_paper(0.95))
    np.testing.assert_allclose(ylim_custom, expected_custom)


def test_plotall_mult_weibull_invalid_range_raises():
    a = Analysis(df=FAILURES, bounds='lrb', bounds_type='2s')
    a.mle()
    with pytest.raises(ValueError):
        PlotAll({'a': a}).mult_weibull(y_min=0.9, y_max=0.1)
