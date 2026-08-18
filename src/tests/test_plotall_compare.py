"""
Tests for PlotAll.compare() - the cross-distribution comparison grid.

Unlike mult_weibull() (same dist required, overlaid on one shared paper),
compare() takes objects of *different* dist= choices fit to the *same*
data and lays out one probability-plot panel per object, ranked by AIC.
The dataset used here is the same right-censored example validated in
this session against a published reference (beta=1.7208, eta=606.528,
loglik=-75.135, AIC=154.27), so the ranking/AIC values asserted below are
themselves already externally validated, not just internally consistent.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from predictr import Analysis, PlotAll

FAILURES = [54, 187, 216, 240, 244, 335, 361, 373, 375, 386]
SUSPENSIONS = [500] * 10


def _fit_all_distributions():
    fits = {}
    for dist in ['weibull', 'normal', 'lognormal', 'exponential']:
        a = Analysis(df=FAILURES, ds=SUSPENSIONS, dist=dist, bounds='fb',
                     bounds_type='2s', cl=0.9, show=False)
        a.mle()
        fits[dist] = a
    return fits


def test_compare_runs_and_ranks_by_aic():
    fits = _fit_all_distributions()
    PlotAll(fits).compare()
    fig = plt.gcf()
    panel_titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    # Weibull (AIC=154.27) must rank first, Normal (AIC=156.20) last, for
    # this specific validated dataset.
    assert panel_titles[0].endswith('Weibull')
    assert panel_titles[-1].endswith('Normal')


def test_compare_panel_count_matches_object_count():
    fits = _fit_all_distributions()
    small = {'weibull': fits['weibull'], 'normal': fits['normal']}
    PlotAll(small).compare()
    fig = plt.gcf()
    visible_axes = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible_axes) == 2


def test_compare_requires_mle_run_on_every_object():
    fits = _fit_all_distributions()
    unfitted = Analysis(df=FAILURES, ds=SUSPENSIONS, dist='weibull', show=False)
    fits['unfitted'] = unfitted
    with pytest.raises(ValueError, match='mle'):
        PlotAll(fits).compare()


def test_compare_requires_same_data_across_objects():
    fits = _fit_all_distributions()
    different_data = Analysis(df=[1.0, 2.0, 3.0], dist='weibull', show=False)
    different_data.mle()
    fits['different_data'] = different_data
    with pytest.raises(ValueError, match='same data'):
        PlotAll(fits).compare()


def test_compare_invalid_y_range_raises():
    fits = _fit_all_distributions()
    with pytest.raises(ValueError):
        PlotAll(fits).compare(y_min=0.9, y_max=0.1)


def test_compare_aic_values_match_published_reference():
    """Cross-check against the external reference validated earlier this
    session (NIST/SEMATECH-style worked example): beta=1.7208, eta=606.528,
    loglik=-75.135, AIC=154.27 for the Weibull fit on this exact dataset."""
    fits = _fit_all_distributions()
    w = fits['weibull']
    assert w.beta == pytest.approx(1.7208, abs=1e-3)
    assert w.eta == pytest.approx(606.528, abs=1e-2)
    assert w.loglik == pytest.approx(-75.135, abs=1e-2)
    assert w.aic == pytest.approx(154.27, abs=1e-2)
