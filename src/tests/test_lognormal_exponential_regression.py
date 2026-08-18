"""
Tests for LogNormal (dist='lognormal') and Exponential (dist='exponential').

LogNormal is fit as Normal on ln(data) - mu/sigma are the Normal
parameters of ln(T) - so its tests mirror test_normal_regression.py's
approach (finite-difference gradient check, closed-form uncensored
standard errors) but on log-transformed data, plus a check that bounds
stay positive after the exp() round-trip.

Exponential has a fully closed-form MLE (theta_hat = total time on test /
number of failures) even under censoring, so there's no root-finder to
validate - instead these check theta_hat against the total-time-on-test
formula directly and the Fisher bounds against the closed-form
Var(theta_hat) = theta_hat^2/r.
"""
import numpy as np
import pytest
from scipy.stats import norm

from predictr import Analysis

RNG = np.random.default_rng(7)
LOGNORMAL_DATA = list(RNG.lognormal(2.0, 0.5, 200))
LOGNORMAL_DF_CENSORED = sorted(RNG.lognormal(2.0, 0.5, 15))
LOGNORMAL_DS_CENSORED = sorted(RNG.lognormal(2.0, 0.5, 5) + 3)
EXP_DATA = list(RNG.exponential(scale=7.5, size=200))
EXP_DF_CENSORED = list(RNG.exponential(scale=7.5, size=15))
EXP_DS_CENSORED = list(RNG.exponential(scale=7.5, size=5))


def _inv_mills(z):
    return norm.pdf(z) / (1 - norm.cdf(z))


def _score(mu, sigma, df, ds):
    resid = np.asarray(df) - mu
    g_mu = np.sum(resid) / sigma ** 2
    g_sigma = -len(df) / sigma + np.sum(resid ** 2) / sigma ** 3
    if ds is not None and len(ds) > 0:
        z = (np.asarray(ds) - mu) / sigma
        h = _inv_mills(z)
        g_mu += np.sum(h) / sigma
        g_sigma += np.sum(z * h) / sigma
    return np.array([g_mu, g_sigma])


# --- LogNormal ---------------------------------------------------------

def test_lognormal_mle_matches_normal_on_log_data():
    a = Analysis(df=LOGNORMAL_DATA, dist='lognormal', show=False)
    a.mle()
    log_data = np.log(LOGNORMAL_DATA)
    assert a.mu == pytest.approx(np.mean(log_data))
    assert a.sigma == pytest.approx(np.std(log_data))


def test_lognormal_censored_score_vanishes_in_log_space():
    a = Analysis(df=list(LOGNORMAL_DF_CENSORED), ds=list(LOGNORMAL_DS_CENSORED),
                 dist='lognormal', show=False)
    a.mle()
    g = _score(a.mu, a.sigma, np.log(LOGNORMAL_DF_CENSORED), np.log(LOGNORMAL_DS_CENSORED))
    np.testing.assert_allclose(g, [0.0, 0.0], atol=1e-6)


def test_lognormal_requires_positive_data():
    with pytest.raises(ValueError):
        Analysis(df=[1.0, -2.0, 3.0], dist='lognormal', show=False).mle()


def test_lognormal_fisher_bounds_uncensored_matches_textbook_variance_in_log_space():
    a = Analysis(df=LOGNORMAL_DATA, dist='lognormal', bounds='fb', show=False)
    a.mle()
    n = len(LOGNORMAL_DATA)
    assert a.se_mu == pytest.approx(a.sigma / np.sqrt(n))
    assert a.se_sigma == pytest.approx(a.sigma / np.sqrt(2 * n))


@pytest.mark.parametrize('bounds_type', ['2s', '1su', '1sl'])
def test_lognormal_fisher_bounds_direction_and_positivity(bounds_type):
    a = Analysis(df=LOGNORMAL_DATA, dist='lognormal', bounds='fb',
                 bounds_type=bounds_type, cl=0.9, show=False)
    a.mle()
    t_p = np.exp(a.mu + a.sigma * norm.ppf(a.unrel))
    if bounds_type == '2s':
        assert np.all(a.bounds_lower < t_p) and np.all(t_p < a.bounds_upper)
        assert np.all(a.bounds_lower > 0)
    elif bounds_type == '1su':
        assert a.bounds_lower is None
        assert np.all(t_p < a.bounds_upper)
    elif bounds_type == '1sl':
        assert a.bounds_upper is None
        assert np.all(a.bounds_lower < t_p)
        assert np.all(a.bounds_lower > 0)


def test_lognormal_plot_runs_without_error():
    import matplotlib
    matplotlib.use('Agg')
    a = Analysis(df=LOGNORMAL_DATA, dist='lognormal', bounds='fb', show=False)
    a.mle()
    a.plot()


def test_lognormal_plot_uses_weibull_style_log_x_paper():
    """LogNormal/Exponential plot on the same log-x, double-log-y paper as
    Weibull (not the alternative linear-x papers considered and rejected -
    see conversation), and the x-axis must be pinned to a sensible,
    data-driven range rather than autoscaling to the MLE line's own wide
    evaluation extremes (evaluated near p=0.0001, which on a log axis can
    extend towards 0 far past anything meaningful)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    a = Analysis(df=LOGNORMAL_DATA, dist='lognormal', bounds='fb', show=False)
    a.mle()
    a.plot()
    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    left, right = ax.get_xlim()
    # The actual data spans roughly exp(mu +/- 3*sigma); the axis should be
    # in the same ballpark, not many orders of magnitude wider.
    data_min, data_max = min(LOGNORMAL_DATA), max(LOGNORMAL_DATA)
    assert left > data_min / 100
    assert right < data_max * 100


# --- Exponential ---------------------------------------------------------

def test_exponential_mle_uncensored_matches_sample_mean():
    a = Analysis(df=EXP_DATA, dist='exponential', show=False)
    a.mle()
    assert a.theta == pytest.approx(np.mean(EXP_DATA))


def test_exponential_mle_censored_matches_total_time_on_test():
    a = Analysis(df=list(EXP_DF_CENSORED), ds=list(EXP_DS_CENSORED),
                 dist='exponential', show=False)
    a.mle()
    expected = (sum(EXP_DF_CENSORED) + sum(EXP_DS_CENSORED)) / len(EXP_DF_CENSORED)
    assert a.theta == pytest.approx(expected)


def test_exponential_fisher_bounds_matches_closed_form_variance():
    a = Analysis(df=EXP_DATA, dist='exponential', bounds='fb', show=False)
    a.mle()
    expected_se = a.theta / np.sqrt(len(EXP_DATA))
    assert a.se_theta == pytest.approx(expected_se)


def test_exponential_chi2_bounds_match_nelson_1982_worked_example():
    """Cross-check bounds='chi2' against Nelson, "Applied Life Data
    Analysis" (1982), p.255: theta=6097.3, t=60, 90% one-sided lower
    confidence interval, r=10 failures, no censoring -> published lower
    confidence limit for R(60) is 0.98612 (see also NCSS PASS Chapter 408,
    Example 2, which reproduces the same figure via the exact chi-square
    pivot 2*r*theta_hat/theta ~ chi2(df=2r)). bounds='fb' (asymptotic) is
    NOT expected to reproduce this exactly - only bounds='chi2' is the
    exact method Nelson's example validates."""
    theta_hat = 6097.3
    r = 10
    synthetic_df = [theta_hat] * r  # sum(df)/r == theta_hat exactly
    a = Analysis(df=synthetic_df, dist='exponential', bounds='chi2',
                 bounds_type='1sl', cl=0.9, show=False)
    a.mle()
    implied_theta_lower = a.bounds_lower[0] / (-np.log(1 - a.unrel[0]))
    r_60_lower = np.exp(-60 / implied_theta_lower)
    assert r_60_lower == pytest.approx(0.98612, abs=1e-5)


def test_exponential_bounds_fb_and_chi2_differ():
    """bounds='fb' (asymptotic Fisher-information/delta-method) and
    bounds='chi2' (exact chi-square pivot) are genuinely different methods
    for dist='exponential' - see fisher_bounds()'s docstring - so they
    must NOT produce the same bounds, even though both are valid ways to
    get a confidence interval around the same theta_hat."""
    a_fb = Analysis(df=EXP_DATA, dist='exponential', bounds='fb', show=False)
    a_fb.mle()
    a_chi2 = Analysis(df=EXP_DATA, dist='exponential', bounds='chi2', show=False)
    a_chi2.mle()
    assert not np.allclose(a_fb.bounds_lower, a_chi2.bounds_lower)
    assert not np.allclose(a_fb.bounds_upper, a_chi2.bounds_upper)


@pytest.mark.parametrize('dist', ['weibull', 'normal', 'lognormal'])
def test_bounds_chi2_rejected_for_non_exponential_dists(dist):
    with pytest.raises(ValueError, match='chi2'):
        Analysis(df=EXP_DATA, dist=dist, bounds='chi2', show=False).mle()


@pytest.mark.parametrize('bounds_type', ['2s', '1su', '1sl'])
def test_exponential_fisher_bounds_direction_and_positivity(bounds_type):
    a = Analysis(df=EXP_DATA, dist='exponential', bounds='fb',
                 bounds_type=bounds_type, cl=0.9, show=False)
    a.mle()
    t_p = a.theta * (-np.log(1 - np.array(a.unrel)))
    if bounds_type == '2s':
        assert np.all(a.bounds_lower < t_p) and np.all(t_p < a.bounds_upper)
        assert np.all(a.bounds_lower > 0)
    elif bounds_type == '1su':
        assert a.bounds_lower is None
        assert np.all(t_p < a.bounds_upper)
    elif bounds_type == '1sl':
        assert a.bounds_upper is None
        assert np.all(a.bounds_lower < t_p)
        assert np.all(a.bounds_lower > 0)


def test_exponential_plot_runs_without_error():
    import matplotlib
    matplotlib.use('Agg')
    a = Analysis(df=EXP_DATA, dist='exponential', bounds='fb', show=False)
    a.mle()
    a.plot()


def test_exponential_plot_uses_weibull_style_log_x_paper():
    """See test_lognormal_plot_uses_weibull_style_log_x_paper() - same
    reasoning applies to Exponential's plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    a = Analysis(df=EXP_DATA, dist='exponential', bounds='fb', show=False)
    a.mle()
    a.plot()
    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    left, right = ax.get_xlim()
    data_min, data_max = min(EXP_DATA), max(EXP_DATA)
    assert left > data_min / 100
    assert right < data_max * 100


def test_exponential_plot_xlim_matches_bounds_extent_when_bounds_widen_it():
    """When the Fisher bounds extend past the point-estimate curve's own
    [0.1%, 99.9%] range (as they legitimately can for small n or a poor
    exponential fit), the x-axis must grow to match - not stay clipped to
    the narrower point-estimate range."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    small_n_data = [0.30481336314657737, 0.5793918872111126, 0.633217732127894,
                     0.7576700925659532, 0.8394342818048925, 0.9118100898948334,
                     1.0110147142055477, 1.0180126386295232, 1.3201853093496474,
                     1.492172669340363]
    a = Analysis(df=small_n_data, dist='exponential', bounds='fb',
                 bounds_type='2s', show=False)
    a.mle()
    a.plot()
    left, right = plt.gca().get_xlim()
    expected_left = 10 ** (np.ceil(np.log10(min(a.bounds_lower))) - 1)
    expected_right = 10 ** (np.ceil(np.log10(max(a.bounds_upper))))
    assert left == pytest.approx(expected_left)
    assert right == pytest.approx(expected_right)


# --- Shared config guards -------------------------------------------------

@pytest.mark.parametrize('dist', ['lognormal', 'exponential'])
def test_unsupported_bounds_method_raises(dist):
    data = LOGNORMAL_DATA if dist == 'lognormal' else EXP_DATA
    with pytest.raises(ValueError):
        Analysis(df=data, dist=dist, bounds='lrb', show=False).mle()


@pytest.mark.parametrize('dist', ['lognormal', 'exponential'])
def test_unsupported_bias_correction_raises(dist):
    data = LOGNORMAL_DATA if dist == 'lognormal' else EXP_DATA
    with pytest.raises(ValueError):
        Analysis(df=data, dist=dist, bcm='c4', show=False).mle()


@pytest.mark.parametrize('dist', ['normal', 'lognormal', 'exponential'])
def test_mrr_rejects_non_weibull_dist(dist):
    data = LOGNORMAL_DATA if dist == 'lognormal' else EXP_DATA if dist == 'exponential' else EXP_DATA
    with pytest.raises(ValueError):
        Analysis(df=data, dist=dist, show=False).mrr()


def test_plotall_mult_weibull_rejects_non_weibull_objects():
    import matplotlib
    matplotlib.use('Agg')
    from predictr import PlotAll

    failures = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
    w = Analysis(df=failures, bounds='fb', show=False)
    w.mle()
    n = Analysis(df=LOGNORMAL_DATA, dist='normal', bounds='fb', show=False)
    n.mle()
    with pytest.raises(ValueError):
        PlotAll({'w': w, 'n': n}).mult_weibull()
