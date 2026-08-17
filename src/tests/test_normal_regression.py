"""
Tests for the Normal distribution support (dist='normal').

Unlike the Weibull characterization tests, these validate actual
correctness, not just "unchanged from before" - there's no prior behavior
to pin down. Key idea: since the MLE equations for Normal are hand-derived
(no scipy.stats.fit), the gradient/Hessian used for the censored Newton
solve and for the Fisher bounds are checked against finite differences of
the log-likelihood itself, and the uncensored closed-form path is checked
against its own known textbook variance formulas
(Var(mu_hat) = sigma^2/n, Var(sigma_hat) = sigma^2/(2n), Cov = 0).
"""
import numpy as np
import pytest
from scipy.stats import norm

from predictr import Analysis

RNG = np.random.default_rng(42)
DATA_UNCENSORED = list(RNG.normal(loc=5.0, scale=2.0, size=200))
DF_CENSORED = sorted(RNG.normal(5.0, 2.0, 15))
DS_CENSORED = sorted(RNG.normal(5.0, 2.0, 5) + 3)


def _inv_mills(z):
    return norm.pdf(z) / (1 - norm.cdf(z))


def _score(mu, sigma, df, ds):
    resid = np.asarray(df) - mu
    g_mu = np.sum(resid) / sigma ** 2
    g_sigma = -len(df) / sigma + np.sum(resid ** 2) / sigma ** 3
    if ds:
        z = (np.asarray(ds) - mu) / sigma
        h = _inv_mills(z)
        g_mu += np.sum(h) / sigma
        g_sigma += np.sum(z * h) / sigma
    return np.array([g_mu, g_sigma])


def _loglik(mu, sigma, df, ds):
    ll = np.sum(norm.logpdf(df, mu, sigma))
    if ds:
        ll += np.sum(np.log(1 - norm.cdf(ds, mu, sigma)))
    return ll


def test_mle_uncensored_matches_closed_form():
    a = Analysis(df=DATA_UNCENSORED, dist='normal', show=False)
    a.mle()
    assert a.mu == pytest.approx(np.mean(DATA_UNCENSORED))
    assert a.sigma == pytest.approx(np.std(DATA_UNCENSORED))


def test_mle_censored_score_vanishes_at_solution():
    """The censored MLE is found by root-finding on the hand-derived score
    equations - confirm the solver actually landed on a root."""
    a = Analysis(df=list(DF_CENSORED), ds=list(DS_CENSORED), dist='normal', show=False)
    a.mle()
    g = _score(a.mu, a.sigma, DF_CENSORED, DS_CENSORED)
    np.testing.assert_allclose(g, [0.0, 0.0], atol=1e-6)


def test_score_matches_finite_difference_of_loglikelihood():
    """Validates the hand-derived score (gradient) equations themselves,
    independent of the solver, by comparing to a numerical gradient of the
    actual censored log-likelihood."""
    mu0, sigma0, eps = 4.7, 2.3, 1e-6
    analytic = _score(mu0, sigma0, DF_CENSORED, DS_CENSORED)
    numeric = np.array([
        (_loglik(mu0 + eps, sigma0, DF_CENSORED, DS_CENSORED)
         - _loglik(mu0 - eps, sigma0, DF_CENSORED, DS_CENSORED)) / (2 * eps),
        (_loglik(mu0, sigma0 + eps, DF_CENSORED, DS_CENSORED)
         - _loglik(mu0, sigma0 - eps, DF_CENSORED, DS_CENSORED)) / (2 * eps),
    ])
    np.testing.assert_allclose(analytic, numeric, atol=1e-3)


def test_fisher_bounds_uncensored_matches_textbook_variance():
    """For uncensored Normal data, Var(mu_hat)=sigma^2/n and
    Var(sigma_hat)=sigma^2/(2n) have known closed forms independent of the
    general (and more complex) observed-information code path - this
    catches sign/indexing mistakes in fisher_bounds()'s matrix handling."""
    a = Analysis(df=DATA_UNCENSORED, dist='normal', bounds='fb', show=False)
    a.mle()
    n = len(DATA_UNCENSORED)
    assert a.se_mu == pytest.approx(a.sigma / np.sqrt(n))
    assert a.se_sigma == pytest.approx(a.sigma / np.sqrt(2 * n))
    assert a.f_inv.item(1) == pytest.approx(0.0, abs=1e-9)  # Cov(mu, sigma) = 0


@pytest.mark.parametrize('bounds_type', ['2s', '1su', '1sl'])
def test_fisher_bounds_direction(bounds_type):
    a = Analysis(df=DATA_UNCENSORED, dist='normal', bounds='fb',
                 bounds_type=bounds_type, cl=0.9, show=False)
    a.mle()
    t_p = a.mu + a.sigma * norm.ppf(a.unrel)
    if bounds_type == '2s':
        assert np.all(a.bounds_lower < t_p)
        assert np.all(t_p < a.bounds_upper)
    elif bounds_type == '1su':
        assert a.bounds_lower is None
        assert np.all(t_p < a.bounds_upper)
    elif bounds_type == '1sl':
        assert a.bounds_upper is None
        assert np.all(a.bounds_lower < t_p)


def test_unsupported_bounds_method_raises():
    with pytest.raises(ValueError):
        Analysis(df=DATA_UNCENSORED, dist='normal', bounds='lrb', show=False).mle()


def test_unsupported_bias_correction_raises():
    with pytest.raises(ValueError):
        Analysis(df=DATA_UNCENSORED, dist='normal', bcm='c4', show=False).mle()


def test_invalid_dist_raises():
    with pytest.raises(ValueError):
        Analysis(df=DATA_UNCENSORED, dist='gamma', show=False)


def test_default_plot_title_is_dist_aware():
    a = Analysis(df=DATA_UNCENSORED, dist='normal', show=False)
    assert a.plot_title == 'Normal Probability Plot'
    w = Analysis(df=DATA_UNCENSORED, show=False)
    assert w.plot_title == 'Weibull Probability Plot'


def test_plot_runs_without_error():
    import matplotlib
    matplotlib.use('Agg')
    a = Analysis(df=DATA_UNCENSORED, dist='normal', bounds='fb', show=False)
    a.mle()
    a.plot()  # must not raise
