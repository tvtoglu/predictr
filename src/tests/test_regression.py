"""
Characterization tests for the lifetime-regression class (Regression):
Weibull AFT and Cox PH, uncensored/right-censored, Wald ('fb') and
profile-likelihood ('lrb') bounds, the two data-input paths, prediction
helpers and the plots.

The pinned numbers are this implementation's own output on the fixed
dataset below (21 failures / 32 units, two tied event times at 28.4).
They exist to catch accidental behaviour changes, not to validate against
an external reference; the analytic score/information are checked against
finite differences in a separate step during development.
"""
import numpy as np
import pandas as pd
import pytest

from predictr import Regression

# --- fixed dataset -------------------------------------------------------
TIME = [72.7, 28.0, 28.4, 37.4, 8.3, 17.6, 27.9, 16.4, 86.5, 19.8, 11.2,
        30.7, 13.9, 8.5, 27.2, 23.8, 19.0, 48.4, 27.7, 22.8, 7.4, 14.6,
        38.8, 16.6, 22.2, 43.1, 21.3, 33.0, 19.6, 28.4, 32.3, 24.7]
EVENT = [0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0,
         1, 1, 0, 1, 0, 1, 1, 1, 1, 1]
TEMP = [60.0, 100.0, 80.0, 80.0, 80.0, 100.0, 60.0, 100.0, 60.0, 60.0, 80.0,
        100.0, 100.0, 100.0, 100.0, 100.0, 80.0, 60.0, 100.0, 80.0, 80.0,
        80.0, 60.0, 100.0, 100.0, 80.0, 80.0, 100.0, 80.0, 80.0, 80.0, 60.0]
LOAD = [1.0, 1.5, 2.0, 1.0, 2.0, 2.0, 1.0, 1.5, 1.0, 2.0, 2.0, 1.5, 1.0,
        2.0, 1.5, 2.0, 2.0, 2.0, 2.0, 1.0, 1.5, 1.5, 1.5, 1.0, 1.5, 1.0,
        2.0, 2.0, 2.0, 2.0, 1.5, 2.0]


@pytest.fixture
def data():
    return pd.DataFrame({'time': TIME, 'event': EVENT,
                         'temp': TEMP, 'load': LOAD})


def _fit(data, **kw):
    kw.setdefault('duration_col', 'time')
    kw.setdefault('event_col', 'event')
    kw.setdefault('covariate_cols', ['temp', 'load'])
    return Regression(data=data, **kw).fit()


# --- Weibull AFT -------------------------------------------------------------
def test_aft_point_estimates(data):
    r = _fit(data, model='weibull_aft', bounds='fb')
    assert r.coef == pytest.approx([-0.0177537, -0.452861], rel=1e-5)
    assert r.se_coef == pytest.approx([0.0047125, 0.1950878], rel=1e-5)
    assert r.intercept == pytest.approx(5.7668593, rel=1e-5)
    assert r.sigma == pytest.approx(0.314438, rel=1e-5)
    assert r.beta == pytest.approx(3.1802772, rel=1e-5)
    assert r.loglik == pytest.approx(-81.6068629, rel=1e-6)
    assert r.aic == pytest.approx(171.2137258, rel=1e-6)
    assert r.lr_stat == pytest.approx(23.6738479, rel=1e-5)
    assert r.n == 32 and r.n_events == 21


def test_aft_wald_bounds(data):
    r = _fit(data, model='weibull_aft', bounds='fb')
    assert r.ci_lower == pytest.approx([-0.025505, -0.7737519], rel=1e-5)
    assert r.ci_upper == pytest.approx([-0.0100023, -0.1319702], rel=1e-5)


def test_aft_profile_bounds_differ_from_wald_same_point(data):
    r = _fit(data, model='weibull_aft', bounds='lrb')
    assert r.coef == pytest.approx([-0.0177537, -0.452861], rel=1e-5)
    assert r.ci_lower == pytest.approx([-0.0263422, -0.7932511], rel=1e-4)
    assert r.ci_upper == pytest.approx([-0.0102558, -0.1263934], rel=1e-4)


def test_aft_predictions(data):
    r = _fit(data, model='weibull_aft', bounds=None)
    prof = pd.DataFrame({'temp': [60.0, 100.0], 'load': [1.0, 2.0]})
    assert r.predict_median(prof) == pytest.approx([62.3979387, 19.5021414],
                                                  rel=1e-5)
    assert r.predict_time_ratio(prof) == pytest.approx([1.9867813, 0.6209578],
                                                       rel=1e-5)
    s = r.predict_survival(prof, times=[40.0]).iloc[0].to_numpy()
    assert s == pytest.approx([0.8449046, 0.0011052], rel=1e-4)


# --- Cox PH ---------------------------------------------------------------
def test_cox_efron_point_estimates(data):
    r = _fit(data, model='cox_ph', ties='efron', bounds='fb')
    assert r.coef == pytest.approx([0.0452352, 0.990813], rel=1e-5)
    assert r.se_coef == pytest.approx([0.0172534, 0.6328163], rel=1e-5)
    assert r.loglik == pytest.approx(-46.5896972, rel=1e-6)
    assert r.aic == pytest.approx(97.1793944, rel=1e-6)
    assert r.lr_stat == pytest.approx(12.8527903, rel=1e-5)
    assert r.concordance == pytest.approx(0.7058824, rel=1e-5)
    assert r.hazard_ratio == pytest.approx(np.exp([0.0452352, 0.990813]),
                                           rel=1e-5)


def test_cox_baseline_cumhaz(data):
    r = _fit(data, model='cox_ph', ties='efron', bounds=None)
    times, H0 = r.baseline_cumhaz
    assert times[0] == pytest.approx(8.5) and times[-1] == pytest.approx(86.5)
    assert len(times) == 20
    assert H0[0] == pytest.approx(0.0244709, rel=1e-5)
    assert H0[-1] == pytest.approx(7.7445152, rel=1e-5)
    assert np.all(np.diff(H0) > 0)


def test_cox_efron_vs_breslow_differ(data):
    e = _fit(data, model='cox_ph', ties='efron', bounds='fb')
    b = _fit(data, model='cox_ph', ties='breslow', bounds='fb')
    assert b.coef == pytest.approx([0.0455926, 0.9744869], rel=1e-5)
    assert not np.allclose(e.coef, b.coef)


def test_cox_profile_bounds(data):
    r = _fit(data, model='cox_ph', bounds='lrb')
    assert r.ci_lower == pytest.approx([0.0183071, -0.0307886], rel=1e-4)
    assert r.ci_upper == pytest.approx([0.0756241, 2.0740074], rel=1e-4)


def test_cox_predictions(data):
    r = _fit(data, model='cox_ph', bounds=None)
    prof = pd.DataFrame({'temp': [60.0, 100.0], 'load': [1.0, 2.0]})
    assert r.predict_hazard_ratio(prof) == pytest.approx([0.1920805, 3.1593948],
                                                         rel=1e-5)
    S = r.predict_survival(prof, times=r.baseline_cumhaz[0])
    assert list(S.index[:1]) == [8.5]
    assert S.iloc[0, 0] == pytest.approx(0.9953107, rel=1e-5)
    assert S.iloc[-1, 1] == pytest.approx(0.0, abs=1e-6)
    # default grid runs back to (almost) t = 0, where S = 1
    Sg = r.predict_survival(prof)
    assert Sg.index[0] < 0.01 * max(TIME)
    assert np.all(Sg.iloc[0].to_numpy() > 0.999)
    # survival is non-increasing in time and in [0, 1]
    assert np.all(np.diff(Sg.to_numpy(), axis=0) <= 1e-12)
    assert Sg.to_numpy().min() >= -1e-12 and Sg.to_numpy().max() <= 1 + 1e-12


# --- data-input paths ---------------------------------------------------
def test_df_ds_path_matches_frame_path(data):
    frame = _fit(data, model='cox_ph', bounds='fb')
    m = np.array(EVENT) == 1
    t = np.array(TIME)
    cov = data[['temp', 'load']]
    split = Regression(df=list(t[m]), ds=list(t[~m]),
                       x_df=cov[m], x_ds=cov[~m],
                       model='cox_ph', bounds='fb').fit()
    assert split.coef == pytest.approx(frame.coef, rel=1e-8)
    assert split.se_coef == pytest.approx(frame.se_coef, rel=1e-8)


def test_x_df_accepts_dict_and_array(data):
    ref = _fit(data, model='weibull_aft', bounds=None)
    m = np.array(EVENT) == 1
    t = np.array(TIME)
    tp, ld = np.array(TEMP), np.array(LOAD)
    as_dict = Regression(
        df=list(t[m]), ds=list(t[~m]),
        x_df={'temp': list(tp[m]), 'load': list(ld[m])},
        x_ds={'temp': list(tp[~m]), 'load': list(ld[~m])},
        model='weibull_aft', bounds=None).fit()
    as_arr = Regression(
        df=list(t[m]), ds=list(t[~m]),
        x_df=np.c_[tp[m], ld[m]], x_ds=np.c_[tp[~m], ld[~m]],
        feature_names=['temp', 'load'], model='weibull_aft', bounds=None).fit()
    assert as_dict.coef == pytest.approx(ref.coef, rel=1e-8)
    assert as_arr.coef == pytest.approx(ref.coef, rel=1e-8)


def test_string_covariate_is_one_hot_encoded(data):
    d = data.copy()
    d['material'] = np.where(np.arange(len(d)) % 2 == 0, 'A', 'B')
    r = Regression(data=d, duration_col='time', event_col='event',
                   covariate_cols=['temp', 'material'],
                   model='cox_ph', bounds=None).fit()
    assert 'material_B' in r.feature_names
    assert len(r.coef) == len(r.feature_names)


# --- summary / plots --------------------------------------------------------
def test_summary_prints_and_returns_frame(data, capsys):
    r = _fit(data, model='weibull_aft', bounds='fb')
    out = r.summary()
    captured = capsys.readouterr().out
    assert isinstance(out, pd.DataFrame)
    assert 'Intercept' in out.index
    assert list(r.feature_names) == ['temp', 'load']
    assert 'Weibull AFT' in captured and 'AIC' in captured
    silent = r.summary(print_report=False)
    assert capsys.readouterr().out == ''
    assert isinstance(silent, pd.DataFrame)


def test_plots_return_figures(data):
    import matplotlib
    matplotlib.use('Agg')
    r = _fit(data, model='cox_ph', bounds='fb')
    prof = pd.DataFrame({'temp': [60.0, 100.0], 'load': [1.0, 2.0]},
                        index=['low', 'high'])
    ra = _fit(data, model='weibull_aft', bounds='fb')
    # show=False -> the Figure is returned for programmatic use
    assert r.plot(show=False) is not None
    assert r.plot_survival(prof, show=False) is not None
    assert ra.plot(show=False) is not None
    assert ra.plot_survival(prof, show=False) is not None
    # plot* methods default to show=True, which draws and returns None
    assert ra.plot() is None
    assert ra.plot_survival(prof) is None


# --- error handling -------------------------------------------------------
def test_errors(data):
    with pytest.raises(ValueError):                       # no covariates
        Regression(df=[1.0, 2.0, 3.0], model='weibull_aft')
    with pytest.raises(ValueError):                       # unknown bounds
        _fit(data, model='cox_ph', bounds='wald')
    with pytest.raises(ValueError):                       # unknown model
        _fit(data, model='logistic')
    with pytest.raises(ValueError):                       # both input paths
        Regression(df=[1.0], x_df=[[1.0]], data=data,
                   duration_col='time', event_col='event')
    with pytest.raises(ValueError):                       # reserved kwarg
        _fit(data, model='cox_ph', strata='temp')
    with pytest.raises(ValueError):                       # x_df length mismatch
        Regression(df=[1.0, 2.0], x_df=[[1.0]], model='weibull_aft')
    with pytest.raises(ValueError):                       # non-positive time
        Regression(df=[0.0, 1.0], x_df=[[1.0], [2.0]], model='weibull_aft')
    with pytest.raises(ValueError):                       # predict before fit
        Regression(data=data, duration_col='time', event_col='event',
                   covariate_cols=['temp']).predict_median(
                       pd.DataFrame({'temp': [80.0]}))


def test_cox_perfect_separation_raises():
    # covariate x orders the (all-failure) event times exactly -> the
    # partial-likelihood maximum is at beta = -inf.
    n = 10
    t = list(np.arange(1.0, n + 1))
    Regression(df=t, x_df=[[float(i)] for i in range(n)],
               model='cox_ph', bounds=None)
    with pytest.raises(ValueError, match='did not converge'):
        Regression(df=t, x_df=[[float(i)] for i in range(n)],
                   model='cox_ph', bounds=None).fit()


# --- confidence bands on the predicted survival function ----------------
_PROF = pd.DataFrame({'temp': [60.0, 100.0], 'load': [1.0, 2.0]},
                     index=['mild', 'hard'])


@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
@pytest.mark.parametrize('bnd', ['fb', 'lrb'])
def test_survival_band_shape_and_order(data, model, bnd):
    r = _fit(data, model=model, bounds=bnd)
    out = r.predict_survival(_PROF, ci=True, simultaneous=True)
    assert set(out) >= {'surv', 'lower', 'upper', 'lower_sim', 'upper_sim',
                        'in_data_range', 'method', 'cl'}
    S, L, U = out['surv'].to_numpy(), out['lower'].to_numpy(), out['upper'].to_numpy()
    assert np.all(L <= S + 1e-9) and np.all(S <= U + 1e-9)
    assert L.min() >= -1e-9 and U.max() <= 1 + 1e-9
    SL, SU = out['lower_sim'].to_numpy(), out['upper_sim'].to_numpy()
    # the simultaneous band encloses the pointwise band ...
    m = np.isfinite(SL) & np.isfinite(SU)
    assert np.all(SL[m] <= L[m] + 1e-6) and np.all(SU[m] >= U[m] - 1e-6)
    # ... and both of its edges are monotone non-increasing in time
    for j in range(SL.shape[1]):
        for edge in (SL[:, j], SU[:, j]):
            fin = edge[np.isfinite(edge)]
            assert np.all(np.diff(fin) <= 1e-9)
    assert out['cl'] == pytest.approx(0.9)


def test_cox_lrb_band_equals_fb_band(data):
    a = _fit(data, model='cox_ph', bounds='fb').predict_survival(_PROF, ci=True)
    b = _fit(data, model='cox_ph', bounds='lrb').predict_survival(_PROF, ci=True)
    assert b['method'] == 'lrb'
    np.testing.assert_allclose(a['lower'].to_numpy(), b['lower'].to_numpy())
    np.testing.assert_allclose(a['upper'].to_numpy(), b['upper'].to_numpy())


def test_band_cl_widens(data):
    r = _fit(data, model='weibull_aft', bounds='fb')
    narrow = r.predict_survival(_PROF, ci=True, cl=0.80)
    wide = r.predict_survival(_PROF, ci=True, cl=0.99)
    w80 = (narrow['upper'] - narrow['lower']).to_numpy()
    w99 = (wide['upper'] - wide['lower']).to_numpy()
    assert np.nanmean(w99) > np.nanmean(w80)


def test_predict_quantile_ci_aft(data):
    r = _fit(data, model='weibull_aft', bounds='fb')
    for method in ('fb', 'lrb'):
        tq, lo, hi = r.predict_quantile(_PROF, q=0.1, ci=True, bounds=method)
        assert np.all(lo < tq) and np.all(tq < hi)
    assert r.predict_quantile(_PROF, q=0.1).shape == (2,)   # ci=False unchanged


def test_cox_bq_life_triplet(data):
    r = _fit(data, model='cox_ph', bounds='fb')
    u, m, o = r._cox_quantile_ci(np.array([60.0, 1.0]), 0.1, 0.9)
    assert u <= m <= o
    assert np.isfinite(m)


@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
def test_plot_survival_band_options(data, model):
    import matplotlib
    matplotlib.use('Agg')
    r = _fit(data, model=model, bounds='fb')
    assert r.plot_survival(_PROF, ci=True, show=False) is not None
    assert r.plot_survival(_PROF, ci=True, simultaneous=True,
                           target_bq=0.1, show=False) is not None
    assert r.plot_survival(_PROF, target_bq=0.05, show=False) is not None


def test_survival_band_rejects_bad_bounds(data):
    r = _fit(data, model='weibull_aft', bounds='fb')
    with pytest.raises(ValueError):
        r.predict_survival(_PROF, ci=True, bounds='wald')


@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
@pytest.mark.parametrize('bnd', ['npbb', 'pbb'])
def test_bootstrap_bounds(data, model, bnd):
    r = _fit(data, model=model, bounds=bnd, n_boot=80, cl=0.9)
    s = r.summary(print_report=False)
    lo_col = [c for c in s.columns if c.startswith('coef lower')][0]
    hi_col = [c for c in s.columns if c.startswith('coef upper')][0]
    assert np.all(np.isfinite(s[lo_col].to_numpy()))
    assert np.all(s[lo_col].to_numpy() <= s[hi_col].to_numpy() + 1e-9)

    band = r.predict_survival(_PROF, ci=True)
    lo = band['lower'].to_numpy()
    hi = band['upper'].to_numpy()
    assert lo.shape == band['surv'].shape
    assert np.all(lo <= hi + 1e-9)
    assert band['lower_sim'] is None          # no simultaneous bootstrap band
    assert r.plot_survival(_PROF, ci=True, target_bq=0.1, show=False) is not None
    if model == 'weibull_aft':
        tq, qlo, qhi = r.predict_quantile(_PROF, q=0.1, ci=True)
        assert np.all(qlo <= tq + 1e-9) and np.all(tq <= qhi + 1e-9)


def test_power_analysis_and_sample_size(data):
    r = _fit(data, model='weibull_aft', bounds=None, cl=0.9)
    pw = r.power_analysis(n_sim=40, seed=0, print_report=False)
    assert list(pw.index[:-1]) == list(r.feature_names)
    assert pw.index[-1] == '(model)'
    assert np.all((pw.to_numpy() >= 0.0) & (pw.to_numpy() <= 1.0))
    # a large hypothetical effect is detected almost always
    strong = np.array(r.coef, dtype=float)
    strong[:] = [v if v != 0 else 0.0 for v in strong]
    hi = r.power_analysis(n=200, n_sim=40, coef=strong * 3.0 - 1.0,
                          seed=0, print_report=False)
    assert hi['(model)'] >= pw['(model)'] - 1e-9
    ss = r.sample_size(target_power=0.5, n_grid=[15, 40], n_sim=25,
                       print_report=False)
    assert ss is None or ss in (15, 40)


# --- goodness of fit ---------------------------------------------------
@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
def test_residuals_and_concordance(data, model):
    r = _fit(data, model=model, bounds='fb')
    assert 0.0 <= r.concordance <= 1.0        # now set for AFT too
    res = r.residuals()
    assert list(res.columns) == ['cox_snell', 'martingale', 'deviance']
    assert len(res) == r.n
    assert np.all(res['cox_snell'].to_numpy() >= 0.0)
    assert np.all(res['martingale'].to_numpy() <= 1.0 + 1e-9)
    # martingale = event - cox_snell, by definition
    np.testing.assert_allclose(res['martingale'].to_numpy(),
                               r._event - res['cox_snell'].to_numpy())
    assert np.all(np.isfinite(res['deviance'].to_numpy()))
    with pytest.raises(ValueError):
        r.residuals('nonsense')


@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
def test_goodness_of_fit(data, model):
    r = _fit(data, model=model, bounds='fb')
    s = r.goodness_of_fit(print_report=False)
    for key in ('concordance', 'aic', 'lr_pvalue', 'cox_snell_slope',
                'cox_snell_max_dev', 'absolute_fit', 'discrimination',
                'verdict'):
        assert key in s.index
    assert np.isfinite(s['cox_snell_slope']) and s['cox_snell_slope'] > 0
    assert s['absolute_fit'] in ('good', 'marginal', 'poor')
    assert s['discrimination'] in ('weak', 'modest', 'good', 'strong')
    assert s['verdict'][:4] in ('GOOD', 'MARG', 'POOR')
    if model == 'cox_ph':
        assert s['proportional_hazards'] in ('good', 'marginal', 'poor')
    assert r.gof.__func__ is r.goodness_of_fit.__func__
    import matplotlib
    matplotlib.use('Agg')
    assert r.plot_gof(show=False) is not None
    assert r.plot_survival(_PROF, km_overlay=True, show=False) is not None


@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
def test_kaplan_meier_nelson_aalen(data, model):
    import matplotlib
    matplotlib.use('Agg')
    r = Regression(data=data, duration_col='time', event_col='event',
                   covariate_cols=['temp', 'load'], model=model)  # no fit()
    km = r.kaplan_meier()
    assert list(km.columns) == ['time', 'n_risk', 'n_event', 'n_censor',
                                'surv', 'surv_se', 'surv_lower', 'surv_upper']
    s = km['surv'].to_numpy()
    assert np.all(np.diff(s) <= 1e-12)                       # non-increasing
    assert s[0] <= 1.0 and s[-1] >= 0.0
    assert np.all(km['surv_lower'] <= km['surv'] + 1e-9)
    assert np.all(km['surv_upper'] >= km['surv'] - 1e-9)
    assert km['n_risk'].iloc[0] == r.n

    na = r.nelson_aalen()
    h = na['cumhaz'].to_numpy()
    assert np.all(np.diff(h) >= -1e-12)                      # non-decreasing

    # KM survival and NA hazard are consistent: S ~ exp(-H) early on
    m = km['surv'] > 0.2
    assert np.allclose(km['surv'][m], np.exp(-na['cumhaz'])[m], atol=0.05)

    # grouping by a covariate column and by an explicit label array
    g = r.kaplan_meier(by='load')
    assert set(g['group']) == {f'load={v}' for v in data['load'].unique()}
    lab = np.where(data['temp'].to_numpy() < 90, 'lo', 'hi')
    assert set(r.nelson_aalen(by=lab)['group']) == {'group=lo', 'group=hi'}

    assert r.plot_km(by='load', show=False) is not None
    assert r.plot_na(ci=False, show=False) is not None
    assert r.plot_km() is None                               # default show=True
    with pytest.raises(ValueError):
        r.kaplan_meier(by='not_a_column')
    with pytest.raises(ValueError):
        r.nelson_aalen(by=[0, 1, 2])                         # wrong length


def test_goodness_of_fit_flags_a_bad_fit():
    # heavy-tailed times with a useless covariate: the Weibull link and the
    # covariate should both look bad
    rng = np.random.default_rng(0)
    n = 80
    t = rng.pareto(0.7, n) + 0.05
    dfb = pd.DataFrame({'time': t, 'event': np.ones(n, int),
                        'x': rng.normal(size=n)})
    r = Regression(data=dfb, duration_col='time', event_col='event',
                   covariate_cols=['x'], model='weibull_aft').fit()
    s = r.goodness_of_fit(print_report=False)
    assert s['absolute_fit'] == 'poor'
    assert s['verdict'].startswith('POOR FIT')


def test_check_ph_cox_and_rejects_aft(data):
    r = _fit(data, model='cox_ph', bounds='fb')
    tab = r.check_ph(print_report=False)
    assert list(tab.index) == list(r.feature_names) + ['GLOBAL']
    p = tab['p'].to_numpy()
    assert np.all((p[np.isfinite(p)] >= 0.0) & (p[np.isfinite(p)] <= 1.0))
    for tr in ('rank', 'log', 'identity'):
        assert 'GLOBAL' in r.check_ph(transform=tr, print_report=False).index
    with pytest.raises(ValueError):
        r.check_ph(transform='bogus')
    with pytest.raises(ValueError):
        _fit(data, model='weibull_aft', bounds='fb').check_ph()


# --- named life-stress (aging) laws ----------------------------------
def _alt_data(seed=7):
    """Small accelerated-test data set: 3 temperatures x 2 voltages,
    generated from Arrhenius (Ea=0.7 eV) + inverse power (n=2.5), Weibull
    shape 2, with random right-censoring."""
    kB = 8.617333262e-5
    rng = np.random.default_rng(seed)
    Ea, nexp, beta = 0.7, 2.5, 2.0
    C = 2000.0 / (np.exp(Ea / (kB * (125 + 273.15))) * 5.0 ** (-nexp))
    rows = []
    for TC in (85.0, 105.0, 125.0):
        for V in (3.0, 5.0):
            eta = C * np.exp(Ea / (kB * (TC + 273.15))) * V ** (-nexp)
            for _ in range(30):
                t = eta * (-np.log(rng.uniform())) ** (1 / beta)
                c = rng.uniform(1500, 25000)
                rows.append((TC, V, round(min(t, c), 1), int(t <= c)))
    return pd.DataFrame(rows, columns=['temp_C', 'volt', 'hours', 'failed'])


def test_stress_model_recovers_physical_parameters():
    df = _alt_data()
    r = Regression(data=df, duration_col='hours', event_col='failed',
                   model='weibull_aft', bounds='lrb', cl=0.9,
                   stress_model={'temp_C': 'arrhenius',
                                 'volt': 'inverse_power'}).fit()
    assert r.feature_names == ['temp_C_invkT', 'volt_ln']
    sp = r.stress_params.set_index('parameter')
    assert sp.loc['Ea_eV', 'ci_lower'] < 0.7 < sp.loc['Ea_eV', 'ci_upper']
    assert sp.loc['n', 'ci_lower'] < 2.5 < sp.loc['n', 'ci_upper']
    assert 1.5 < r.beta < 2.6
    # summary prints the life-stress block and still returns the table
    tbl = r.summary(print_report=True)
    assert 'coef' in tbl.columns


def test_stress_model_predict_and_af_in_raw_units():
    df = _alt_data()
    r = Regression(data=df, duration_col='hours', event_col='failed',
                   model='weibull_aft', bounds='fb',
                   stress_model={'temp_C': 'arrhenius',
                                 'volt': 'inverse_power'}).fit()
    use = {'temp_C': 55, 'volt': 3.3}
    tq, lo, hi = r.predict_quantile(use, q=0.1, ci=True)
    assert lo[0] < tq[0] < hi[0]
    # milder field condition -> much longer life than the harshest test cell
    f, fl, fh = r.acceleration_factor(from_={'temp_C': 125, 'volt': 5.0},
                                      to_=use)
    assert f > 1.0 and fl < f < fh
    assert np.isfinite(r.predict_median(use)[0])
    # a DataFrame of raw conditions works too
    grid = pd.DataFrame({'temp_C': [55, 85], 'volt': [3.3, 3.3]},
                        index=['field', 'mild_test'])
    S = r.predict_survival(grid, times=np.linspace(1, 5e4, 50))
    assert list(S.columns) == ['field', 'mild_test']


def test_stress_model_diagnostics_and_plots():
    import matplotlib
    matplotlib.use('Agg')
    df = _alt_data()
    r = Regression(data=df, duration_col='hours', event_col='failed',
                   model='weibull_aft',
                   stress_model={'temp_C': 'arrhenius',
                                 'volt': 'inverse_power'}).fit()
    tab = r.check_shape(print_report=False)
    assert list(tab.index)[-1] == 'pooled'
    assert 'beta' in tab.columns
    assert r.plot_stress_life(show=False) is not None
    assert r.plot_survival({'temp_C': 55, 'volt': 3.3}, ci=True,
                           target_bq=0.1, show=False,
                           times=np.linspace(1, 2e5, 80)) is not None


def test_stress_model_errors_and_cox():
    df = _alt_data()
    with pytest.raises(ValueError):                       # unknown law
        Regression(data=df, duration_col='hours', event_col='failed',
                   stress_model={'temp_C': 'bogus'})
    with pytest.raises(ValueError):                       # column not present
        Regression(data=df, duration_col='hours', event_col='failed',
                   stress_model={'nope': 'arrhenius'})
    with pytest.raises(ValueError):                       # bad unit
        Regression(data=df, duration_col='hours', event_col='failed',
                   stress_model={'temp_C': 'arrhenius'},
                   stress_units={'temp_C': 'F'})
    # single stress level -> constant term -> clear error
    one = df[(df.temp_C == 85.0) & (df.volt == 3.0)]
    with pytest.raises(ValueError):
        Regression(data=one, duration_col='hours', event_col='failed',
                   stress_model={'temp_C': 'arrhenius'})

    # Cox: transformed columns fit, but no physical-parameter table
    rc = Regression(data=df, duration_col='hours', event_col='failed',
                    model='cox_ph',
                    stress_model={'temp_C': 'arrhenius',
                                  'volt': 'inverse_power'}).fit()
    assert rc.stress_params is None
    assert np.isfinite(rc.predict_hazard_ratio({'temp_C': 55, 'volt': 3.3})[0])
    with pytest.raises(ValueError):
        rc.acceleration_factor(from_={'temp_C': 125, 'volt': 5.0},
                               to_={'temp_C': 55, 'volt': 3.3})
    with pytest.raises(ValueError):
        rc.check_shape()


# --- bounds default / stress_model shorthand ---------------------------
@pytest.mark.parametrize('model', ['weibull_aft', 'cox_ph'])
def test_bounds_default_is_none(data, model, capsys):
    r = Regression(data=data, duration_col='time', event_col='event',
                   covariate_cols=['temp', 'load'], model=model).fit()
    assert r.bounds is None and r.bounds_type is None
    assert np.all(np.isnan(r.ci_lower)) and np.all(np.isnan(r.ci_upper))
    r.summary()
    out = capsys.readouterr().out
    assert 'bounds:' not in out
    # explicitly asking for bounds still works
    r2 = Regression(data=data, duration_col='time', event_col='event',
                    covariate_cols=['temp', 'load'], model=model,
                    bounds='fb').fit()
    assert np.all(np.isfinite(r2.ci_lower))


def test_stress_model_string_shorthand():
    df = _alt_data()
    one = df[df.volt == 5.0][['temp_C', 'hours', 'failed']]
    r_str = Regression(data=one, duration_col='hours', event_col='failed',
                       model='weibull_aft', bounds='fb',
                       stress_model='arrhenius').fit()
    r_dict = Regression(data=one, duration_col='hours', event_col='failed',
                        model='weibull_aft', bounds='fb',
                        stress_model={'temp_C': 'arrhenius'}).fit()
    assert r_str.feature_names == ['temp_C_invkT'] == r_dict.feature_names
    np.testing.assert_allclose(r_str.coef, r_dict.coef)
    # ambiguous with more than one covariate
    with pytest.raises(ValueError):
        Regression(data=df, duration_col='hours', event_col='failed',
                   stress_model='arrhenius')
