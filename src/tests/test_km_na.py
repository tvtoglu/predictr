"""
Non-parametric Kaplan-Meier / Nelson-Aalen estimators, shared by the
Analysis class (df/ds lists, no DataFrame) and the Regression class
(single life-table engine). These pin the table layout, the monotonicity
and the S ~ exp(-H) consistency, and check that show= is handled the
Analysis way (constructor flag, draws like mle()).
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import pytest

from predictr import (Analysis, Regression, PREDICTR_PALETTE,
                      PREDICTR_FIT_COLOR)

FAILURES = [0.45, 0.68, 0.77, 0.88, 1.48, 1.63, 2.10, 2.90]
SUSPENSIONS = [1.0, 1.2, 3.3]

KM_COLS = ['time', 'n_risk', 'n_event', 'n_censor',
           'surv', 'surv_se', 'surv_lower', 'surv_upper']
NA_COLS = ['time', 'n_risk', 'n_event', 'n_censor',
           'cumhaz', 'cumhaz_se', 'cumhaz_lower', 'cumhaz_upper']


def test_analysis_kaplan_meier_table():
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, show=False)
    km = a.kaplan_meier()
    assert list(km.columns) == KM_COLS
    # one row per distinct observed time
    assert len(km) == len(set(FAILURES + SUSPENSIONS))
    assert km['n_risk'].iloc[0] == len(FAILURES) + len(SUSPENSIONS)
    s = km['surv'].to_numpy()
    assert np.all(np.diff(s) <= 1e-12) and s[0] <= 1.0 and s[-1] >= 0.0
    # survival only drops at event times, flat across a pure censor time
    flat = km.loc[km['n_event'] == 0, 'surv'].to_numpy()
    assert np.all(np.isin(flat, s))
    assert np.all(km['surv_lower'] <= km['surv'] + 1e-9)
    assert np.all(km['surv_upper'] >= km['surv'] - 1e-9)


def test_analysis_nelson_aalen_and_consistency():
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, show=False)
    na = a.nelson_aalen(cl=0.95)
    assert list(na.columns) == NA_COLS
    h = na['cumhaz'].to_numpy()
    assert np.all(np.diff(h) >= -1e-12)
    km = a.kaplan_meier(cl=0.95)
    m = km['surv'].to_numpy() > 0.3
    assert np.allclose(km['surv'].to_numpy()[m], np.exp(-h)[m], atol=0.05)


def test_analysis_no_dataframe_needed_and_no_suspensions():
    a = Analysis(df=FAILURES, show=False)          # no ds at all
    km = a.kaplan_meier()
    assert km['n_censor'].sum() == 0
    assert km['surv'].iloc[-1] == pytest.approx(0.0, abs=1e-12)


def test_analysis_km_needs_a_failure():
    with pytest.raises(ValueError):
        Analysis(ds=SUSPENSIONS, show=False).kaplan_meier()


def test_analysis_show_draws_like_mle():
    # show=True -> the method draws (Agg backend: no window, must not raise)
    # and still returns the life table
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, show=True)
    km = a.kaplan_meier()
    na = a.nelson_aalen()
    assert list(km.columns) == KM_COLS
    assert list(na.columns) == NA_COLS


def test_grouped_km_uses_categorical_palette():
    rng = np.random.default_rng(1)
    n = 60
    frame = pd.DataFrame({'t': rng.uniform(1, 40, n),
                          'e': rng.integers(0, 2, n),
                          'g': rng.choice(list('ABC'), n)})
    r = Regression(data=frame, duration_col='t', event_col='e',
                   covariate_cols=['g'], model='weibull_aft')
    # one curve -> single-result blue
    fig1 = r.plot_km(show=False)
    l1 = [l for l in fig1.axes[0].get_lines()
          if not l.get_label().startswith('_')][0]
    assert l1.get_color() == PREDICTR_FIT_COLOR
    # >= 2 groups -> categorical palette, but with slots 2 and 3 swapped
    # for KM/NA (teal, purple, blue, ...), all solid
    fig3 = r.plot_km(by='g', ci=False, show=False)
    cols = [l.get_color() for l in fig3.axes[0].get_lines()
            if not l.get_label().startswith('_')]
    p = PREDICTR_PALETTE
    assert cols == [p[0], p[2], p[1]]


def test_matches_regression_engine():
    # same numbers whichever class computes them
    t = np.array(FAILURES + SUSPENSIONS, dtype=float)
    e = np.array([1] * len(FAILURES) + [0] * len(SUSPENSIONS))
    frame = pd.DataFrame({'t': t, 'e': e, 'x': np.arange(len(t), dtype=float)})
    r = Regression(data=frame, duration_col='t', event_col='e',
                   covariate_cols=['x'], model='weibull_aft')
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, show=False)
    np.testing.assert_allclose(r.kaplan_meier()['surv'].to_numpy(),
                               a.kaplan_meier()['surv'].to_numpy())
    np.testing.assert_allclose(r.nelson_aalen()['cumhaz'].to_numpy(),
                               a.nelson_aalen()['cumhaz'].to_numpy())
