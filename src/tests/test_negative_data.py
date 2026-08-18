"""
Tests for how negative df/ds values are handled per distribution.

Weibull/LogNormal/Exponential are all undefined for negative t (Weibull's
power terms, LogNormal's log, Exponential's rate model), so Analysis raises
for those. Normal is defined on the whole real line, so it's exempt.
PlotAll.compare() fits every distribution to the same data, so it falls
back to fitting Normal alone - rather than letting the other three raise -
whenever the data contains negative values, and says so in the subtitle.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from predictr import Analysis, PlotAll

NEGATIVE_DATA = [-2.5, 1.2, 3.4, 5.1, 2.2, 0.5, 4.4, -0.8, 3.9, 2.7]
POSITIVE_DATA = [54, 187, 216, 240, 244, 335, 361, 373, 375, 386]


@pytest.mark.parametrize('dist', ['weibull', 'lognormal', 'exponential'])
def test_analysis_rejects_negative_df_for_non_normal_dists(dist):
    with pytest.raises(ValueError, match='non-negative'):
        Analysis(df=NEGATIVE_DATA, dist=dist, show=False)


@pytest.mark.parametrize('dist', ['weibull', 'lognormal', 'exponential'])
def test_analysis_rejects_negative_ds_for_non_normal_dists(dist):
    with pytest.raises(ValueError, match='non-negative'):
        Analysis(df=POSITIVE_DATA, ds=[-1.0, 500.0], dist=dist, show=False)


def test_analysis_normal_accepts_negative_df():
    a = Analysis(df=NEGATIVE_DATA, dist='normal', show=False)
    a.mle()
    assert a.mu is not None and a.sigma is not None


def test_compare_falls_back_to_normal_only_on_negative_data():
    fig = PlotAll().compare(df=NEGATIVE_DATA, show=False, plot_pdf=False)
    assert len(fig.axes) == 1
    assert fig.axes[0].get_title().endswith('Normal')
    plt.close(fig)


def test_compare_subtitle_explains_normal_only_fallback():
    fig = PlotAll().compare(df=NEGATIVE_DATA, show=False, plot_pdf=False)
    subtitle_texts = [t.get_text() for t in fig.texts]
    assert any('Only Normal shown' in t for t in subtitle_texts)
    plt.close(fig)


def test_compare_still_fits_all_distributions_without_negative_data():
    fig = PlotAll().compare(df=POSITIVE_DATA, show=False, plot_pdf=False)
    assert len(fig.axes) == len(Analysis.SUPPORTED_DISTRIBUTIONS)
    plt.close(fig)
