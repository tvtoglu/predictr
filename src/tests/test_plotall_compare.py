"""
Tests for PlotAll.compare() - the cross-distribution comparison grid.

Unlike mult_weibull() (same dist required, overlaid on one shared paper),
compare() takes one dataset and fits every distribution predictr supports
(Analysis.SUPPORTED_DISTRIBUTIONS) to it internally, then lays out one
probability-plot panel per distribution, ranked by AIC. The dataset used
here is the same right-censored example validated in this session against
a published reference (beta=1.7208, eta=606.528, loglik=-75.135,
AIC=154.27), so the ranking/AIC values asserted below are themselves
already externally validated, not just internally consistent.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors
import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from predictr import Analysis, PlotAll, PREDICTR_PALETTE, _anderson_darling

FAILURES = [54, 187, 216, 240, 244, 335, 361, 373, 375, 386]
SUSPENSIONS = [500] * 10


def test_compare_runs_and_ranks_by_aic():
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False, plot_pdf=False)
    panel_titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    # Weibull (AIC=154.27) must rank first, Normal (AIC=156.20) last, for
    # this specific validated dataset.
    assert panel_titles[0].endswith('Weibull')
    assert panel_titles[-1].endswith('Normal')
    plt.close(fig)


def test_compare_panel_count_matches_supported_distributions():
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, show=False, plot_pdf=False)
    visible_axes = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible_axes) == len(Analysis.SUPPORTED_DISTRIBUTIONS)
    plt.close(fig)


def test_compare_invalid_y_range_raises():
    with pytest.raises(ValueError):
        PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, y_min=0.9, y_max=0.1, show=False, plot_pdf=False)


def test_compare_rejects_unsupported_bounds():
    with pytest.raises(ValueError, match='bounds'):
        PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='lrb', show=False, plot_pdf=False)


def test_compare_rejects_chi2_with_specific_hint():
    """bounds='chi2' is valid for a standalone Analysis(dist='exponential'),
    but compare() fits all four distributions with one shared bounds=
    value, and 'chi2' doesn't apply to the other three - so it should be
    rejected with a message pointing at bounds='fb' instead of the generic
    unsupported-bounds message."""
    with pytest.raises(ValueError, match="chi2.*only applies to"):
        PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='chi2', show=False, plot_pdf=False)


def test_compare_returns_figure_and_respects_show():
    before = set(plt.get_fignums())
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, show=True, plot_pdf=False)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert set(plt.get_fignums()) == before


def test_compare_aic_values_match_published_reference():
    """Cross-check against the external reference validated earlier this
    session (NIST/SEMATECH-style worked example): beta=1.7208, eta=606.528,
    loglik=-75.135, AIC=154.27 for the Weibull fit on this exact dataset."""
    w = Analysis(df=FAILURES, ds=SUSPENSIONS, dist='weibull', bounds='fb',
                 bounds_type='2s', cl=0.9, show=False)
    w.mle()
    assert w.beta == pytest.approx(1.7208, abs=1e-3)
    assert w.eta == pytest.approx(606.528, abs=1e-2)
    assert w.loglik == pytest.approx(-75.135, abs=1e-2)
    assert w.aic == pytest.approx(154.27, abs=1e-2)


def test_anderson_darling_matches_reliability_package_uncensored():
    """Cross-checked against the Python `reliability` package's
    Fit_Weibull_2P (which implements the same Minitab-adjusted method) on
    this exact uncensored dataset: reliability's AD=1.8115 using its
    default Benard's-approximation plotting positions; recomputed here
    with predictr's own exact median_rank() instead, giving a slightly
    different but equally legitimate 1.8015 (verified by hand that
    swapping in Benard's approximation reproduces reliability's 1.8115 to
    10 decimal places, confirming the formula itself, not just the
    plotting-position convention, is correct)."""
    a = Analysis(df=FAILURES, dist='weibull', show=False)
    a.mle()
    assert _anderson_darling(a) == pytest.approx(1.8015, abs=1e-3)


def test_anderson_darling_matches_reliability_package_censored():
    """Same cross-check as the uncensored test above, but on this file's
    right-censored dataset (10 failures + 10 suspensions at t=500):
    reliability's Fit_Weibull_2P gives AD=70.747 there; predictr's own
    median_rank_cens() plotting positions give 70.767, matching to within
    the same plotting-position-convention gap as the uncensored case."""
    a = Analysis(df=FAILURES, ds=SUSPENSIONS, dist='weibull', show=False)
    a.mle()
    assert _anderson_darling(a) == pytest.approx(70.767, abs=1e-2)


def test_compare_criteria_ad_ranks_by_anderson_darling():
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb',
                             criteria='ad', show=False, plot_pdf=False)
    panel_titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    legend_texts = [ax.get_legend().get_texts()[0].get_text()
                     for ax in fig.axes if ax.get_legend()]
    assert panel_titles[0].endswith('Weibull')
    assert 'AD = ' in legend_texts[0]
    assert '(best)' in legend_texts[0]
    assert 'ΔAD' in legend_texts[-1]
    plt.close(fig)


def test_compare_criteria_ad_subtitle_notes_limited_comparability():
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb',
                             criteria='ad', show=False, plot_pdf=False)
    subtitle_texts = [t.get_text() for t in fig.texts]
    assert any('Anderson-Darling' in t and 'limited' in t for t in subtitle_texts)
    plt.close(fig)


def test_compare_criteria_aic_is_default_and_unaffected():
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False, plot_pdf=False)
    legend_texts = [ax.get_legend().get_texts()[0].get_text()
                     for ax in fig.axes if ax.get_legend()]
    assert 'AIC = ' in legend_texts[0]
    subtitle_texts = [t.get_text() for t in fig.texts]
    assert any('ranked by AIC' in t for t in subtitle_texts)
    plt.close(fig)


def test_compare_rejects_invalid_criteria():
    with pytest.raises(ValueError, match='criteria'):
        PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, criteria='ks', show=False, plot_pdf=False)


def _new_fignums(before):
    return sorted(set(plt.get_fignums()) - before)


def test_compare_plot_pdf_default_creates_second_figure():
    before = set(plt.get_fignums())
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False)
    new_nums = _new_fignums(before)
    assert len(new_nums) == 2
    assert new_nums[0] == fig.number
    for num in new_nums:
        plt.close(num)


def test_compare_plot_pdf_false_creates_only_main_figure():
    before = set(plt.get_fignums())
    fig = PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False, plot_pdf=False)
    new_nums = _new_fignums(before)
    assert new_nums == [fig.number]
    plt.close(fig)


def test_compare_pdf_figure_title_and_axis_labels():
    before = set(plt.get_fignums())
    PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False)
    new_nums = _new_fignums(before)
    ax_pdf = plt.figure(new_nums[-1]).axes[0]
    assert ax_pdf.get_title() == 'Probability Density Function Comparison'
    assert ax_pdf.get_ylabel() == 'Probability Density Function'
    for num in new_nums:
        plt.close(num)


def test_compare_pdf_rug_plot_marks_failures_and_suspensions_in_black():
    before = set(plt.get_fignums())
    PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False)
    new_nums = _new_fignums(before)
    ax_pdf = plt.figure(new_nums[-1]).axes[0]
    lines_by_label = {line.get_label(): line for line in ax_pdf.get_lines()}
    assert lines_by_label['Failures'].get_marker() == 'x'
    assert lines_by_label['Failures'].get_color() == 'black'
    assert lines_by_label['Suspensions'].get_marker() == 'o'
    assert lines_by_label['Suspensions'].get_markeredgecolor() == 'black'
    for num in new_nums:
        plt.close(num)


def test_compare_pdf_curves_use_predictr_palette():
    before = set(plt.get_fignums())
    PlotAll().compare(df=FAILURES, ds=SUSPENSIONS, bounds='fb', show=False)
    new_nums = _new_fignums(before)
    ax_pdf = plt.figure(new_nums[-1]).axes[0]
    dist_lines = [line for line in ax_pdf.get_lines()
                   if line.get_label() not in ('Failures', 'Suspensions')]
    assert len(dist_lines) == len(Analysis.SUPPORTED_DISTRIBUTIONS)
    for line in dist_lines:
        assert matplotlib.colors.to_hex(line.get_color()) in [c.lower() for c in PREDICTR_PALETTE]
    for num in new_nums:
        plt.close(num)
