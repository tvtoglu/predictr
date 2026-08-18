"""
Tests for the matplotlib figure lifecycle across predictr's plotting
methods. Every plotting method captures its own Figure (instead of relying
on bare pyplot module-level state), returns it to the caller, and closes it
exactly when show=True or save=True - but leaves it open when both are
False, since PlotAll.simple_weibull() depends on Analysis.plot() leaving a
show=False, save=False figure open so it can savefig() onto it afterwards.
Before this fix, plt.show() displayed every currently-registered open
figure (predictr never called plt.close()), so running several plots in one
script could resurrect old, already-dismissed figures.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from predictr import Analysis, PlotAll

FAILURES = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]


def _open_fignums():
    return set(plt.get_fignums())


@pytest.mark.parametrize('dist', ['weibull', 'normal', 'lognormal', 'exponential'])
def test_analysis_plot_closes_figure_on_show(dist):
    a = Analysis(df=FAILURES, dist=dist, show=True)
    a.mle()
    before = _open_fignums()
    fig = a.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before


@pytest.mark.parametrize('dist', ['weibull', 'normal', 'lognormal', 'exponential'])
def test_analysis_plot_leaves_figure_open_without_show_or_save(dist):
    a = Analysis(df=FAILURES, dist=dist, show=False)
    a.mle()
    before = _open_fignums()
    fig = a.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_analysis_plot_save_writes_file_and_closes(tmp_path):
    path = tmp_path / 'weibull.png'
    a = Analysis(df=FAILURES, show=False, save=True, path=str(path))
    a.mle()
    before = _open_fignums()
    fig = a.plot()
    assert path.exists()
    assert _open_fignums() == before


def test_analysis_plot_mrr_closes_figure_on_show():
    a = Analysis(df=FAILURES, bounds='bbb', show=True)
    a.mrr()
    before = _open_fignums()
    fig = a.plot_mrr()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before


def test_analysis_plot_mrr_leaves_figure_open_without_show_or_save():
    a = Analysis(df=FAILURES, bounds='bbb', show=False)
    a.mrr()
    before = _open_fignums()
    fig = a.plot_mrr()
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_plotall_mult_weibull_closes_on_show():
    a = Analysis(df=FAILURES, show=False)
    a.mle()
    before = _open_fignums()
    fig = PlotAll({'a': a}).mult_weibull(show=True)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before


def test_plotall_mult_weibull_leaves_figure_open_without_show():
    a = Analysis(df=FAILURES, show=False)
    a.mle()
    before = _open_fignums()
    fig = PlotAll({'a': a}).mult_weibull(show=False)
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_plotall_contour_plot_closes_on_show():
    a = Analysis(df=FAILURES, bounds='lrb', bounds_type='2s', show=False, cl=0.9)
    a.mle()
    before = _open_fignums()
    fig = PlotAll({'a': a}).contour_plot(show=True)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before


def test_plotall_contour_plot_leaves_figure_open_without_show():
    a = Analysis(df=FAILURES, bounds='lrb', bounds_type='2s', show=False, cl=0.9)
    a.mle()
    before = _open_fignums()
    fig = PlotAll({'a': a}).contour_plot(show=False)
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_plotall_contour_plot_show_weibull_returns_figure():
    a = Analysis(df=FAILURES, bounds='lrb', bounds_type='2s', show=False, cl=0.9)
    a.mle()
    fig = PlotAll({'a': a}).contour_plot(show_weibull=True)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plotall_weibull_pdf_closes_on_show():
    before = _open_fignums()
    fig = PlotAll().weibull_pdf(beta=[1.5], eta=[100.0], linestyle=['-'],
                                 x_bounds=[0, 200, 50], show=True)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before


def test_plotall_weibull_pdf_leaves_figure_open_without_show():
    before = _open_fignums()
    fig = PlotAll().weibull_pdf(beta=[1.5], eta=[100.0], linestyle=['-'],
                                 x_bounds=[0, 200, 50], show=False)
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_plotall_simple_weibull_returns_figure():
    before = _open_fignums()
    fig = PlotAll().simple_weibull(beta=1.5, eta=100.0)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert _open_fignums() == before | {fig.number}
    plt.close(fig)


def test_plotall_simple_weibull_save_writes_file(tmp_path):
    path = tmp_path / 'simple_weibull.png'
    fig = PlotAll().simple_weibull(beta=1.5, eta=100.0, save=True, path=str(path))
    assert path.exists()
    plt.close(fig)
