#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tamertevetoglu
"""

from array import array
from logging import raiseExceptions
from math import floor, ceil
import copy
import colorsys
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
try:
    from matplotlib.cm import get_cmap
except ImportError:
    get_cmap = None
import pandas as pd
from scipy import optimize
from scipy.special import gamma, xlogy
from scipy.stats import norm, chi2, beta, linregress, trim_mean, expon
from scipy.stats.distributions import weibull_min
from scipy.spatial import ConvexHull


# predictr's categorical palette (up to 6 clearly distinct series) and the
# linestyles used to extend it beyond 6 - see _categorical_style() below.
PREDICTR_PALETTE = ['#008b8b', '#0077b2', '#810f7c', '#8b8b8b', '#fc4f30', '#e5c494']
PREDICTR_LINESTYLES = ['-', '--', ':', '-.']

# Single-result colour convention (a lone fit, not a set of series to
# tell apart): the MLE line and its markers on a probability plot. The
# confidence bounds and median-rank markers on a probability plot
# reuse this same colour (matching how PlotAll draws one dataset); a
# reference line stays 'grey'.
PREDICTR_FIT_COLOR = '#0077b2'

# Fixed candidate tick values for probability plots' y-axis. Which of these
# actually get drawn/labeled is controlled by Analysis(y_min=, y_max=) -
# only the candidates that fall within [y_min, y_max] are used, so the
# step spacing stays the same no matter how the visible range is narrowed.
PROBABILITY_PLOT_TICKS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01,
                                    0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3,
                                    0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
                                    0.99, 0.999])


def _categorical_style(n):
    """
    Returns n (color, linestyle) pairs for plotting n series/datasets.
    The first len(PREDICTR_PALETTE) series get one of predictr's palette
    colors each, all solid. Past that, colors repeat from the start of the
    palette but the linestyle advances to the next one in PREDICTR_LINESTYLES
    for every full pass through the palette - so up to
    len(PREDICTR_PALETTE) * len(PREDICTR_LINESTYLES) series stay
    distinguishable by color+shape before any (color, linestyle)
    combination repeats. A generated hue is never invented for a 7th+
    series, as that would silently blur into an existing one instead of
    reusing an already-distinct combination.
    """
    n_colors = len(PREDICTR_PALETTE)
    colors = [PREDICTR_PALETTE[i % n_colors] for i in range(n)]
    linestyles = [PREDICTR_LINESTYLES[(i // n_colors) % len(PREDICTR_LINESTYLES)] for i in range(n)]
    return colors, linestyles


# predictr's default matplotlib style, embedded here as a plain dict so
# 'predictr' works everywhere right after `pip install predictr`, with no
# extra file to ship/discover and no per-machine registration in
# matplotlib's stylelib directory required.
PREDICTR_STYLE = {
    'lines.linewidth': 4,
    'lines.solid_capstyle': 'butt',

    'legend.fancybox': True,
    'legend.frameon': True,
    'legend.facecolor': '#f7f7f7',
    'legend.edgecolor': '#bababa',
    'legend.framealpha': 0.9,
    'legend.borderpad': 0.5,
    'legend.title_fontsize': 10.0,

    'axes.prop_cycle': plt.cycler('color', PREDICTR_PALETTE),
    'axes.facecolor': '#fefefe',
    'axes.labelsize': 12,
    'axes.axisbelow': True,
    'axes.grid': True,
    'axes.edgecolor': '#c7c7c7',
    'axes.linewidth': 1.0,
    'axes.titlesize': 'x-large',
    'axes.titlepad': 15.0,
    'axes.labelpad': 8.0,

    'patch.edgecolor': '#f0f0f0',
    'patch.linewidth': 0.5,

    'svg.fonttype': 'path',

    'grid.linestyle': '-',
    'grid.linewidth': 1.0,
    'grid.color': '#c7c7c7',

    'xtick.major.size': 5,
    'xtick.major.width': 1.2,
    'xtick.major.pad': 6,
    'xtick.minor.size': 3,
    'xtick.minor.width': 0.8,
    'xtick.minor.visible': False,
    'xtick.color': '#333333',
    'xtick.labelsize': 10,

    'ytick.major.size': 5,
    'ytick.major.width': 1.2,
    'ytick.major.pad': 6,
    'ytick.minor.size': 0,
    'ytick.color': '#333333',
    'ytick.labelsize': 10,

    'font.size': 11.0,

    'figure.subplot.left': 0.08,
    'figure.subplot.right': 0.95,
    'figure.subplot.bottom': 0.07,
    'figure.facecolor': '#ffffff',
    'figure.dpi': 120,

    'savefig.dpi': 300,
}


# Every plotting method passes this to plt.figure(num=...)/plt.subplots(
# num=...) instead of leaving figure numbering to matplotlib's own default
# (lowest unused integer): since every plot closes its figure right after
# showing it (see _finish_plot() below), the default would hand out the
# same freed-up low number - usually 1 - to every subsequent plot in a
# script, so consecutive figures in one run would all be titled "Figure 1"
# in their window instead of counting up. This counter is process-wide
# (module-level) and monotonically increasing for the life of the process,
# so each new figure gets a number none of its predecessors ever had, even
# after they've been closed.
_figure_counter = itertools.count(1)


def _finish_plot(fig, show, save):
    """
    Shared show/close/return tail for every plotting method: shows fig if
    requested, then closes it whenever it was shown or saved (never
    unconditionally - PlotAll.simple_weibull() depends on Analysis.plot()
    leaving a show=False, save=False figure open so it can savefig() onto
    it afterwards) and returns it either way.

    The plt.pause(0.001) after close() matters on interactive GUI backends
    (e.g. MacOSX, Qt): plt.show() there only schedules a window to be
    drawn/raised and plt.close() only schedules it to be torn down: without
    pumping the event loop in between, a figure closed here can still be
    sitting in the backend's queue and visibly reappear - briefly, behind
    the next figure's window - once the *next* plot's plt.show() call
    finally pumps the event loop. It's a no-op on the non-interactive Agg
    backend used in tests.
    """
    if show:
        plt.show()
    if show or save:
        plt.close(fig)
        if show:
            plt.pause(0.001)
    return fig


def _km_na_life_table(t, d, cl):
    """Kaplan-Meier survival S(t) and Nelson-Aalen cumulative hazard H(t)
    for a right-censored sample - event/censor times t, event indicator d
    in {0, 1} - as one life table, one row per distinct observed time,
    with pointwise confidence limits at level cl.

    At each distinct time tau, with n at risk and k events there:
        S(tau) = prod (1 - k / n)        product-limit (Kaplan-Meier)
        H(tau) = sum  k / n              Nelson-Aalen
    Pointwise variances:
        Var(S) / S^2 = sum k / (n (n - k))          Greenwood
        Var(H)       = sum k / n^2
    The KM interval is formed on the scale ln(-ln S) (so it stays inside
    (0, 1)); the NA interval on ln H (so it stays positive). Returns a
    DataFrame: time, n_risk, n_event, n_censor, surv, surv_se, surv_lower,
    surv_upper, cumhaz, cumhaz_se, cumhaz_lower, cumhaz_upper.
    """
    t = np.asarray(t, dtype=float)
    d = np.asarray(d, dtype=float)
    order = np.argsort(t, kind='mergesort')
    t, d = t[order], d[order]
    n0 = len(t)
    z = norm.ppf(1.0 - (1.0 - cl) / 2.0)
    S, gw, H, vH = 1.0, 0.0, 0.0, 0.0
    rows = []
    for tau in np.unique(t):
        lo = int(np.searchsorted(t, tau, side='left'))
        hi = int(np.searchsorted(t, tau, side='right'))
        n_risk = n0 - lo
        blk = d[lo:hi]
        k = int(blk.sum())
        n_cens = len(blk) - k
        if k and n_risk:
            S *= 1.0 - k / n_risk
            if n_risk - k > 0:
                gw += k / (n_risk * (n_risk - k))
            H += k / n_risk
            vH += k / n_risk ** 2
        se_S = S * np.sqrt(gw)
        if 0.0 < S < 1.0 and gw > 0.0:
            se_th = np.sqrt(gw) / abs(np.log(S))
            S_lo = S ** np.exp(z * se_th)
            S_hi = S ** np.exp(-z * se_th)
        else:
            S_lo = S_hi = S
        se_H = np.sqrt(vH)
        if H > 0.0 and vH > 0.0:
            f = np.exp(z * se_H / H)
            H_lo, H_hi = H / f, H * f
        else:
            H_lo = H_hi = H
        rows.append((tau, n_risk, k, n_cens, S, se_S, S_lo, S_hi,
                     H, se_H, H_lo, H_hi))
    return pd.DataFrame(rows, columns=[
        'time', 'n_risk', 'n_event', 'n_censor',
        'surv', 'surv_se', 'surv_lower', 'surv_upper',
        'cumhaz', 'cumhaz_se', 'cumhaz_lower', 'cumhaz_upper'])


def _plot_km_na(obj, tables, kind, ci):
    """Step plot of the Kaplan-Meier survival (kind='km') or Nelson-Aalen
    cumulative hazard (kind='na'). `tables` is a list of (label, life
    table) pairs - one entry per group, label None meaning "no split".
    Uses obj's shared predictr plot attributes (fig_size, plot_style,
    fonts, legend, x_label, save/save_path). Returns the Figure; the
    caller passes it through the class's own show/save tail."""
    _apply_plot_style(obj.plot_style)
    fig = plt.figure(figsize=obj.fig_size, num=next(_figure_counter))
    ax = plt.gca()
    if kind == 'km':
        ycol, lo_c, hi_c = 'surv', 'surv_lower', 'surv_upper'
        y0, ylab, title = 1.0, 'survival  S(t)', 'Kaplan-Meier estimate'
    else:
        ycol, lo_c, hi_c = 'cumhaz', 'cumhaz_lower', 'cumhaz_upper'
        y0, ylab, title = 0.0, 'cumulative hazard  H(t)', \
            'Nelson-Aalen estimate'
    # One group: PREDICTR_FIT_COLOR (the single-result blue). Several
    # groups: the categorical palette (colour, and linestyle past the
    # 6th) as in PlotAll.mult_*() / Regression.plot_survival(), but with
    # palette slots 2 and 3 swapped for KM/NA - so a two-group split
    # reads teal vs. purple instead of teal vs. blue (the blue is too
    # close to PREDICTR_FIT_COLOR / the single-curve look).
    if len(tables) > 1:
        pal = list(PREDICTR_PALETTE)
        if len(pal) >= 3:
            pal[1], pal[2] = pal[2], pal[1]
        n = len(tables)
        colours = [pal[i % len(pal)] for i in range(n)]
        linestyles = [PREDICTR_LINESTYLES[(i // len(pal)) % len(PREDICTR_LINESTYLES)]
                      for i in range(n)]
    else:
        colours, linestyles = [PREDICTR_FIT_COLOR], ['-']
    for k, (label, tab) in enumerate(tables):
        c, ls = colours[k], linestyles[k]
        tt = np.concatenate(([0.0], tab['time'].to_numpy()))
        yy = np.concatenate(([y0], tab[ycol].to_numpy()))
        ax.step(tt, yy, where='post', color=c, linestyle=ls, linewidth=1.6,
                label=(label if label is not None else 'all units'))
        if ci:
            lo = np.concatenate(([y0], tab[lo_c].to_numpy()))
            hi = np.concatenate(([y0], tab[hi_c].to_numpy()))
            ax.fill_between(tt, lo, hi, step='post', color=c, alpha=0.15,
                            linewidth=0)
        m = tab['n_censor'].to_numpy() > 0
        if m.any():
            xt = tab['time'].to_numpy()[m]
            yt = np.interp(xt, tab['time'].to_numpy(), tab[ycol].to_numpy())
            ax.plot(xt, yt, '|', color=c, ms=7, mew=1.2, zorder=4)
    if kind == 'km':
        ax.set_ylim(0.0, 1.02)
    ax.set_xlabel(getattr(obj, 'x_label', None) or 'time',
                  fontsize=obj.xy_fontsize, color='black')
    ax.set_ylabel(ylab, fontsize=obj.xy_fontsize, color='black')
    ax.set_title(title, fontsize=obj.plot_title_fontsize, color='black')
    ax.tick_params(labelsize=obj.tick_fontsize)
    ax.grid(True, which='major')
    if getattr(obj, 'show_legend', True) and len(tables) > 1:
        ax.legend(fontsize=obj.legend_fontsize)
    plt.tight_layout()
    if obj.save:
        try:
            plt.savefig(obj.save_path)
        except Exception:
            raise ValueError('Path is faulty.')
    return fig


def _anderson_darling(obj):
    """
    Used by PlotAll.compare(criteria='ad'): an "adjusted" Anderson-Darling
    goodness-of-fit statistic for one fitted Analysis object - lower means
    the data sit closer to the fitted distribution.

    The classical Anderson-Darling statistic
        A^2 = -r - (1/r) * sum_i (2i-1) * [ln F(t_(i)) + ln(1 - F(t_(r+1-i)))]
    only works for a complete, uncensored sample of r observations: the
    weights (2i-1)/r are baked-in plotting positions for an *exact* rank i
    of r, and censoring breaks that assumption (some ranks are never
    observed as failures at all). This generalizes it by using whatever
    plotting positions (empirical CDF) actually apply - Kaplan-Meier, or
    here, predictr's own median_rank()/median_rank_cens() - instead of
    hard-coding (2i-1)/r, so the identical approach handles censored and
    uncensored data alike, with no distribution-specific or
    censoring-specific correction table needed.

    Derivation: A^2 is the definition
        A^2 = r * integral_0^1 [Fn(F^-1(u)) - u]^2 / [u * (1-u)] du
    (u = F(t), so Fn(F^-1(u)) is the empirical CDF re-expressed as a
    function of the fitted-CDF value instead of t). Fn is a step function
    that jumps at each failure's fitted-CDF value z_i = F(t_i); sorting
    those gives z_1 <= ... <= z_r, and between z_(i-1) and z_i (z_0 = 0),
    Fn is constant at q_(i-1) - the plotting position of the (i-1)-th
    failure (q_0 = 0, nothing observed yet). Expanding the integrand via
    partial fractions, (p-u)^2/[u(1-u)] = p^2/u + (1-p)^2/(1-u) - 1, which
    integrates in closed form to p^2*ln(u) - (1-p)^2*ln(1-u) - u. Summing
    that antiderivative's change across all r+1 intervals (the last one
    running up to a numerical stand-in for 1, since ln(1-u) diverges there)
    and multiplying by r gives A^2 directly - no numerical integration
    needed. xlogy(x, y) (= x*ln(y), defined as 0 when x == 0 even if y ==
    0) keeps the very first interval's p=0 term finite at u=0.

    Cross-checked numerically against an independent reference
    implementation of this same generalized-plotting-position approach, on
    both censored and uncensored data: matches to float precision once both
    use the same plotting-position convention (predictr's exact median rank
    vs. that implementation's default Benard's-approximation plotting
    positions give slightly different, but both legitimate, numbers -
    predictr uses its own convention here for consistency with what its
    probability plots already show).

    Like the AIC compare() otherwise ranks by, this is a *comparative*
    number only - no p-value is computed, since exact null distributions
    under estimated parameters and arbitrary censoring aren't tractable in
    general. Unlike AIC, this statistic is not reliably comparable *across
    different distributions* the way AIC's likelihood-based penalty is
    designed to be - PlotAll.compare() surfaces that caveat in its subtitle
    whenever criteria='ad'.
    """
    t = np.asarray(obj.df, dtype=float)
    if obj.dist == 'weibull':
        fitted_cdf = weibull_min.cdf(t, obj.beta, scale=obj.eta)
    elif obj.dist == 'exponential':
        fitted_cdf = expon.cdf(t, scale=obj.theta)
    elif obj.dist == 'normal':
        fitted_cdf = norm.cdf(t, obj.mu, obj.sigma)
    elif obj.dist == 'lognormal':
        fitted_cdf = norm.cdf(np.log(t), obj.mu, obj.sigma)

    plotting_positions = np.array(obj.median_rank() if obj.ds is None
                                    else obj.median_rank_cens())

    z = np.sort(np.clip(fitted_cdf, 1e-12, 1 - 1e-12))
    q = np.sort(plotting_positions)
    r = len(z)

    b_lo = np.insert(z, 0, 0.0)
    b_hi = np.append(z, 1 - 1e-12)
    p = np.insert(q, 0, 0.0)

    def antiderivative(u, p):
        return xlogy(p ** 2, u) - xlogy((1 - p) ** 2, 1 - u) - u

    return r * np.sum(antiderivative(b_hi, p) - antiderivative(b_lo, p))


def _apply_plot_style(style):
    """
    Applies a matplotlib plot style. 'predictr' is handled directly via
    PREDICTR_STYLE (no file/registration needed); anything else (a built-in
    style name, a custom style name registered in matplotlib's stylelib, or
    a path to a .mplstyle file) is passed through to plt.style.use().
    """
    if style == 'predictr':
        plt.rcParams.update(PREDICTR_STYLE)
    else:
        plt.style.use(style)


def _set_log_minor_ticks(ax):
    """
    Used by PlotAll.compare(), whose panels run smaller/lower-fontsize than
    a standalone plot(). Major ticks keep matplotlib's default log-axis
    "10^N" mathtext formatting (LogFormatterSciNotation - a coefficient of
    1 is dropped, so "1x10^2" reads as "10^2", but "2x10^2" keeps its "2").
    Minor ticks default to hidden (NullFormatter): LogFormatterSciNotation
    renders its "Nx10^M" minor-tick form visibly larger than a plain
    "10^N" major tick even at an identical fontsize (a matplotlib mathtext
    quirk, reproducible in a bare script with no predictr involved), which
    looks glaringly inconsistent once the panel is shrunk. But whenever the
    current view has at most one major (decade) tick in range - e.g. data
    clustered within one decade - there's no competing "10^N" label for the
    minor ones to look inconsistent next to, and hiding them there would
    leave the panel with no labeled tick at all if zoomed/panned
    interactively; minor labels are switched back on for that case only.
    """
    xlim = ax.get_xlim()
    n_major = floor(np.log10(xlim[1])) - ceil(np.log10(xlim[0])) + 1
    if n_major <= 1:
        ax.xaxis.set_minor_formatter(mpl.ticker.LogFormatterSciNotation())
    else:
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())


class Analysis:
    """
    Analysis provides parameter estimations, confidence bounds
    computations, bias corrections, and plotting of the data.
    """

    # Distributions Analysis can be configured for via dist=. Only 'weibull'
    # is implemented so far; this set is the single place future
    # distributions get registered as they're added.
    SUPPORTED_DISTRIBUTIONS = {'weibull', 'normal', 'lognormal', 'exponential'}

    @staticmethod
    def _get_cmap(name):
        """
        Version-safe colormap lookup. matplotlib >= 3.6 exposes
        mpl.colormaps as a subscriptable registry; older versions only
        have the (now-deprecated) get_cmap() function, where
        mpl.colormaps/plt.colormaps is a plain, non-subscriptable function.
        """
        try:
            return mpl.colormaps[name]
        except (TypeError, AttributeError):
            if get_cmap is None:
                raise
            return get_cmap(name)

    @staticmethod
    def _inv_mills_ratio(z):
        """
        Inverse Mills ratio h(z) = phi(z) / (1 - Phi(z)) of the standard
        normal distribution - the term that appears in the score and
        Hessian of the Normal/LogNormal log-likelihood for right-censored
        (suspended) observations, since d/dz[ln(1-Phi(z))] = -h(z) (see
        e.g. Lawless, "Statistical Models and Methods for Lifetime Data",
        2nd ed., ch. 3.2). Its derivative satisfies the identity
        h'(z) = h(z) * (h(z) - z), used by fisher_bounds()/mle() to build
        the observed Fisher information without a separate closed-form
        second derivative.

        Uses norm.sf(z) (scipy's own numerically stable tail probability,
        via erfc rather than a naive 1 - cdf(z) subtraction) instead of
        `1 - norm.cdf(z)` directly: the naive form catastrophically cancels
        to exactly 0.0 once z is only ~9 or so (norm.cdf(z) already rounds
        to 1.0 there), while norm.sf(z) stays accurate out to z ~ 37 before
        the true tail probability itself underflows float64's range.

        Beyond that point, norm.pdf(z) itself has *already* underflowed to
        0.0 (its exp(-z**2/2) term decays faster than norm.sf's erfc-based
        tail does), so the direct ratio collapses to 0/0 = nan right where
        h(z) is needed most - a suspension several tens of sigma out from
        the fitted mu/sigma is not a pathological input for MLE fitting on
        small, heavily censored samples (see e.g. mle()'s docstring for
        why this matters for the Normal/LogNormal Fisher bounds). Rather
        than let that silently poison the whole confidence-bound
        computation with a nan, values of z above the direct formula's
        breakdown point fall back to the classical asymptotic expansion of
        the Mills ratio, h(z) ~ z + 1/z (leading terms of Abramowitz &
        Stegun 26.2.13's continued-fraction expansion for the normal
        hazard function), which stays finite and within ~1e-4 relative
        error of the true value out there - letting the Fisher-information
        matrix stay finite (rather than nan) for exactly the inputs where
        the direct ratio would otherwise wipe it out.
        """
        z_in = np.asarray(z, dtype=float)
        z_flat = np.atleast_1d(z_in)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = norm.pdf(z_flat) / norm.sf(z_flat)
        unstable = ~np.isfinite(ratio) & (z_flat > 0)
        if np.any(unstable):
            z_u = z_flat[unstable]
            ratio[unstable] = z_u + 1.0 / z_u
        return ratio.reshape(z_in.shape) if z_in.ndim else ratio.item()

    @staticmethod
    def _point_estimate(val):
        """
        Returns the (beta, eta) point estimate of an Analysis object,
        preferring the bias-corrected values when a bias-correction method
        (bcm) was applied - mirrors the bcm dispatch used in fisher_bounds()
        and lrb(), where the bias-corrected point is also what the
        confidence region is centered on.
        """
        bcm = getattr(val, 'bcm')
        if bcm is None:
            return getattr(val, 'beta'), getattr(val, 'eta')
        elif bcm == 'c4':
            return getattr(val, 'beta_c4'), getattr(val, 'eta_c4')
        elif bcm == 'hrbu':
            return getattr(val, 'beta_hrbu'), getattr(val, 'eta_hrbu')
        elif bcm == 'np_bs':
            return getattr(val, 'beta_np_bs'), getattr(val, 'eta_np_bs')
        elif bcm == 'p_bs':
            return getattr(val, 'beta_p_bs'), getattr(val, 'eta_p_bs')
        else:
            raise ValueError(f'"{bcm}" is not a supported bias-correction method.')

    @staticmethod
    def _lrb_param_names(dist):
        """
        Which pair of attribute names holds the likelihood-ratio bounds'
        parameter-space coordinates for a given dist. Used by
        contour_plot() (in the PlotAll class below) to work with Weibull's
        (beta_lrb, eta_lrb) and Normal/LogNormal's (mu_lrb, sigma_lrb)
        through the same generic code, instead of hard-coding beta/eta.
        """
        if dist == 'weibull':
            return 'beta_lrb', 'eta_lrb'
        elif dist in ('normal', 'lognormal'):
            return 'mu_lrb', 'sigma_lrb'
        else:
            raise ValueError(f'contour_plot() has no likelihood-ratio '
                              f'parameterization for dist="{dist}" '
                              f'(Exponential has only one free parameter - '
                              f'a 2D contour doesn\'t apply - and bounds='
                              f'"lrb" isn\'t offered for it; see '
                              f'_exact_bounds_exponential()\'s docstring).')

    @classmethod
    def _lrb_arrays(cls, val):
        """
        Returns the (x, y) likelihood-ratio-bounds contour point arrays for
        an already-fitted Analysis object (val.lrb() must have been run),
        raising a clear error instead of contour_plot() crashing on a
        `None` if it wasn't.
        """
        x_name, y_name = cls._lrb_param_names(val.dist)
        x_lrb, y_lrb = getattr(val, x_name), getattr(val, y_name)
        if x_lrb is None or y_lrb is None:
            raise ValueError(f'This dist="{val.dist}" object has no '
                              f'likelihood-ratio bounds computed yet - run '
                              f'.mle() with bounds="lrb" (or call .lrb() '
                              f'directly) before contour_plot().')
        return x_lrb, y_lrb

    @classmethod
    def _lrb_point(cls, val):
        """
        Returns the (x, y) point estimate in the same parameter space as
        _lrb_arrays(), for the point-estimate marker: Weibull defers to
        _point_estimate() (bias-correction aware); Normal/LogNormal have
        no bias-correction methods (bcm is always None for them - see
        mle()'s config checks), so mu/sigma are used directly.
        """
        if val.dist == 'weibull':
            return cls._point_estimate(val)
        return val.mu, val.sigma

    @staticmethod
    def _lrb_axis_labels(dist):
        """
        Default (x, y) axis labels for contour_plot(), one pair per dist:
        Weibull's (beta, eta) are unitless-shape/real-time-scale; Normal's
        (mu, sigma) are both in the data's own units; LogNormal's (mu,
        sigma) are the Normal parameters of ln(t) - i.e. in log-space, not
        the data's units - so they're labeled distinctly from Normal's to
        avoid implying they're on the same numeric scale (see
        contour_plot()'s docstring on why dist families can't be mixed on
        one plot in the first place).

        LogNormal's subscript is written as \\ln(t) (what's being logged),
        not a bare \\ln, precisely so it can't be misread as "ln applied to
        mu" instead of its intended meaning "the mu that belongs to ln(t)"
        - mu itself is never logged, it's already estimated directly on
        ln(t) (see mle()'s LogNormal branch). This mirrors how the
        literature handles it: Meeker & Escobar's mu/sigma are always
        explicitly "of ln(T)" in text, and NIST's Engineering Statistics
        Handbook defines mu = ln(t50) outright - mu *is* a log-scale
        quantity by construction, not a log taken of something else - the
        same idea this subscript spells out explicitly instead of leaving
        implicit.
        """
        if dist == 'weibull':
            return r'$\widehat\beta$', r'$\widehat\eta$'
        elif dist == 'normal':
            return r'$\widehat\mu$', r'$\widehat\sigma$'
        elif dist == 'lognormal':
            return r'$\widehat\mu_{\ln(t)}$', r'$\widehat\sigma_{\ln(t)}$'
        raise ValueError(f'No axis labels defined for dist="{dist}".')

    @staticmethod
    def _vectorized_weibull_mle(samples):
        """
        Vectorized Weibull MLE for many independent uncensored samples at
        once. Mirrors Analysis.mle()'s beta_init() (analytic initial
        estimate) followed by a Newton-Raphson refinement using the same
        log-likelihood derivatives as ll_weib_beta_uncen / ll_weib_beta_uncen_2,
        but applied to a whole batch of (e.g. bootstrap) samples
        simultaneously instead of looping over them one at a time.

        Parameters
        ----------
        samples : ndarray, shape (B, n)
            Each row is one independent uncensored sample.

        Returns
        -------
        beta_vec, eta_vec : ndarray, shape (B,)
            MLE estimates of beta and eta for every row. Rows for which no
            valid estimate could be found are set to np.nan.
        """
        n = samples.shape[1]

        # Analytic initial estimate (row-wise equivalent of beta_init())
        t_N = samples.max(axis=1, keepdims=True)
        beta1 = (1.0 / n * np.sum(np.log(t_N / samples), axis=1)) ** -1

        x_k = (samples / t_N) ** beta1[:, None]
        ln_x_k = np.log(x_k)
        ln_x_k = np.where(ln_x_k == 0, 1e-8, ln_x_k)

        denom = np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1), axis=1)

        def sig(p):
            num = np.sum(x_k * ln_x_k ** (p - 1) * (ln_x_k ** 2 + p * ln_x_k + p), axis=1)
            return num / denom

        sig_0, sig_2, sig_3, sig_4, sig_5 = sig(0), sig(2), sig(3), sig(4), sig(5)

        beta_k = beta1 * (1 - sig_0 - sig_2 / 2 * sig_0 ** 2
                          + ((sig_3 / 6) - (sig_2 ** 2) / 2) * sig_0 ** 3
                          + ((5 * sig_3 * sig_2 / 12)
                             - (5 * sig_2 ** 3 / 8 - sig_4 / 24))
                          * sig_0 ** 4 + ((sig_5 / 120)
                                          + (7 * sig_2 ** 2 * sig_3 / 8)
                                          - (7 * sig_2 ** 4 / 8) - (sig_3 ** 2 / 12)
                                          - 3 * sig_2 * sig_4 / 24) * sig_0 ** 5)

        # Newton-Raphson refinement, vectorized across all rows, only
        # iterating on rows that have not converged yet.
        log_samples = np.log(samples)
        active = np.ones(len(beta_k), dtype=bool)
        for _ in range(100):
            if not np.any(active):
                break
            bk = beta_k[active, None]
            s = samples[active]
            ls = log_samples[active]
            pw = s ** bk
            a = ls.sum(axis=1)
            b = (pw * ls).sum(axis=1)
            c = pw.sum(axis=1)
            f = (1.0 / beta_k[active]) + (a / n) - (b / c)

            eta = (c / n) ** (1.0 / beta_k[active])
            ratio = s / eta[:, None]
            f2 = np.sum((-1.0 / bk ** 2) - ratio ** bk * (np.log(ratio) ** 2), axis=1)

            step = f / f2
            beta_k[active] = beta_k[active] - step

            idx = np.where(active)[0]
            active[idx[np.abs(step) < 1.0e-5]] = False

        beta_k = np.where(np.isfinite(beta_k) & (beta_k > 0), beta_k, np.nan)
        eta_k = (np.sum(samples ** beta_k[:, None], axis=1) / n) ** (1.0 / beta_k)
        return beta_k, eta_k

    @staticmethod
    def _vectorized_weibull_mrr(samples):
        """
        Vectorized Weibull median rank regression (MRR) for many
        independent uncensored samples at once, mirroring Analysis.mrr()
        for the uncensored case.

        The median ranks only depend on the sample size n (not on the data
        values), so they are computed once and reused as the fixed
        regressor for the OLS fit of every row.

        Parameters
        ----------
        samples : ndarray, shape (B, n)
            Each row is one independent uncensored sample, ascending sorted.

        Returns
        -------
        beta_vec, eta_vec : ndarray, shape (B,)
            MRR estimates of beta and eta for every row.
        """
        n = samples.shape[1]
        ranks = beta.ppf(0.5, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1)
        y = np.log(-np.log(1 - ranks))

        x = np.log(samples)
        x_mean = x.mean(axis=1, keepdims=True)
        y_mean = y.mean()
        y_centered = y - y_mean

        x_centered = x - x_mean
        cov = np.sum(x_centered * y_centered[None, :], axis=1)
        var_x = np.sum(x_centered ** 2, axis=1)

        beta_k = cov / var_x
        intercept = y_mean - beta_k * x_mean[:, 0]
        eta_k = np.exp(intercept / (-1 * beta_k))
        return beta_k, eta_k

    def __init__(self, df: list = None, ds: list = None, dist='weibull',
                 show: bool = False,
                 plot_style='predictr', bounds=None, bounds_type='2s',
                 cl=0.9, bcm=None, bs_size=5000, est_type='median',
                 unit='-', x_label = 'Time to Failure',
                 y_label = 'Unreliability', xy_fontsize=12, tick_fontsize=10,
                 plot_title_fontsize=14,
                 plot_title=None, plot_ranks=True, y_min=0.01, y_max=0.99,
                 fig_size=(6, 7), show_legend=True, legend_fontsize=9, save=False, **kwargs):
        """
        Parameters
        ----------
        df : list
            Contains failures. The default is None.
        ds : list
            Contains suspensions. The default is None.
        dist : string, optional
            Sets the distribution to fit: 'weibull', 'normal', 'lognormal',
            or 'exponential'. More distributions will be added over time.
            The default is 'weibull'. bcm=None is required for all but
            'weibull'. Bounds are limited too: 'normal'/'lognormal' only
            support bounds='fb' (Fisher-information-based). 'exponential'
            supports both bounds='fb' (the same asymptotic Fisher-
            information approach) and bounds='chi2' (an exact chi-square
            pivot - see _exact_bounds_exponential()'s docstring); the two
            give genuinely different bounds, not the same result under two
            names. bounds='chi2' is only accepted for dist='exponential',
            not the other three.
        show : bool, optional
            If True, plot will be shown. The default is False.
        plot_style : string, optional
            Sets. The default is 'predictr'.
        bounds : string, optional
            Sets the bounds method. The default is None.
        bounds_type : string, optional
            Sets the bounds type, e.g. two-sided, one-sided upper bound etc. The default is '2s'.
        cl : float, optional
            Stes confidence level between 0 and 1.0. The default is 0.9.
        bcm : string, optional
            Sets bias-correction method. The default is None.
        bs_size : int, optional
            Sets number of bootstrap samlpes. The default is 5000.
        est_type : float, optional
            Sets which statistic to compute from bootstrap samples. The default is 'median'.
        unit : string, optional
            Unit shown in the Weibull plot on the x-axis. The default is '-'.
        x_label : string, optional
            Label for the x-axis. The default is 'Time to Failure'.
        y_label : string, optional
            Label for the y-axis. The default is 'Unreliability'.
        xy_fontsize : float, optional
            Fontsize for the axes labels (xlabel/ylabel). The default is 12.
        tick_fontsize : float, optional
            Fontsize for the tick labels (the numbers on the axes). The default is 10.
        legend_fontsize : float, optional
            Fontsize for the legend. The default is 9.
        plot_title : string, optional
            Title for the plot. Defaults to "<Dist> Probability Plot" for
            the chosen dist (e.g. 'Weibull Probability Plot',
            'Normal Probability Plot') unless explicitly set.
        plot_title_fontsize : float, optional
            Fontsize of the plot title. The default is 14.
        y_min : float, optional
            Lower y-axis limit (unreliability, as a fraction) shown on the
            probability plot. Must satisfy 0 < y_min < y_max < 1. The
            default is 0.01. Only the fixed step values in
            PROBABILITY_PLOT_TICKS that fall within [y_min, y_max] are
            drawn, so narrowing the range doesn't change the tick spacing.
        y_max : float, optional
            Upper y-axis limit (unreliability, as a fraction) shown on the
            probability plot. Must satisfy 0 < y_min < y_max < 1. The
            default is 0.99.
        fig_size : tuple of floats, optional
            Sets width and height in inches: (width, height)
        save : boolean, optional
            If True, the plot is saved according to the path. The default is False.
        plot_style : string, optional
            Matplotlib plot style. The default is 'predictr'.
        plot_ranks : boolean, optional
            If True, median ranks will be plotted.
        show_legend : boolean, optional
            If True, the legend will be plotted. The default is True.
        **kwargs :
            path: string
                Path defines the directory and format of the figure E.g. r'var/user/.../test.pdf'

        """

        # Raise error if no data is given
        if (df is None and ds is None):
            raise ValueError('No data given. Please enter failures \
                             and/or suspensions.')

        # Raise error if show argument is not bool
        if not isinstance(show, bool):
            raise ValueError('Argument show must be of type bool.')

        # Raise error if dist is not supported
        if dist not in self.SUPPORTED_DISTRIBUTIONS:
            raise ValueError(f'"{dist}" is not a supported distribution. '
                              f'Supported distributions: '
                              f'{sorted(self.SUPPORTED_DISTRIBUTIONS)}')
        self.dist = dist

        # Raise error if the y-axis display range is not valid
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        self.y_min, self.y_max = y_min, y_max

        # Raise error if bounds type is not supported
        if bounds is None and (bounds_type != '2s'
                               or bounds_type != '1su'
                               or bounds_type != '1sl'):
            bounds_type = None

        if df is not None:
            self.df = sorted(df)
        else:
            self.df = df

        if ds is not None:
            self.ds = sorted(ds)
        else:
            self.ds = ds

        # Raise error if df/ds contain negative values for a distribution
        # that isn't defined there: Weibull's power terms, LogNormal's log,
        # and Exponential's rate model are all undefined for negative t.
        # Normal is the only one defined on the whole real line.
        if dist != 'normal':
            negative_data = [x for x in (self.df or []) + (self.ds or []) if x < 0]
            if negative_data:
                raise ValueError(
                    f'dist="{dist}" requires non-negative failure/suspension '
                    f'times, but negative values were found: {negative_data}. '
                    f'Only dist="normal" supports negative values.')

        # Plot related attributes
        self.show = show
        self.plot_style = plot_style
        if plot_title is None:
            _dist_titles = {'weibull': 'Weibull', 'normal': 'Normal',
                             'lognormal': 'Log-Normal', 'exponential': 'Exponential'}
            plot_title = f'{_dist_titles.get(dist, dist.capitalize())} Probability Plot'
        self.x_label, self.y_label, self.plot_title = x_label, y_label, plot_title
        self.xy_fontsize, self.plot_title_fontsize = xy_fontsize, plot_title_fontsize
        self.tick_fontsize = tick_fontsize
        self.show_legend, self.legend_fontsize = show_legend, legend_fontsize
        self.fig_size, self.plot_ranks, self.save = fig_size, plot_ranks, save
        if self.save:
            try:
                self.save_path = kwargs['path']
            except:
                raise ValueError('Path is not defined.')

        self.bounds_type = bounds_type
        self.bounds = bounds
        self.cl = cl
        self.bcm = bcm
        self.unit = unit
        self.bs_size = bs_size
        self.est_type = est_type
        self.beta_init = None
        self.censoring = None
        self.beta, self.eta = None, None
        self.mu, self.sigma = None, None
        self.theta = None
        self.loglik, self.aic = None, None
        self.beta_hrbu, self.eta_hrbu = None, None
        self.f, self.f_inv= None, None
        self.k_a_bound, self.se_beta,self.se_eta = None, None, None
        self.se_mu, self.se_sigma = None, None
        self.se_theta = None
        self.tmin_plot, self.tmax_plot, self.xplot = None, None, None
        self.unrel = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01,
                               0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4,
                               0.5, 0.6, 0.632, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        self.adj_ranks = None
        self.beta_c4, self.eta_c4 = None, None
        self.beta_p_bs, self.eta_p_bs = None, None
        self.beta_np_bs, self.eta_np_bs = None, None
        self.rvalue, self.pvalue = None, None
        self.bounds_upper, self.bounds_lower = None, None
        self.title = None
        self.mins, self.maxes, self.sol, self._z, self.z= None, None, None, None, None
        self.beta_lrb, self.eta_lrb, self.cl_lrb = None, None, None
        self.mu_lrb, self.sigma_lrb = None, None
        self.beta_pairs_lrb, self.eta_pairs_lrb = None, None
        self.sol_post, self.sol_b, self.sol_eta = None, None, None
        self.eta_range_init, self.beta_range_init = None, None
        self.beta_f_range, self.eta_f_range = None, None

    def mle(self):
        """
        mle() conducts a maximum likelihood estimation for the
        distribution set via dist= in the constructor. It handles
        uncensored and censored data.
        """
        # Check for configuration errors
        if (self.bounds is not None
            and self.bounds != 'fb'
            and self.bounds != 'lrb'
            and self.bounds != 'pbb'
            and self.bounds != 'npbb'
            and self.bounds != 'chi2'):
            raise ValueError(f'"{self.bounds}" is not supported by mle')
        if self.dist == 'exponential':
            if self.bounds not in (None, 'fb', 'chi2'):
                raise ValueError(f'"{self.bounds}" is not yet supported for '
                                  f'dist="exponential". Only bounds="fb" '
                                  f'(asymptotic Fisher-information bounds) '
                                  f'or bounds="chi2" (exact chi-square '
                                  f'bounds - see '
                                  f'_exact_bounds_exponential()\'s '
                                  f'docstring) are currently implemented.')
        elif self.dist in ('normal', 'lognormal') and self.bounds not in (None, 'fb', 'lrb'):
            raise ValueError(f'"{self.bounds}" is not yet supported for '
                              f'dist="{self.dist}". Only bounds="fb" '
                              f'(Fisher-information bounds) or bounds="lrb" '
                              f'(likelihood-ratio bounds - see lrb()\'s '
                              f'docstring) are currently implemented.')
        elif self.dist not in ('weibull', 'normal', 'lognormal') and self.bounds not in (None, 'fb'):
            raise ValueError(f'"{self.bounds}" is not yet supported for '
                              f'dist="{self.dist}". Only bounds="fb" is '
                              f'currently implemented.')
        elif self.dist == 'weibull' and self.bounds == 'chi2':
            raise ValueError('bounds="chi2" is only supported for '
                              'dist="exponential".')
        if self.dist != 'weibull' and self.bcm is not None:
            raise ValueError(f'Bias-correction methods (bcm) are not yet '
                              f'supported for dist="{self.dist}".')
        # Needed Log-Likelihood equations for uncensored data
        def ll_weib_beta_uncen(beta_, df):
            """
            First derivative of the Weibull log-likelihood function
            with respect to the shape parameter beta.
            """
            a = np.sum(np.fromiter((np.log(x) for x in df), float))
            b = np.sum(np.fromiter((x ** beta_ * np.log(x) for x in df), float))
            c = np.sum(np.fromiter((x ** beta_ for x in df), float))

            return (1 / beta_) + ((1 / len(df)) * a) - (b / c)

        def ll_weib_beta_uncen_2(beta_, df):
            """
            Second derivate of the Weibull log-likelihood function
            with respect to the shape parameter beta.
            """
            eta = (1 / len(df) * np.sum([x ** beta_ for x in df])) ** (1 / beta_)
            iterable = (((-1 / beta_ ** 2)
                         - (x / eta) ** beta_ * (np.log(x / eta) ** 2)) for x in df)
            return np.sum(np.fromiter((iterable), float))

        # Needed Log-Likelihood equations for right-censored data
        def ll_weib_beta_cen_r(beta_, df, ds):
            """
            First derivative of the Weibull log-likelihood function
            with respect to the shape parameter beta.
            """
            dat = df + ds
            a = np.sum(np.fromiter((np.log(x) for x in df), float))
            b = np.sum(np.fromiter((x ** beta_ * np.log(x) for x in dat), float))
            c = np.sum(np.fromiter((x ** beta_ for x in dat), float))

            return (1 / beta_) + ((1 / len(df)) * a) - (b / c)

        def ll_weib_beta_cen_2_r(beta_, df, ds):
            """
            Second derivate of the Weibull log-likelihood function
            with respect to the shape parameter beta.
            """
            dat = df + ds
            iter_eta = (x ** beta_ for x in dat)
            eta = ((1 / len(df)) * np.sum(np.fromiter(iter_eta, float))) ** (1 / beta_)
            iterable = ((-1 / beta_ ** 2)
                        - (x / eta) ** beta_ * (np.log(x / eta) ** 2) for x in df)
            iterable2 = (-(x / eta) ** beta_ * np.log(x / eta) ** 2 for x in ds)

            return np.sum(np.fromiter(iterable, float)) + np.sum(np.fromiter(iterable2, float))

        def np_bootstrap(sample, bs_size, est_type='median'):
            """
            Conducts a non-parametric bootstrap bias correction

            Parameters
            ----------
            sample : list
                Input sample from which bootstrap samples are drawn.
            bs_size : int
                Sets number of bootstraps samples to be drawn.
            est_type : str, optional
                Sets the statistic type:
                    1. median
                    2. mean
                    3. trimmed mean with alpha=0.1
                The default is 'median'.

            Returns
            -------
            beta_bs : float
                Bias-corrected Weibull shape parameter beta.
            eta_bs : float
                Bias-corrected Weibull shape parameter eta.
            """
            n = len(sample)
            sample_arr = np.asarray(sample, dtype=float)

            with np.errstate(divide='ignore', invalid='ignore'):
                # Draw all bootstrap resamples at once
                bs_samples = np.random.choice(sample_arr, size=(bs_size, n), replace=True)
                weib_beta, _ = self._vectorized_weibull_mle(bs_samples)

                # Safety net: a handful of resamples can fail to yield a
                # valid MLE estimate (e.g. degenerate resamples for small
                # sample sizes). Redraw and refit only those rows.
                invalid = ~np.isfinite(weib_beta)
                while np.any(invalid):
                    idx = np.where(invalid)[0]
                    redrawn = np.random.choice(sample_arr, size=(len(idx), n), replace=True)
                    bs_samples[idx] = redrawn
                    weib_beta[idx], _ = self._vectorized_weibull_mle(redrawn)
                    invalid = ~np.isfinite(weib_beta)

                # Calculate Bootstrap beta
                if est_type == 'median':
                    beta_bs = 2.0 * self.beta - np.median(weib_beta)
                elif est_type == 'mean':
                    beta_bs = 2.0 * self.beta - np.mean(weib_beta)
                elif est_type == 'trimmed_mean':
                    beta_bs = 2.0 * self.beta - trim_mean(weib_beta, 0.1)
                it = (i ** beta_bs for i in sample)
                eta_bs = (1 / len(sample) * np.sum(np.fromiter(it, float))) ** (1 / beta_bs)
                return beta_bs, eta_bs


        def p_bootstrap(sample, bs_size, est_type='median'):
            """
            Conducts a parametric bootstrap bias correction

            Parameters
            ----------
            sample : list
                Input sample from which bootstrap samples are drawn.
            bs_size : int
                Sets number of bootstraps samples to be drawn.
            est_type : str, optional
                Sets the statistic type:
                    1. median
                    2. mean
                    3. trimmed mean with alpha=0.1
                The default is 'median'.

            Returns
            -------
            beta_bs : float
                Bias-corrected Weibull shape parameter beta.
            eta_bs : float
                Bias-corrected Weibull shape parameter eta.

            """
            # Assign initial MLE of Weibull parameters
            beta_0 = self.beta
            eta_0 = self.eta
            n = len(sample)

            # Draw all bs_samples from initial estimation at once and
            # compute beta for every resample in one vectorized pass
            with np.errstate(divide='ignore', invalid='ignore'):
                bs_samples = weibull_min.rvs(beta_0, loc=0, scale=eta_0, size=(bs_size, n))
                weib_beta, _ = self._vectorized_weibull_mle(bs_samples)

            # Calculate Bootstrap beta
            if est_type == 'median':
                beta_bs = 2.0 * self.beta - np.median(weib_beta)
            elif est_type == 'mean':
                beta_bs = 2.0 * self.beta - np.mean(weib_beta)
            elif est_type == 'trimmed_mean':
                beta_bs = 2.0 * self.beta - trim_mean(weib_beta, 0.2)
            it = (i ** beta_bs for i in sample)
            eta_bs = (1 / len(sample) * np.sum(np.fromiter(it, float))) ** (1 / beta_bs)
            return beta_bs, eta_bs

        if self.dist == 'weibull':
            def beta_init(df, ds=None):
                """
                Analytic method for the initial estimation of beta.
                Best method so far, and may even replace Newton-Raphson Optimizer.
                """
                if ds is None:
                    ds = []
                data = df + ds
                t_N = max(data)
                n = len(df)

                term1 = np.array([np.log(t_N / t) for t in df])
                beta1 = ((1 / n) * np.sum(term1)) ** (-1)
                x_k = np.array([(t / t_N) ** beta1 for t in data])

                # Avoid dividing by zero
                ln_x_k = np.log(x_k)
                ln_x_k[ln_x_k == 0] = 1e-8

                # Calculate needed params
                sig_0 = np.sum(x_k * (ln_x_k ** (-1)) * (ln_x_k ** 2 + 0 * ln_x_k + 0)) \
                        / np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1))
                sig_2 = np.sum(x_k * ln_x_k ** (2 - 1) * (ln_x_k ** 2 + 2 * ln_x_k + 2)) \
                        / np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1))
                sig_3 = np.sum(x_k * ln_x_k ** (3 - 1) * (ln_x_k ** 2 + 3 * ln_x_k + 3)) \
                        / np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1))
                sig_4 = np.sum(x_k * ln_x_k ** (4 - 1) * (ln_x_k ** 2 + 4 * ln_x_k + 4)) \
                        / np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1))
                sig_5 = np.sum(x_k * ln_x_k ** (5 - 1) * (ln_x_k ** 2 + 5 * ln_x_k + 5)) \
                        / np.sum(x_k * (ln_x_k ** 2 + ln_x_k + 1))

                # Beta estimation
                beta_ = beta1 * (1 - sig_0 - sig_2 / 2 * sig_0 ** 2
                                 + ((sig_3 / 6) - (sig_2 ** 2) / 2) * sig_0 ** 3
                                 + ((5 * sig_3 * sig_2 / 12)
                                    - (5 * sig_2 ** 3 / 8 - sig_4 / 24))
                                 * sig_0 ** 4 + ((sig_5 / 120)
                                                 + (7 * sig_2 ** 2 * sig_3 / 8)
                                                 - (7 * sig_2 ** 4 / 8) - (sig_3 ** 2 / 12)
                                                 - 3 * sig_2 * sig_4 / 24) * sig_0 ** 5)
                return beta_

            # Check if input data is uncensored or censored
            if self.ds is None:
                # Compute beta_init
                self.beta_init = beta_init(self.df)

                # Define needed log-likelihood functions
                func = ll_weib_beta_uncen
                #fprime = ll_weib_beta_uncen_2

                # Conduct Secant method for finding the root
                # Add fprime = ll_weib_beta_uncen_2 for NR method
                self.beta = optimize.newton(func=func, x0=self.beta_init,
                                            args=(self.df,), tol=1.0e-5,
                                            maxiter=100)
                if self.beta <= 0:
                    print('b_init: {}'.format(self.beta_init))
                    raise ValueError('Beta estimation is false. \
                                 Check initial estimation beta_init')

                # Compute the Weibull scale parameter eta
                it = (x ** self.beta for x in self.df)
                self.eta = (1 / len(self.df)
                            * np.sum(np.fromiter(it, float))) ** (1 / self.beta)
            else:
                # Data is censored
                # Compute beta_init
                self.beta_init = beta_init(self.df, self.ds)

                # Determine the type of censoring
                self.censoring = 'right-censored'

                # Define needed log-likelihood functions
                if self.censoring == 'right-censored':
                    func = ll_weib_beta_cen_r
                    #fprime = ll_weib_beta_cen_2_r

                # Conduct Secant method for finding the root
                # Add fprime = ll_weib_beta_uncen_2 for NR method
                self.beta = optimize.newton(func=func, x0=self.beta_init,
                                            args=(self.df, self.ds,), tol=1.0e-5,
                                            maxiter=100)
                if self.beta <= 0:
                    print('b_init: {}'.format(self.beta_init))
                    raise ValueError('Beta estimation is false. \
                                 Check initial estimation beta_init')

                # Compute the Weibull scale parameter eta
                iter_eta = (x ** self.beta for x in self.df + self.ds)
                self.eta = ((1 / len(self.df))
                            * np.sum(np.fromiter(iter_eta, float))) ** (1 / self.beta)

            # Log-likelihood at the MLE and AIC (k=2 parameters: beta, eta).
            # weibull_min.logpdf/logsf evaluate the same density/survival
            # functions used to derive beta/eta above, just via scipy's
            # (well-tested, numerically robust) implementation rather than
            # a hand-written log/exp formula - this is evaluation of an
            # already-fitted density, not the estimation itself.
            self.loglik = float(np.sum(weibull_min.logpdf(self.df, self.beta, scale=self.eta)))
            if self.ds is not None:
                self.loglik += float(np.sum(weibull_min.logsf(self.ds, self.beta, scale=self.eta)))
            self.aic = 2 * 2 - 2 * self.loglik

        elif self.dist in ('normal', 'lognormal'):
            # LogNormal is fit as Normal on ln(data): if T is LogNormal(mu,
            # sigma) then ln(T) is Normal(mu, sigma) by definition, so mu/
            # sigma below are always the parameters of the underlying
            # Normal distribution of ln(T), same convention scipy.stats
            # uses. fisher_bounds() mirrors this (fit in log space, then
            # exponentiate the resulting time bounds back).
            if self.dist == 'lognormal':
                if any(x <= 0 for x in self.df + (self.ds or [])):
                    raise ValueError('LogNormal requires strictly positive '
                                     'data (df/ds).')
                fit_df = list(np.log(self.df))
                fit_ds = list(np.log(self.ds)) if self.ds is not None else None
            else:
                fit_df, fit_ds = self.df, self.ds

            # Score equations (d ln L / d mu, d ln L / d sigma) for the
            # Normal distribution with right-censored (suspended) data,
            # following Lawless ch. 3.2. For failures x_i and suspensions
            # s_j, with r = number of failures and z_j = (s_j - mu) / sigma:
            #   d lnL/d mu    = sum(x_i - mu)/sigma^2 + sum(h(z_j))/sigma
            #   d lnL/d sigma = -r/sigma + sum((x_i - mu)^2)/sigma^3
            #                   + sum(z_j * h(z_j))/sigma
            # where h is the inverse Mills ratio (_inv_mills_ratio). With
            # no suspensions both equations reduce to the closed-form
            # sample mean/variance MLE used directly below.
            def normal_score(params, df, ds):
                mu_, sigma_ = params
                resid = np.asarray(df, dtype=float) - mu_
                g_mu = np.sum(resid) / sigma_ ** 2
                g_sigma = -len(df) / sigma_ + np.sum(resid ** 2) / sigma_ ** 3
                if ds:
                    z = (np.asarray(ds, dtype=float) - mu_) / sigma_
                    h = self._inv_mills_ratio(z)
                    g_mu += np.sum(h) / sigma_
                    g_sigma += np.sum(z * h) / sigma_
                return [g_mu, g_sigma]

            def normal_score_jac(params, df, ds):
                """
                Jacobian of normal_score, i.e. the Hessian of ln L, built
                from h(z) via the identity h'(z) = h(z) * (h(z) - z) - see
                _inv_mills_ratio(). Reused by fisher_bounds() as the
                observed Fisher information at the MLE.
                """
                mu_, sigma_ = params
                resid = np.asarray(df, dtype=float) - mu_
                r = len(df)
                h_mumu = -r / sigma_ ** 2
                h_musigma = -2 * np.sum(resid) / sigma_ ** 3
                h_sigmasigma = r / sigma_ ** 2 - 3 * np.sum(resid ** 2) / sigma_ ** 4
                if ds:
                    z = (np.asarray(ds, dtype=float) - mu_) / sigma_
                    h = self._inv_mills_ratio(z)
                    h_prime = h * (h - z)
                    g = z * h
                    g_prime = h + z * h_prime
                    h_mumu += -np.sum(h * (h - z)) / sigma_ ** 2
                    h_musigma += (-np.sum(h) / sigma_ ** 2
                                  - np.sum(z * h_prime) / sigma_ ** 2)
                    h_sigmasigma += (-np.sum(g) / sigma_ ** 2
                                      - np.sum(z * g_prime) / sigma_ ** 2)
                return [[h_mumu, h_musigma], [h_musigma, h_sigmasigma]]

            if fit_ds is None:
                # Closed-form solution of the score equations above with
                # no suspension terms.
                data = np.asarray(fit_df, dtype=float)
                self.mu = data.mean()
                self.sigma = np.sqrt(np.mean((data - self.mu) ** 2))
            else:
                self.censoring = 'right-censored'
                mu_init = float(np.mean(fit_df))
                sigma_init = float(np.std(fit_df)) if len(fit_df) > 1 else 1.0

                # Maximize ln L directly (bounded L-BFGS-B, gradient =
                # normal_score) rather than root-find normal_score == 0
                # via fsolve: fsolve has no notion of "this is a
                # likelihood to climb," so on small, heavily censored
                # samples it can walk the unconstrained (mu, sigma) plane
                # straight through the true root and off to a spurious
                # stationary point of the *score* equations at
                # nonsensical values (e.g. sigma ~ 1e55) instead of
                # settling on the MLE - confirmed empirically on 5
                # failures + 3 suspensions, where fsolve's reported
                # "solution" varied wildly (including diverging) across
                # otherwise similar suspension times. Bounding sigma away
                # from 0 and maximizing the likelihood itself keeps the
                # search anchored to values that actually increase ln L,
                # so it can't wander into a region with a worse fit than
                # the starting guess.
                def neg_loglik(params):
                    mu_, sigma_ = params
                    ll = float(np.sum(norm.logpdf(fit_df, mu_, sigma_)))
                    if fit_ds:
                        ll += float(np.sum(norm.logsf(fit_ds, mu_, sigma_)))
                    return -ll

                def neg_loglik_jac(params):
                    g_mu, g_sigma = normal_score(params, fit_df, fit_ds)
                    return [-g_mu, -g_sigma]

                sigma_floor = max(sigma_init, 1.0) * 1e-6
                result = optimize.minimize(
                    neg_loglik, x0=[mu_init, sigma_init], jac=neg_loglik_jac,
                    method='L-BFGS-B', bounds=[(None, None), (sigma_floor, None)],
                    options={'ftol': 1e-15, 'gtol': 1e-12})
                if not result.success or result.x[1] <= 0:
                    raise ValueError('Sigma estimation is not valid. '
                                     'Check the input data.')

                # L-BFGS-B's own convergence tolerance leaves the score a
                # little short of machine-precision zero even once it's
                # unambiguously in the true root's basin (L-BFGS-B
                # approximates the Hessian; it doesn't use the exact one).
                # A couple of exact Newton corrections with the analytical
                # Hessian (normal_score_jac) - safe here precisely because
                # L-BFGS-B has already done the hard part of finding that
                # basin - polish the estimate the rest of the way, matching
                # the tight convergence fsolve gave on well-behaved inputs.
                params = np.array(result.x, dtype=float)
                for _ in range(4):
                    grad = np.asarray(normal_score(params, fit_df, fit_ds))
                    hess = np.asarray(normal_score_jac(params, fit_df, fit_ds))
                    try:
                        step = np.linalg.solve(hess, grad)
                    except np.linalg.LinAlgError:
                        break
                    candidate = params - step
                    if candidate[1] <= sigma_floor or not np.all(np.isfinite(candidate)):
                        break
                    params = candidate
                    if np.max(np.abs(grad)) < 1e-10:
                        break
                self.mu, self.sigma = params

            # Log-likelihood at the MLE and AIC (k=2: mu, sigma). For
            # LogNormal this is the log-likelihood of the *original* data T,
            # not of ln(T): since f_T(t) = f_lnT(ln t) / t (density of a
            # monotonic transform picks up a 1/|dt_lnT/dt| Jacobian term),
            # ln f_T(t_i) = norm.logpdf(ln t_i, mu, sigma) - ln(t_i) for
            # every failure. Without this term LogNormal's likelihood would
            # be evaluated in the wrong units and not be comparable via AIC
            # to Weibull/Normal/Exponential, which are evaluated in T
            # directly. The survival terms (suspensions) don't need the
            # Jacobian: P(T > s) = P(ln T > ln s) exactly, no density
            # involved.
            self.loglik = float(np.sum(norm.logpdf(fit_df, self.mu, self.sigma)))
            if self.dist == 'lognormal':
                self.loglik -= float(np.sum(np.log(self.df)))
            if fit_ds is not None:
                self.loglik += float(np.sum(norm.logsf(fit_ds, self.mu, self.sigma)))
            self.aic = 2 * 2 - 2 * self.loglik

        elif self.dist == 'exponential':
            # Closed-form MLE, valid for both complete and right-censored
            # (Type I/II) data (see e.g. Meeker & Escobar, "Statistical
            # Methods for Reliability Data", ch. 7.3): the log-likelihood
            #   lnL(theta) = -r*ln(theta) - (1/theta) * (sum(x_i) + sum(s_j))
            # (r = number of failures) has a single stationary point at
            #   theta_hat = total time on test / r,
            # i.e. total observed time divided by the number of failures -
            # no root-finding needed.
            r = len(self.df)
            total_time = sum(self.df) + (sum(self.ds) if self.ds is not None else 0)
            self.theta = total_time / r
            if self.ds is not None:
                self.censoring = 'right-censored'

            # Log-likelihood at the MLE and AIC (k=1: theta).
            self.loglik = float(np.sum(expon.logpdf(self.df, scale=self.theta)))
            if self.ds is not None:
                self.loglik += float(np.sum(expon.logsf(self.ds, scale=self.theta)))
            self.aic = 2 * 1 - 2 * self.loglik

        # Bias corrections
        if self.bcm is not None:
            if self.bcm == 'c4':
                if len(self.df) < 2:
                    raise ValueError('C4 bias correction requires two or more failures.')
                self.beta_c4 = self.beta * (np.sqrt(2 / (len(self.df) - 1)) *
                                            gamma(len(self.df) / 2)
                                            / gamma((len(self.df) - 1) / 2)) ** 3.52
                it = (x ** self.beta_c4 for x in self.df)
                self.eta_c4 = self.eta #(1 / len(self.df) * np.sum(np.fromiter(it, float))) ** (1 / self.beta_c4)
            elif self.bcm == 'hrbu':
                if self.ds is None:
                    self.beta_hrbu = (self.beta / (1.0115 + (1.278 / len(self.df))
                                                   + (2.001 / len(self.df) ** 2)
                                                   + (20.35 / len(self.df) ** 3)
                                                   - (46.98 / len(self.df) ** 4)))
                    it = (x ** self.beta_hrbu for x in self.df)
                    self.eta_hrbu = (1 / len(self.df)
                                     * np.sum(np.fromiter(it, float))) ** (1 / self.beta_hrbu)
                else:
                    self.beta_hrbu = (self.beta
                                      / (1 + (1.37 / ((len(self.df)) - 1.92)
                                              * np.sqrt((len(self.df) + len(self.ds))
                                                        / len(self.df)))))
                    iter_eta = (x ** self.beta_hrbu for x in self.df + self.ds)
                    self.eta_hrbu = ((1 / len(self.df))
                                     * np.sum(np.fromiter(iter_eta, float))) ** (1 / self.beta_hrbu)
            elif self.bcm == 'np_bs':
                if self.ds is None:
                    self.beta_np_bs, self.eta_np_bs = np_bootstrap(self.df,
                                                                   self.bs_size,
                                                                   self.est_type)
                else:
                    raise ValueError(f'"{self.bcm}" does not support suspensions yet')
            elif self.bcm == 'p_bs':
                if self.ds is None:
                    self.beta_p_bs, self.eta_p_bs = p_bootstrap(self.df,
                                                                self.bs_size,
                                                                self.est_type)
                else:
                    raise ValueError(f'"{self.bcm}" does not support suspensions yet')
            else:
                raise ValueError(f'"{self.bcm}" is not supported by mle')
        # Compute confidence bounds
        if self.bounds in ('fb', 'chi2'):
            self.fisher_bounds()
        elif self.bounds == 'lrb':
            self.lrb()
        elif self.bounds == 'npbb':
            self.npbb_bounds('mle')
        elif self.bounds == 'pbb':
            self.pb_bounds('mle')
        
        # Check save and show parameter
        if self.show or self.save:
            self.plot()

    # ------------------------------------------------------------------
    # non-parametric, distribution-free descriptions of the sample
    # ------------------------------------------------------------------
    def _km_na_input(self):
        """(times, event indicator) built from the failures (df, event)
        and right-censored suspensions (ds, censored). No DataFrame - the
        two lists the class is constructed with are enough."""
        fails = list(self.df) if self.df is not None else []
        susp = list(self.ds) if self.ds is not None else []
        if not fails:
            raise ValueError('kaplan_meier()/nelson_aalen() need at least '
                             'one failure in df.')
        t = np.asarray(fails + susp, dtype=float)
        d = np.array([1.0] * len(fails) + [0.0] * len(susp))
        return t, d

    def kaplan_meier(self, cl=None):
        """Non-parametric Kaplan-Meier survival estimate S(t) from the
        failures (df) and right-censored suspensions (ds) - no
        distribution assumed, no DataFrame needed.

        S(tau) = prod over event times <= tau of (1 - k / n), with k
        events out of n at risk; the pointwise interval uses the
        Greenwood standard error on the ln(-ln S) scale at level cl
        (default self.cl).

        Returns the life table as a DataFrame: time, n_risk, n_event,
        n_censor, surv, surv_se, surv_lower, surv_upper (a
        right-continuous step function - surv holds from `time` to the
        next row). With show=True (constructor) the step plot is drawn as
        well, in predictr style, exactly like mle()/mrr().
        """
        cl = float(self.cl if cl is None else cl)
        if not 0.0 < cl < 1.0:
            raise ValueError('cl must be between 0 and 1.')
        t, d = self._km_na_input()
        full = _km_na_life_table(t, d, cl)
        if self.show or self.save:
            _finish_plot(_plot_km_na(self, [(None, full)], 'km', True),
                         self.show, self.save)
        return full[['time', 'n_risk', 'n_event', 'n_censor',
                     'surv', 'surv_se', 'surv_lower', 'surv_upper']]

    def nelson_aalen(self, cl=None):
        """Non-parametric Nelson-Aalen cumulative-hazard estimate H(t)
        from df/ds - companion to kaplan_meier().

        H(tau) = sum over event times <= tau of k / n; the pointwise
        interval uses Var(H) = sum k / n^2 on the ln H scale at level cl.

        Returns a DataFrame: time, n_risk, n_event, n_censor, cumhaz,
        cumhaz_se, cumhaz_lower, cumhaz_upper. Draws the step plot too
        when show=True.
        """
        cl = float(self.cl if cl is None else cl)
        if not 0.0 < cl < 1.0:
            raise ValueError('cl must be between 0 and 1.')
        t, d = self._km_na_input()
        full = _km_na_life_table(t, d, cl)
        if self.show or self.save:
            _finish_plot(_plot_km_na(self, [(None, full)], 'na', True),
                         self.show, self.save)
        return full[['time', 'n_risk', 'n_event', 'n_censor',
                     'cumhaz', 'cumhaz_se', 'cumhaz_lower', 'cumhaz_upper']]

    def median_rank(self, cl=0.5):
        """
        Mediran ranks for uncensored data. Returns a list with
        median ranks.
        """
        ranks = []
        n = len(self.df)
        for i in range(1, n+1):
            ranks.append(beta.ppf(cl, i, n-i+1))
        return ranks

    def median_rank_cens(self, cl=0.5):
        """
        Returns adjusted median ranks as described in the
        Weibull Handbook. Returns a list with adjusted median ranks.
        """

        def bernard(adj_r, n, cl):
            """
            Returns Bernards Approximation for the adjusted ranks
            """
            #return (np.array(i) - 0.3) / (len(self.df+self.ds) + 0.4)
            return [beta.ppf(cl, i, n-i+1) for i in adj_r]

        n = len(self.df + self.ds)
        # Reverse ranks need to consider suspensions and their order
        all_ = self.df + self.ds
        rev_rank = []
        prev = 0
        for j in self.df:
            # Check if failure time is entered multiple times
            if self.df.count(j) > 1:
                # Ignore same elements after first time
                if prev == j:
                    pass
                else:
                    # Number of times element is in df
                    count_element = self.df.count(j)
                    # Loop through identical failure time
                    for i in range(count_element):
                        count = sum(map(lambda x : x < j, all_)) + i
                        rev_rank.append(len(all_) - count)
                prev = j
            else:
                count = sum(map(lambda x : x < j, all_))
                rev_rank.append(len(all_) - count)

        #Calculate adjusted rank
        adj_ranks = []
        prev_rank = 0
        for i in range(1, len(self.df)+1):
            adj_ranks.append((rev_rank[i-1] * prev_rank + n + 1) / (rev_rank[i-1] +1))
            prev_rank = adj_ranks[-1]
        self.adj_ranks = bernard(adj_ranks, n, cl)
        return self.adj_ranks

    def mrr(self):
        """
        mrr conducts the median rank regression

        Parameters
        ----------
        dist : str, optional
            Sets the distribution. The default is 'weibull'.

        Returns
        -------
        Median ranks and Binomial confidence bounds

        """
        if self.dist != 'weibull':
            raise ValueError(f'mrr() only supports dist="weibull" so far, '
                              f'not dist="{self.dist}". Use mle() instead.')

        # Check for configuration errors
        if (self.bounds is not None
            and self.bounds != 'bbb'
            and self.bounds != 'pbb'
            and self.bounds != 'npbb'
            and self.bounds != 'mcpb'):
            raise ValueError(f'"{self.bounds}" is not supported by mrr')

        # Compute the regression line using least squared method
        if self.ds is None:
            y_median_lnln = np.log(-np.log(1 - np.array(self.median_rank(0.5))))
            x_ = np.log(self.df)

            # Use linear regression
            ret = linregress(x_, y_median_lnln)
            self.beta = ret[0]
            intercept = ret[1]
            self.eta = np.exp(intercept / (-1 * self.beta))
            self.rvalue = ret[2] ** 2
            self.pvalue = ret[3]
        else:
            y_median_lnln = np.log(-np.log(1 - np.array(self.median_rank_cens(0.5))))
            x_ = np.log(self.df)

            # Use linear regression
            ret = linregress(x_, y_median_lnln)
            self.beta = ret[0]
            intercept= ret[1]
            self.eta = np.exp(intercept / (-1 * self.beta))
            self.rvalue = ret[2] ** 2
            self.pvalue = ret[3]

        # Check if bounds need to be computed
        # Compute confidence bounds
        if self.bounds == 'bbb':
            self.bb_bounds()
        elif self.bounds == 'pbb':
            self.pb_bounds('mrr')
        elif self.bounds == 'npbb':
            self.npbb_bounds('mrr')
        elif self.bounds == 'mcpb':
            self.mcp_bounds()

        # Show and/or save Weibull plot
        if self.show or self.save:
            self.plot_mrr()

    def bb_bounds(self):
        """
        Computes Beta-Binomial confidence bounds.
        """
        # Compute confidence bounds
        if self.ds is None:
            # Compute lower and upper bounds depending on bounds_type
            # and confidence level
            if self.bounds_type == '2s':
                self.bounds_lower = self.median_rank(0.5 - self.cl/2)
                self.bounds_upper = self.median_rank(0.5 + self.cl/2)
            elif self.bounds_type == '1su':
                self.bounds_upper = self.median_rank(self.cl)
            elif self.bounds_type == '1sl':
                self.bounds_lower= self.median_rank(1 - self.cl)
        else:
            # Compute lower and upper bounds depending on bounds_type
            # and confidence level
            if self.bounds_type == '2s':
                self.bounds_lower = self.median_rank_cens(0.5 - self.cl/2)
                self.bounds_upper = self.median_rank_cens(0.5 + self.cl/2)
            elif self.bounds_type == '1su':
                self.bounds_upper = self.median_rank_cens(self.cl)
            elif self.bounds_type == '1sl':
                self.bounds_lower= self.median_rank_cens(1 - self.cl)

    def pb_bounds(self, method_call):
        """
        Computes parametric bootstrap confidence bounds.
        """
        # Use the initial estimation of Weibull parameters
        beta_0 = self.beta
        eta_0 = self.eta
        n = len(self.df)

        # Draw all bootstrap resamples at once and fit every one of them
        # in a single vectorized pass instead of looping bs_size times.
        bs_samples = weibull_min.rvs(beta_0, loc=0, scale=eta_0, size=(self.bs_size, n))

        if method_call == 'mrr':
            bs_samples.sort(axis=1)
            beta_bs, eta_bs = self._vectorized_weibull_mrr(bs_samples)
        elif method_call == 'mle':
            beta_bs, eta_bs = self._vectorized_weibull_mle(bs_samples)
        else:
            raise ValueError(f'pb_bounds() does not support {method_call}')

        # Percentile life for every bootstrap fit and every unreliability
        # level in self.unrel, then sort each unreliability column
        neg_log = -np.log(1 - self.unrel)
        life_matrix = eta_bs[:, None] * (neg_log[None, :] ** (1.0 / beta_bs[:, None]))
        life_matrix.sort(axis=0)

        # Compute iloc position of lower and upper coonfidence bounds
        # -1 necessary, since python indexing starts at 0
        if self.bounds_type == '2s':
            lower_perc_position = ceil(self.bs_size * ((1 - self.cl) / 2)) - 1
            upper_perc_position = floor(self.bs_size * (1 - ((1 - self.cl) / 2))) - 1
            self.bounds_lower = life_matrix[lower_perc_position].tolist()
            self.bounds_upper = life_matrix[upper_perc_position].tolist()
        elif self.bounds_type == '1su':
            upper_perc_position = floor(self.bs_size * self.cl) - 1
            self.bounds_upper = life_matrix[upper_perc_position].tolist()
        elif self.bounds_type == '1sl':
            lower_perc_position = ceil(self.bs_size * (1 - self.cl)) - 1
            self.bounds_lower = life_matrix[lower_perc_position].tolist()

    def npbb_bounds(self, method_call):
        """
        Computes non-parametric bootstrap confidence bounds.
        """
        def cen_index(df, ds):
            'Solely censored samples need this index for resampling'
            # Create dat with length: df+ds and input tuples
            # check if only censored data is available
            if df is not None:
                # index 1 -> failure
                df_w_index = [(i, 1) for i in df]
                # index 0 -> suspension
                ds_w_index = [(i, 0) for i in ds]

            # create new list of tuples with all information
            dat = df_w_index + ds_w_index
            return dat

        if method_call not in ('mrr', 'mle'):
            raise ValueError(f'npbb_bounds() does not support {method_call}')

        neg_log = -np.log(1 - self.unrel)

        if self.ds is None:
            # Uncensored data: every bootstrap resample has the same
            # size, so all resamples can be drawn and fitted in a single
            # vectorized pass instead of looping bs_size times.
            n = len(self.df)
            sample_arr = np.asarray(self.df, dtype=float)

            def fit(samples):
                if method_call == 'mrr':
                    return self._vectorized_weibull_mrr(np.sort(samples, axis=1))
                return self._vectorized_weibull_mle(samples)

            with np.errstate(divide='ignore', invalid='ignore'):
                bs_samples = np.random.choice(sample_arr, size=(self.bs_size, n), replace=True)
                beta_bs, eta_bs = fit(bs_samples)

                # Safety net: redraw and refit any resample that failed
                # to produce a valid estimate (mirrors the try/except
                # retry of the original loop-based implementation).
                invalid = ~(np.isfinite(beta_bs) & np.isfinite(eta_bs))
                while np.any(invalid):
                    idx = np.where(invalid)[0]
                    redrawn = np.random.choice(sample_arr, size=(len(idx), n), replace=True)
                    bs_samples[idx] = redrawn
                    beta_bs[idx], eta_bs[idx] = fit(redrawn)
                    invalid = ~(np.isfinite(beta_bs) & np.isfinite(eta_bs))

            life_matrix = eta_bs[:, None] * (neg_log[None, :] ** (1.0 / beta_bs[:, None]))
            life_matrix.sort(axis=0)
        else:
            # Censored data: each bootstrap resample can contain a
            # different number of failures/suspensions, which does not
            # vectorize as cleanly, so this path keeps the original
            # per-resample loop.
            dat = cen_index(self.df, self.ds)
            rows = []
            j = 0
            with np.errstate(divide='ignore', invalid='ignore'):
                while j < self.bs_size:
                    try:
                        # np.random.choice requires a 1darray,
                        # which dat is not after transforming it to an array
                        # Randomly draw indices instead
                        bs_samples_idx = np.random.choice(len(dat),
                                                            size=len(dat),
                                                            replace=True)
                        # Use randomly drawn indices to generate random bootstrap sample
                        bs_samples = np.array(dat)[bs_samples_idx]

                        # Filter failures and suspenions by ID (0 or 1)
                        df_temp = [i for i, k in bs_samples if k == 1]
                        ds_temp = [i for i, k in bs_samples if k == 0]

                        # Conduct MRR/MLE to compute Weibull parameters
                        y = Analysis(df=df_temp, ds=ds_temp)
                        if method_call == 'mrr':
                            y.mrr()
                        else:
                            y.mle()
                        rows.append(np.array(y.eta) *
                                    ((-np.log(1 - self.unrel)) ** (1 / np.array(y.beta))))
                        j += 1
                    except Exception:
                        pass
            life_matrix = np.sort(np.array(rows), axis=0)

        # Compute iloc position of lower and upper coonfidence bounds
        # -1 necessary, since python indexing starts at 0
        if self.bounds_type == '2s':
            lower_perc_position = ceil(self.bs_size * ((1 - self.cl) / 2)) - 1
            upper_perc_position = floor(self.bs_size * (1 - ((1 - self.cl) / 2))) - 1
            self.bounds_lower = life_matrix[lower_perc_position].tolist()
            self.bounds_upper = life_matrix[upper_perc_position].tolist()
        elif self.bounds_type == '1su':
            upper_perc_position = floor(self.bs_size * self.cl) - 1
            self.bounds_upper = life_matrix[upper_perc_position].tolist()
        elif self.bounds_type == '1sl':
            lower_perc_position = ceil(self.bs_size * (1 - self.cl)) - 1
            self.bounds_lower = life_matrix[lower_perc_position].tolist()

    def mcp_bounds(self):
        """
        Computes confidence bounds for mrr() using Monte Carlo pivotals.
        Not to be used with mle()!

        """
        n = len(self.df)

        # Fixed params are needed in this method: beta=eta=1.0. Draw all
        # bs_size samples at once and fit every one of them in a single
        # vectorized pass instead of looping bs_size times.
        samples = weibull_min.rvs(1.0, loc=0, scale=1.0, size=(self.bs_size, n))
        samples.sort(axis=1)
        beta_k, eta_k = self._vectorized_weibull_mrr(samples)

        z_p = ((np.log(eta_k)[:, None] - np.log(np.log(1 / (1 - self.unrel)))[None, :])
               * beta_k[:, None])
        z_p.sort(axis=0)

        # Compute iloc position of lower and upper coonfidence bounds
        # -1 necessary, since python indexing starts at 0
        # Compute time intervals from z_p
        if self.bounds_type == '2s':
            lower_perc_position = ceil(self.bs_size * ((1 - self.cl) / 2)) - 1
            upper_perc_position = floor(self.bs_size * (1 - ((1 - self.cl) / 2))) - 1

            # Get lower and upper z_p per percentile and compute time intervals
            bounds_lower_z_p = z_p[lower_perc_position]
            bounds_upper_z_p = z_p[upper_perc_position]

            # Actual bounds as timestamps
            # ATTENTION: upper_z_p results in lower_t_p bounds
            self.bounds_upper = list(np.exp(np.log(self.eta) - bounds_lower_z_p / self.beta))
            self.bounds_lower = list(np.exp(np.log(self.eta) - bounds_upper_z_p / self.beta))
        elif self.bounds_type == '1su':
            lower_perc_position = ceil(self.bs_size * ((1 - self.cl) / 2)) - 1

            # Get lower and upper z_p per percentile and compute time intervals
            bounds_lower_z_p = z_p[lower_perc_position]

            # Actual bounds as timestamps
            # ATTENTION: upper_z_p results in lower_t_p bounds
            self.bounds_upper = list(np.exp(np.log(self.eta) - bounds_lower_z_p / self.beta))
        elif self.bounds_type == '1sl':
            upper_perc_position = floor(self.bs_size * (1 - ((1 - self.cl) / 2))) - 1

            # Get lower and upper z_p per percentile and compute time intervals
            bounds_upper_z_p = z_p[upper_perc_position]

            # Actual bounds as timestamps
            # ATTENTION: upper_z_p results in lower_t_p bounds
            self.bounds_lower = list(np.exp(np.log(self.eta) - bounds_upper_z_p / self.beta))

    def plot_mrr(self):
        """
        Plots Weibull Probability Plot for Median Rank Regression
        """
        # Some needed functions
        def weibull_prob_paper(x):
            """
            Needed to adjust figure to the Weibull probability plot.
            """
            x = np.asarray(x)

            # Prevent np.log(0) error raise
            x[x > .9999] = np.nan
            return np.log(-np.log(1 - x))

        # Just for y_tickslabel on the y-axis
        def weibull_ticks(y_i, _):
            """
            Adjusts the y-axis tick labels
            """
            return '{:.1f}'.format((100 * (1 - np.exp(-np.exp(y_i)))))

        def unrel_func(x_est, beta_, eta):
            if type(x_est) == list:
                x_est = np.asarray(x_est)
            y_est = (1 - np.exp(-(x_est / eta) ** beta_))
            y_est_lnln = weibull_prob_paper(y_est)
            return y_est_lnln


        def inverse_weibull(perc, beta, eta):
            """
            Computes time to failure data points.
            This function is being used to plot Weibull lines.

            Parameters
            ----------
            perc : float
                Percentage points fo which time to failure data should be computed.
            beta : float
                Weibull shape parameter.
            eta : float
                Weibull scale parameter.

            Returns
            -------
            float
                Time to failure data points.

            """
            return ((-np.log(1 -perc)) ** (1 / beta)) * eta

        # Generate Weibull Plot Figure
        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= self.y_min)
                                          & (PROBABILITY_PLOT_TICKS <= self.y_max)]
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.yticks(lny_ticks, color='black')
        ax.set_yticks([weibull_prob_paper(0.632)], minor=True)

        # Plots the horizontal dashed line for 63.2%
        plt.grid(True, which='minor', axis='y', linestyle='--')

        # X-Axis scaling
        if self.bounds is None:
            self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
            self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta
        elif self.bounds == 'bbb':
            self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
            self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta
        elif self.bounds == 'pbb' :
            if self.bounds_type == '2s':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1su':
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1sl':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta
        elif self.bounds == 'npbb' :
            if self.bounds_type == '2s':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1su':
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1sl':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta
        elif self.bounds == 'mcpb' :
            if self.bounds_type == '2s':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1su':
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
                self.tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1sl':
                self.tmin_plot = min(self.bounds_lower)
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta

        self.xplot = np.linspace(self.tmin_plot, self.tmax_plot, 100)
        left = (10 ** (np.ceil(np.log10(self.tmin_plot)) - 1))
        right = (10 ** (np.ceil(np.log10(self.tmax_plot))))
        plt.xlim(left, right)
        plt.tick_params(axis='x', colors='black')

        # Set labels and legends
        plt.title(self.plot_title, color='black', fontsize=self.plot_title_fontsize)
        plt.xlabel(f'{self.x_label}{" in "+self.unit if self.unit!="-" else ""}', color='black', fontsize=self.xy_fontsize)
        plt.xticks(fontsize=self.tick_fontsize)
        plt.ylabel(self.y_label + ' in %', color='black', fontsize=self.xy_fontsize)
        plt.yticks(fontsize=self.tick_fontsize)

        # Plot legend
        if self.ds is None:
            susp_num = 0
        else:
            susp_num = len(self.ds)

        # Plot median MRR line
        xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
        plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                     color=PREDICTR_FIT_COLOR, linestyle='-',
                     linewidth=1.5, zorder=2)

        # Adapt legend strings
        if self.bounds == 'bbb':
            bounds_legend = 'Beta-Binomial bounds'
        elif self.bounds == 'pbb':
            bounds_legend = 'Par. Bootstrap bounds'
        elif self.bounds == 'npbb':
            bounds_legend = 'Non-Par. Bootstrap bounds'
        elif self.bounds == 'mcpb':
            bounds_legend = 'MC pivotal bounds'
        if self.ds is not None:
            self.title = 'Adj. MRR'
        else:
            self.title = 'MRR'


        # Check for bounds and plot
        if self.bounds is not None:
            if self.bounds == 'bbb':
                if self.ds is not None:
                    med_ranks = self.median_rank_cens()
                else:
                    med_ranks = self.median_rank()
            # Plot legend
            if self.bounds_type == '2s':
                if self.bounds == 'bbb':
                    yerr_lower = (weibull_prob_paper(med_ranks)
                                  - weibull_prob_paper(self.bounds_lower))
                    yerr_upper = (weibull_prob_paper(self.bounds_upper)
                                  -  weibull_prob_paper(med_ranks))
                    plt.errorbar(x = self.df, y=weibull_prob_paper(med_ranks),
                                 yerr=[yerr_lower, yerr_upper], fmt='none', ecolor= PREDICTR_FIT_COLOR,
                                 capsize=5, markeredgewidth=2, alpha=0.5)
                    plt.fill_between(x = self.df,
                                     y1=weibull_prob_paper(self.bounds_lower),
                                     y2=weibull_prob_paper(self.bounds_upper),
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n2s @{}% (pctl)'.format((bounds_legend), self.cl * 100)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'pbb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n2s @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'npbb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n2s @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'mcpb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
            elif self.bounds_type == '1su':
                if self.bounds == 'bbb':
                    yerr_upper = (weibull_prob_paper(self.bounds_upper)
                                  -  weibull_prob_paper(med_ranks))
                    plt.errorbar(x = self.df, y=weibull_prob_paper(med_ranks),
                                 yerr=[len(self.df) * [0],yerr_upper], fmt='none',
                                 ecolor= PREDICTR_FIT_COLOR, capsize=5, markeredgewidth=2, alpha=0.5)
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1su @{}% (pctl)'.format((bounds_legend), self.cl * 100)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'pbb':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1su @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'npbb':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1su @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'mcpb':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1su @{}%'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
            elif self.bounds_type == '1sl':
                if self.bounds == 'bbb':
                    yerr_upper = (weibull_prob_paper(self.bounds_lower)
                                  -  weibull_prob_paper(med_ranks))
                    plt.errorbar(x = self.df, y=weibull_prob_paper(med_ranks),
                                 yerr=[len(self.df) * [0],yerr_upper], fmt='none',
                                 ecolor= PREDICTR_FIT_COLOR, capsize=5, markeredgewidth=2, alpha=0.5)
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                      len(self.df), susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                    + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                    + 'pval={:.2e}'.format(self.pvalue),
                                    '\n{}:\n1sl @{}% (pctl)'.format((bounds_legend),
                                                                    self.cl * 100)],
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'pbb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1sl @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'npbb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1sl @{}% (pctl)'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
                elif self.bounds == 'mcpb':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)

                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                                  len(self.df), susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                                + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                                + 'pval={:.2e}'.format(self.pvalue),
                                '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)
                                +'\nBS samples: {}'.format(self.bs_size)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=self.title)
        else:
            plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                          len(self.df), susp_num)
                        + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                        + r'$\widehat\eta={:.3f}$ '.format(self.eta)
                        + '\n$r^2={:.3f}$ | '.format(self.rvalue)
                        + 'pval={:.2e}'.format(self.pvalue)],
                        loc='lower left', bbox_to_anchor=(0.65, 0.0),
                        fontsize=self.legend_fontsize, title=self.title)

        # Plot discrete median ranks
        if self.plot_ranks:
            if self.ds is None:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank()), marker='o',
                             markerfacecolor=PREDICTR_FIT_COLOR, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)
            else:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank_cens()), marker='o',
                             markerfacecolor=PREDICTR_FIT_COLOR, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder= 3)

        # Pin the y-axis to [y_min, y_max] as the very last step: fill_between/
        # fill_betweenx above auto-expand the view when they're added
        # (matplotlib's add_collection autoscales regardless of an earlier
        # plt.ylim() call), so setting this any earlier would get silently
        # overridden once the bounds shading is drawn.
        ax.set_ylim(bottom=weibull_prob_paper(self.y_min),
                    top=weibull_prob_paper(self.y_max))
        plt.tight_layout()

        # bbox_to_anchor=(0.65, 0.0) above can push the legend past the
        # right edge of the figure depending on bounds/bounds_type; widen
        # the figure so the legend stays fully visible instead of being
        # clipped. Measured via a throwaway Agg renderer so this works
        # reliably regardless of the active display backend (interactive
        # backends like macosx don't reliably populate a renderer via a
        # plain canvas.draw()). Must run AFTER tight_layout(): tight_layout
        # tightens the axes margins, shifting the (axes-relative) legend
        # anchor further right - measuring before it would just get
        # squeezed away again. Deliberately not calling tight_layout()
        # again afterwards, so the added width stays as real margin instead
        # of being re-absorbed.
        legend = plt.gca().get_legend()
        if legend is not None:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            original_canvas = fig.canvas
            FigureCanvasAgg(fig)
            fig.canvas.draw()
            legend_bbox = legend.get_window_extent(fig.canvas.get_renderer())
            fig.canvas = original_canvas
            fig_px_width = fig.get_size_inches()[0] * fig.dpi
            if legend_bbox.x1 > fig_px_width:
                extra_in = (legend_bbox.x1 - fig_px_width) / fig.dpi
                width_in, height_in = fig.get_size_inches()
                fig.set_size_inches(width_in + extra_in + 0.25, height_in)
        plt.grid(True, which='both')

        # Save plot
        if self.save:
            try:
                plt.savefig(self.save_path)
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, self.show, self.save)

    def fisher_bounds(self):
        """
        Computes confidence bounds for the fitted distribution (self.dist).
        Called by mle() for bounds='fb' (and, for dist='exponential' only,
        also for bounds='chi2'). Weibull/Normal/LogNormal only have one
        bounds='fb' method: the shared Fisher-information/delta-method
        approach. Exponential has two genuinely different methods since it
        has a closed-form exact alternative that the other three don't:
        bounds='fb' is the same asymptotic Fisher-information/delta-method
        approach as the other three (see _fisher_bounds_exponential()'s
        docstring), while bounds='chi2' is an exact chi-square pivot (see
        _exact_bounds_exponential()'s docstring) - the two do NOT produce
        the same bounds, unlike an earlier version of this code where 'fb'
        was repointed at the chi-square method entirely.
        """
        if self.dist == 'normal':
            self._fisher_bounds_normal()
            return
        elif self.dist == 'lognormal':
            self._fisher_bounds_lognormal()
            return
        elif self.dist == 'exponential':
            if self.bounds == 'chi2':
                self._exact_bounds_exponential()
            else:
                self._fisher_bounds_exponential()
            return

        # Check if parameters are bias-corrected
        if self.bcm is None:
            b = self.beta
            eta = self.eta
        else:
            if self.bcm == 'c4':
                b = self.beta_c4
                eta = self.eta_c4
            elif self.bcm == 'hrbu':
                b = self.beta_hrbu
                eta = self.eta_hrbu
            elif self.bcm == 'np_bs':
                b = self.beta_np_bs
                eta = self.eta_np_bs
            elif self.bcm == 'p_bs':
                b = self.beta_p_bs
                eta = self.eta_p_bs
            else:
                raise ValueError('No valid bias-correction method is defined')

        # Fisher Information Matrix calculation
        if self.ds is None:
            element_1 = np.sum([(-1 / b ** 2) - (x / eta) ** b
                                * (np.log(x / eta) ** 2) for x in self.df])
            element_2 = np.sum([b / (eta ** 2) - ((x / eta) ** b)
                                * (b / (eta ** 2)) * (b + 1) for x in self.df])
            element_3 = np.sum([-1 / eta + ((x / eta) ** b)
                                * (1 / eta) * ((b * np.log(x / eta)) + 1) for x in self.df])
        else:  # does not work for interval censoring -> to do
            element_1 = (np.sum([(-1 / b ** 2) - (x / eta) ** b
                                * (np.log(x / eta) ** 2) for x in self.df])
                         + np.sum([-1 * (x / eta) ** b * np.log(x / eta) ** 2 for x in self.ds]))
            element_2 = (np.sum([b / (eta ** 2) - ((x / eta) ** b)
                                * (b / (eta ** 2)) * (b + 1) for x in self.df])
                         + np.sum([-1 * ((x / eta) ** b) * b / (eta ** 2)
                                   * (b + 1) for x in self.ds]))
            element_3 = (np.sum([-1 / eta + ((x / eta) ** b) * (1 / eta)
                                 * ((b * np.log(x / eta)) + 1) for x in self.df])
                         + np.sum([(x / eta) ** b / eta
                                   * (b * np.log(x / eta) + 1) for x in self.ds]))

        # f_inv equals cov_matrix, i.e. [var(beta), covar()],[covar(), var(eta)]
        self.f = np.matrix([[-element_1, -element_3], [-element_3, -element_2]])
        self.f_inv = np.linalg.inv(self.f)

        # K_a needed for bounds
        if self.bounds_type == '2s':
            self.k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
        elif self.bounds_type == '1su':
            # 1-sided upper
            self.k_a_bound = norm.ppf(self.cl)
        elif self.bounds_type == '1sl':
            # 1-sided lower
            self.k_a_bound = norm.ppf(1.0 - self.cl)

        # Standard error for each parameter
        self.se_beta = (self.f_inv.item(0)) ** 0.5
        self.se_eta = (self.f_inv.item(3)) ** 0.5

        k_a_p_bound = self.k_a_bound
        # qq.Caluclating t_lower and t_upper according to B. Dodson et. al
        var_y = (self.f_inv.item(3) / eta ** 2
                  + (np.log(-np.log(1 - np.array(self.unrel)))) ** 2
                  * self.f_inv.item(0) / b ** 4
                  - (2 * np.log(-np.log(1 - np.array(self.unrel)))
                    * self.f_inv.item(1) / ((b ** 2) * eta)))

        # Output bounds depending on k_a
        if self.bounds_type == '2s':
            y_l = (np.log(eta) + (np.log(-np.log(1 - np.array(self.unrel))) / b)
                   - k_a_p_bound * np.sqrt(var_y))
            self.bounds_lower = np.exp(y_l)

            y_u = (np.log(eta) + (np.log(-np.log(1 - np.array(self.unrel))) / b)
                   + k_a_p_bound * np.sqrt(var_y))
            self.bounds_upper = np.exp(y_u)
        elif self.bounds_type == '1su':
            # 1-sided upper
            y_u = (np.log(eta) + (np.log(-np.log(1 - np.array(self.unrel))) / b)
                   + k_a_p_bound * np.sqrt(var_y))
            self.bounds_upper = np.exp(y_u)
        elif self.bounds_type == '1sl':
            # 1-sided lower
            y_l = (np.log(eta) + (np.log(-np.log(1 - np.array(self.unrel))) / b)
                   + k_a_p_bound * np.sqrt(var_y))
            self.bounds_lower = np.exp(y_l)

    def _fisher_bounds_normal(self):
        """
        Fisher bounds for the Normal distribution. Same idea as
        fisher_bounds() above (invert the observed Fisher information at
        the MLE to get Var/Cov(mu_hat, sigma_hat), then propagate that to
        the time quantile t_p via the delta method), but simpler: since
        t_p = mu + sigma * z_p is already linear in (mu, sigma), the
        quantile variance is a plain quadratic form and no log/exp
        round-trip is needed the way Weibull needs it to keep eta positive.
        """
        mu, sigma = self.mu, self.sigma

        # Observed Fisher information = -Hessian of ln L at the MLE. Reuses
        # the same building blocks as mle()'s normal_score_jac.
        resid = np.asarray(self.df, dtype=float) - mu
        r = len(self.df)
        h_mumu = -r / sigma ** 2
        h_musigma = -2 * np.sum(resid) / sigma ** 3
        h_sigmasigma = r / sigma ** 2 - 3 * np.sum(resid ** 2) / sigma ** 4
        if self.ds:
            z = (np.asarray(self.ds, dtype=float) - mu) / sigma
            h = self._inv_mills_ratio(z)
            h_prime = h * (h - z)
            g = z * h
            g_prime = h + z * h_prime
            h_mumu += -np.sum(h * (h - z)) / sigma ** 2
            h_musigma += -np.sum(h) / sigma ** 2 - np.sum(z * h_prime) / sigma ** 2
            h_sigmasigma += -np.sum(g) / sigma ** 2 - np.sum(z * g_prime) / sigma ** 2

        self.f = np.matrix([[-h_mumu, -h_musigma], [-h_musigma, -h_sigmasigma]])
        self.f_inv = np.linalg.inv(self.f)

        if self.bounds_type == '2s':
            self.k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
        elif self.bounds_type == '1su':
            self.k_a_bound = norm.ppf(self.cl)
        elif self.bounds_type == '1sl':
            self.k_a_bound = norm.ppf(1.0 - self.cl)

        self.se_mu = (self.f_inv.item(0)) ** 0.5
        self.se_sigma = (self.f_inv.item(3)) ** 0.5
        k_a_p_bound = self.k_a_bound

        # Delta method: Var(t_p) = Var(mu) + z_p^2*Var(sigma) + 2*z_p*Cov(mu,sigma)
        z_p = norm.ppf(np.array(self.unrel))
        var_t = (self.f_inv.item(0) + z_p ** 2 * self.f_inv.item(3)
                 + 2 * z_p * self.f_inv.item(1))
        t_p = mu + sigma * z_p

        if self.bounds_type == '2s':
            self.bounds_lower = t_p - k_a_p_bound * np.sqrt(var_t)
            self.bounds_upper = t_p + k_a_p_bound * np.sqrt(var_t)
        elif self.bounds_type == '1su':
            self.bounds_upper = t_p + k_a_p_bound * np.sqrt(var_t)
        elif self.bounds_type == '1sl':
            # k_a_bound = norm.ppf(1 - cl) is negative for cl > 0.5, so
            # adding it (as fisher_bounds()'s Weibull branch also does for
            # 1sl) correctly pulls the bound below t_p.
            self.bounds_lower = t_p + k_a_p_bound * np.sqrt(var_t)

    def _fisher_bounds_lognormal(self):
        """
        Fisher bounds for the LogNormal distribution: identical to
        _fisher_bounds_normal() above but computed on ln(data), since mu/
        sigma are the Normal parameters of ln(T) (see mle()). The bounds
        are built on ln(t_p) - same reasoning as Weibull's log/exp
        round-trip - and only exponentiated back to the original time
        scale as the very last step, rather than doing the delta method
        directly in the (always-positive) time domain.
        """
        mu, sigma = self.mu, self.sigma
        log_df = np.log(self.df)
        log_ds = np.log(self.ds) if self.ds else None

        resid = log_df - mu
        r = len(self.df)
        h_mumu = -r / sigma ** 2
        h_musigma = -2 * np.sum(resid) / sigma ** 3
        h_sigmasigma = r / sigma ** 2 - 3 * np.sum(resid ** 2) / sigma ** 4
        if log_ds is not None:
            z = (log_ds - mu) / sigma
            h = self._inv_mills_ratio(z)
            h_prime = h * (h - z)
            g = z * h
            g_prime = h + z * h_prime
            h_mumu += -np.sum(h * (h - z)) / sigma ** 2
            h_musigma += -np.sum(h) / sigma ** 2 - np.sum(z * h_prime) / sigma ** 2
            h_sigmasigma += -np.sum(g) / sigma ** 2 - np.sum(z * g_prime) / sigma ** 2

        self.f = np.matrix([[-h_mumu, -h_musigma], [-h_musigma, -h_sigmasigma]])
        self.f_inv = np.linalg.inv(self.f)

        if self.bounds_type == '2s':
            self.k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
        elif self.bounds_type == '1su':
            self.k_a_bound = norm.ppf(self.cl)
        elif self.bounds_type == '1sl':
            self.k_a_bound = norm.ppf(1.0 - self.cl)

        self.se_mu = (self.f_inv.item(0)) ** 0.5
        self.se_sigma = (self.f_inv.item(3)) ** 0.5
        k_a_p_bound = self.k_a_bound

        z_p = norm.ppf(np.array(self.unrel))
        var_ln_t = (self.f_inv.item(0) + z_p ** 2 * self.f_inv.item(3)
                    + 2 * z_p * self.f_inv.item(1))
        ln_t_p = mu + sigma * z_p

        if self.bounds_type == '2s':
            self.bounds_lower = np.exp(ln_t_p - k_a_p_bound * np.sqrt(var_ln_t))
            self.bounds_upper = np.exp(ln_t_p + k_a_p_bound * np.sqrt(var_ln_t))
        elif self.bounds_type == '1su':
            self.bounds_upper = np.exp(ln_t_p + k_a_p_bound * np.sqrt(var_ln_t))
        elif self.bounds_type == '1sl':
            self.bounds_lower = np.exp(ln_t_p + k_a_p_bound * np.sqrt(var_ln_t))

    def _fisher_bounds_exponential(self):
        """
        Confidence bounds for the Exponential distribution's theta, used
        for bounds='fb'. Single parameter theta, so the observed Fisher
        information is a scalar: I(theta) = r/theta^2 (see mle()'s
        docstring for the log-likelihood), giving Var(theta_hat) =
        theta_hat^2/r and, by the delta method, Var(ln(theta_hat)) = 1/r
        exactly (no p-dependence, since ln(t_p) = ln(theta) +
        ln(-ln(1-p)) is additively shifted by a constant that doesn't
        depend on theta). Bounds are built in log space, like
        Weibull/LogNormal, to keep theta positive.

        This is only asymptotically correct (exact coverage as r ->
        infinity, like Weibull/Normal/LogNormal's own Fisher bounds) - see
        _exact_bounds_exponential() (bounds='chi2') for an exact
        alternative available for this one distribution specifically. The
        two are genuinely different methods and do not produce the same
        bounds.
        """
        r = len(self.df)
        self.se_theta = self.theta / np.sqrt(r)

        if self.bounds_type == '2s':
            self.k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
        elif self.bounds_type == '1su':
            self.k_a_bound = norm.ppf(self.cl)
        elif self.bounds_type == '1sl':
            self.k_a_bound = norm.ppf(1.0 - self.cl)
        k_a_p_bound = self.k_a_bound

        se_ln_theta = 1.0 / np.sqrt(r)
        ln_t_p = np.log(self.theta) + np.log(-np.log(1 - np.array(self.unrel)))

        if self.bounds_type == '2s':
            self.bounds_lower = np.exp(ln_t_p - k_a_p_bound * se_ln_theta)
            self.bounds_upper = np.exp(ln_t_p + k_a_p_bound * se_ln_theta)
        elif self.bounds_type == '1su':
            self.bounds_upper = np.exp(ln_t_p + k_a_p_bound * se_ln_theta)
        elif self.bounds_type == '1sl':
            self.bounds_lower = np.exp(ln_t_p + k_a_p_bound * se_ln_theta)

    def _exact_bounds_exponential(self):
        """
        Confidence bounds for the Exponential distribution's theta, used
        for bounds='chi2', via the exact chi-square pivot (Nelson, "Applied
        Life Data Analysis", 1982, p.255; Mathews, "Sample Size
        Calculations", 2010): for Type-II censored (or complete)
        exponential data,
            2 * r * theta_hat / theta ~ chi2(df=2r)
        exactly (r = number of failures) - not just asymptotically, unlike
        _fisher_bounds_exponential()'s Fisher-information/delta-method
        normal approximation (bounds='fb'), which only has the right
        coverage as r -> infinity. Rearranging the pivot for a two-sided
        100*cl% interval:
            theta_lower = 2*r*theta_hat / chi2.ppf(1 - alpha/2, df=2r)
            theta_upper = 2*r*theta_hat / chi2.ppf(alpha/2, df=2r)
        (alpha = 1 - cl; one-sided bounds replace alpha/2 by alpha, same as
        Nelson's convention). theta_lower/upper are still just a single CI
        on the one scalar parameter theta, so - same as bounds='fb' above -
        the resulting bounds on t_p = theta*(-ln(1-p)) are a constant
        multiple of the MLE line (parallel to it on Weibull's probability
        paper), not p-dependent the way Weibull's own beta/eta bounds are;
        that behavior is inherent to Exponential having only one free
        parameter, not a property of which CI method is used.

        self.se_theta (Var(theta_hat) = theta_hat^2/r, from the observed
        Fisher information I(theta) = r/theta^2) is set here too, purely as
        a standard-error summary - it plays no part in building these
        bounds.
        """
        r = len(self.df)
        self.se_theta = self.theta / np.sqrt(r)
        df_chi2 = 2 * r
        total_time = r * self.theta

        if self.bounds_type == '2s':
            alpha = 1.0 - self.cl
            theta_lower = 2 * total_time / chi2.ppf(1 - alpha / 2, df_chi2)
            theta_upper = 2 * total_time / chi2.ppf(alpha / 2, df_chi2)
        elif self.bounds_type == '1su':
            theta_upper = 2 * total_time / chi2.ppf(1 - self.cl, df_chi2)
        elif self.bounds_type == '1sl':
            theta_lower = 2 * total_time / chi2.ppf(self.cl, df_chi2)

        t_p_per_unit_theta = -np.log(1 - np.array(self.unrel))

        if self.bounds_type == '2s':
            self.bounds_lower = theta_lower * t_p_per_unit_theta
            self.bounds_upper = theta_upper * t_p_per_unit_theta
        elif self.bounds_type == '1su':
            self.bounds_upper = theta_upper * t_p_per_unit_theta
        elif self.bounds_type == '1sl':
            self.bounds_lower = theta_lower * t_p_per_unit_theta

    def _lrb_normal(self):
        """
        Likelihood-ratio confidence bounds for Normal/LogNormal, called by
        lrb()
        """
        fit_df = list(np.log(self.df)) if self.dist == 'lognormal' else list(self.df)
        fit_ds = None
        if self.ds:
            fit_ds = list(np.log(self.ds)) if self.dist == 'lognormal' else list(self.ds)

        def ll_func(mu_, sigma_):
            """Vectorized right-censored Normal log-likelihood (mu_, sigma_
            broadcastable arrays/mesh), reusing scipy's norm.logpdf/logsf
            for the same numerical robustness reasons mle() does."""
            ll = np.zeros(np.broadcast(mu_, sigma_).shape, dtype=np.float64)
            for x in fit_df:
                ll = ll + norm.logpdf(x, mu_, sigma_)
            if fit_ds:
                for x in fit_ds:
                    ll = ll + norm.logsf(x, mu_, sigma_)
            return ll

        def lr_z(mu_, sigma_):
            return 2 * (ll_func(mu_, sigma_) - self.ll_ref) + self.chi2_val

        def crossing_idx(z, axis):
            diff = np.diff(np.sign(z), axis=axis)
            diff = np.ma.array(diff, mask=np.isnan(diff))
            idx = np.argwhere(diff)
            return idx[:, 0], idx[:, 1]

        def collapse_to_envelope(fixed_pts, scan_pts, group_idx):
            if len(group_idx) == 0:
                return fixed_pts, scan_pts, group_idx
            order = np.argsort(group_idx, kind='stable')
            group_sorted = group_idx[order]
            fixed_sorted = fixed_pts[order]
            scan_sorted = scan_pts[order]
            out_fixed, out_scan, out_group = [], [], []
            start = 0
            n = len(group_sorted)
            for end in range(1, n + 1):
                if end == n or group_sorted[end] != group_sorted[start]:
                    block_scan = scan_sorted[start:end]
                    if block_scan.size > 2:
                        lo, hi = block_scan.min(), block_scan.max()
                        out_fixed.extend([fixed_sorted[start], fixed_sorted[start]])
                        out_scan.extend([lo, hi])
                        out_group.extend([group_sorted[start], group_sorted[start]])
                    else:
                        out_fixed.extend(fixed_sorted[start:end])
                        out_scan.extend(scan_sorted[start:end])
                        out_group.extend(group_sorted[start:end])
                    start = end
            return np.array(out_fixed), np.array(out_scan), np.array(out_group, dtype=int)

        def bisect_refine(fixed_vals, lo, hi, is_mu_fixed, n_iter=40):
            if is_mu_fixed:
                f_lo = lr_z(fixed_vals, lo)
            else:
                f_lo = lr_z(lo, fixed_vals)
            for _ in range(n_iter):
                mid = (lo + hi) / 2
                f_mid = lr_z(fixed_vals, mid) if is_mu_fixed else lr_z(mid, fixed_vals)
                keep_lo = (np.sign(f_mid) == np.sign(f_lo)) | (f_mid == 0)
                lo = np.where(keep_lo, mid, lo)
                hi = np.where(keep_lo, hi, mid)
                f_lo = np.where(keep_lo, f_mid, f_lo)
            return (lo + hi) / 2

        def scan_axis(z, mu_range, sigma_range, along):
            if along == 'sigma':
                row, col = crossing_idx(z, axis=1)
                mu_pts = mu_range[row]
                sigma_pts = bisect_refine(mu_pts, sigma_range[col], sigma_range[col + 1],
                                           is_mu_fixed=True)
                mu_pts, sigma_pts, group_idx = collapse_to_envelope(mu_pts, sigma_pts, row)
                n_groups = len(mu_range)
            else:
                row, col = crossing_idx(z, axis=0)
                sigma_pts = sigma_range[col]
                mu_pts = bisect_refine(sigma_pts, mu_range[row], mu_range[row + 1],
                                        is_mu_fixed=False)
                sigma_pts, mu_pts, group_idx = collapse_to_envelope(sigma_pts, mu_pts, col)
                n_groups = len(sigma_range)
            return mu_pts, sigma_pts, group_idx, n_groups

        def refine_lone_crossings(mu_pts, sigma_pts, group_idx, n_groups,
                                   fixed_range, other_range, is_mu_fixed,
                                   n_local=4000, expand=3.0):
            counts = np.bincount(group_idx, minlength=n_groups)
            lone = np.where(counts == 1)[0]
            if lone.size == 0:
                return np.array([]), np.array([])
            span = other_range.max() - other_range.min()
            lo = other_range.min() - expand * span
            hi = other_range.max() + expand * span
            if not is_mu_fixed:
                # other_range == sigma here (mu is the scanned axis's
                # fixed_range instead), so sigma must stay positive.
                lo = max(lo, 1e-9)
            local_vals = np.linspace(lo, hi, n_local)
            extra_mu, extra_sigma = [], []
            for i in lone:
                fixed_val = fixed_range[i]
                if is_mu_fixed:
                    mu_local, sigma_local = np.full(n_local, fixed_val), local_vals
                else:
                    mu_local, sigma_local = local_vals, np.full(n_local, fixed_val)
                z_local = lr_z(mu_local, sigma_local)
                diff_local = np.diff(np.sign(z_local))
                cross = np.argwhere(np.ma.array(diff_local, mask=np.isnan(diff_local))).ravel()
                if cross.size >= 2:
                    lo_c, hi_c = cross.min(), cross.max()
                    for c in (lo_c, hi_c):
                        mid = (local_vals[c] + local_vals[c + 1]) / 2
                        if is_mu_fixed:
                            extra_mu.append(fixed_val)
                            extra_sigma.append(mid)
                        else:
                            extra_mu.append(mid)
                            extra_sigma.append(fixed_val)
            return np.array(extra_mu), np.array(extra_sigma)

        def zerofinder(mu_range_init, sigma_range_init, z):
            row_s, col_s = crossing_idx(z, axis=1)
            row_m, col_m = crossing_idx(z, axis=0)
            mu_rough = np.concatenate([mu_range_init[row_s], mu_range_init[row_m]])
            sigma_rough = np.concatenate([sigma_range_init[col_s], sigma_range_init[col_m]])
            if mu_rough.size == 0:
                raise RuntimeError('lrb(): no solutions found on the initial mesh; '
                                    'the parameter range may need to be widened.')
            delta_mu = max(mu_rough.max() - mu_rough.min(), 1e-6)
            delta_sigma = max(sigma_rough.max() - sigma_rough.min(), 1e-6)
            mu_range = np.linspace(mu_rough.min() - delta_mu, mu_rough.max() + delta_mu, 1000)
            sigma_range = np.linspace(max(sigma_rough.min() - delta_sigma, 1e-9),
                                       sigma_rough.max() + delta_sigma, 1000)

            mm, ss = np.meshgrid(mu_range, sigma_range, indexing='ij')
            z_fine = lr_z(mm, ss)

            mu_sigma_pts, sigma_sigma_pts, group_sigma, n_mu = scan_axis(
                z_fine, mu_range, sigma_range, along='sigma')
            mu_mu_pts, sigma_mu_pts, group_mu, n_sigma = scan_axis(
                z_fine, mu_range, sigma_range, along='mu')

            extra_mu_1, extra_sigma_1 = refine_lone_crossings(
                mu_sigma_pts, sigma_sigma_pts, group_sigma, n_mu,
                fixed_range=mu_range, other_range=sigma_range, is_mu_fixed=True)
            extra_mu_2, extra_sigma_2 = refine_lone_crossings(
                mu_mu_pts, sigma_mu_pts, group_mu, n_sigma,
                fixed_range=sigma_range, other_range=mu_range, is_mu_fixed=False)

            mu_pairs = np.concatenate([mu_sigma_pts, mu_mu_pts, extra_mu_1, extra_mu_2])
            sigma_pairs = np.concatenate(
                [sigma_sigma_pts, sigma_mu_pts, extra_sigma_1, extra_sigma_2])
            return mu_pairs, sigma_pairs

        # Initial (mu, sigma) mesh range from the observed Fisher
        # information (same building blocks as _fisher_bounds_normal/
        # _fisher_bounds_lognormal) - only a starting point for the mesh,
        # not the final bounds themselves.
        mu, sigma = self.mu, self.sigma
        resid = np.asarray(fit_df, dtype=float) - mu
        r = len(fit_df)
        h_mumu = -r / sigma ** 2
        h_musigma = -2 * np.sum(resid) / sigma ** 3
        h_sigmasigma = r / sigma ** 2 - 3 * np.sum(resid ** 2) / sigma ** 4
        if fit_ds:
            z_ds = (np.asarray(fit_ds, dtype=float) - mu) / sigma
            h = self._inv_mills_ratio(z_ds)
            h_prime = h * (h - z_ds)
            g = z_ds * h
            g_prime = h + z_ds * h_prime
            h_mumu += -np.sum(h * (h - z_ds)) / sigma ** 2
            h_musigma += -np.sum(h) / sigma ** 2 - np.sum(z_ds * h_prime) / sigma ** 2
            h_sigmasigma += -np.sum(g) / sigma ** 2 - np.sum(z_ds * g_prime) / sigma ** 2

        f = np.array([[-h_mumu, -h_musigma], [-h_musigma, -h_sigmasigma]])
        f_inv = np.linalg.inv(f)
        se_mu = f_inv[0, 0] ** 0.5
        se_sigma = f_inv[1, 1] ** 0.5

        if self.bounds_type == '2s':
            k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
            self.cl_lrb = self.cl
        elif self.bounds_type == '1su':
            k_a_bound = norm.ppf(self.cl)
            self.cl_lrb = 2 * self.cl - 1
        elif self.bounds_type == '1sl':
            k_a_bound = norm.ppf(1.0 - self.cl)
            self.cl_lrb = 2 * self.cl - 1

        # Widen the Fisher-based interval (factor 2) so the initial mesh
        # comfortably contains the (generally wider, non-elliptical) LR
        # contour before zerofinder() refines around it.
        mu_lower = mu - 2 * abs(k_a_bound) * se_mu
        mu_upper = mu + 2 * abs(k_a_bound) * se_mu
        sigma_lower = max(sigma / np.exp(2 * abs(k_a_bound) * se_sigma / sigma), 1e-9)
        sigma_upper = sigma * np.exp(2 * abs(k_a_bound) * se_sigma / sigma)

        self.mu_range_init = np.linspace(mu_lower, mu_upper, 400)
        self.sigma_range_init = np.linspace(sigma_lower, sigma_upper, 400)

        mm_init, ss_init = np.meshgrid(self.mu_range_init, self.sigma_range_init, indexing='ij')

        self.chi2_val = chi2.ppf(self.cl_lrb, 1)
        self.ll_ref = ll_func(mu, sigma)

        with np.errstate(divide='ignore', invalid='ignore'):
            z_init = lr_z(mm_init, ss_init)
            self.mu_lrb, self.sigma_lrb = zerofinder(self.mu_range_init,
                                                      self.sigma_range_init,
                                                      z_init)

        z_p = norm.ppf(np.array(self.unrel))
        # t_p(mu, sigma) for every contour point x every unrel level.
        t_p_grid = self.mu_lrb[:, None] + self.sigma_lrb[:, None] * z_p[None, :]
        if self.dist == 'lognormal':
            t_p_grid = np.exp(t_p_grid)
        lower = t_p_grid.min(axis=0)
        upper = t_p_grid.max(axis=0)

        if self.bounds_type == '2s':
            self.bounds_lower, self.bounds_upper = lower, upper
        elif self.bounds_type == '1su':
            self.bounds_upper = upper
        elif self.bounds_type == '1sl':
            self.bounds_lower = lower

    def lrb(self):
        """
        # Goal: Find all solution pairs (beta, eta) for
        # L(beta, eta) = exp(chi ** 2) / -2) * L(beta_mle, eta_mle)
        """
        if self.dist in ('normal', 'lognormal'):
            self._lrb_normal()
            return

        def t_bounds_from_pars(beta_, eta, unreliability):
            """
            finds the minimum and maximum plausible time parameter for all
            given combination of "solutions" and "unreliability"
            """
            mins = np.zeros(len(unreliability))
            maxes = np.zeros_like(mins)
            for idx, unrel in enumerate(unreliability):
                ret = np.array(eta) * ((-np.log(1 - unrel)) ** (1 / np.array(beta_)))
                mins[idx] = min(ret)
                maxes[idx] = max(ret)
                self.mins = mins
                self.maxes = maxes
            return mins, maxes

        # Precompute the parts of the log-likelihood that do not depend
        # on (beta, eta) once, instead of re-deriving them on every
        # mesh/refinement evaluation. log(x/eta) is computed as
        # log(x) - log(eta) (log(x) precomputed, log(eta) computed once
        # per mesh) and x**beta_ as exp(log(x/eta) * beta_): both are
        # noticeably faster here than a fresh division/power per sample.
        log_df = np.log(self.df)
        sum_log_df = np.sum(log_df)
        n_df = len(self.df)
        if self.ds is not None:
            log_ds = np.log(self.ds)

        def ll_full(beta_, eta):
            """
            Vectorized right-censored Weibull log-likelihood (up to a
            (beta, eta)-independent constant, not needed for lrb).
            beta_, eta: broadcastable ndarrays (or plain floats).
            """
            log_eta = np.log(eta)
            f_sum = np.zeros_like(eta, dtype=np.float64)
            for lx in log_df:
                f_sum = f_sum + np.exp((lx - log_eta) * beta_)
            s_sum = np.zeros_like(eta, dtype=np.float64)
            for lx in log_ds:
                s_sum = s_sum + np.exp((lx - log_eta) * beta_)
            return (n_df * (np.log(beta_) - beta_ * log_eta)
                    + (beta_ - 1) * sum_log_df - f_sum - s_sum)

        def ll_full_no_cens(beta_, eta):
            """
            Vectorized uncensored Weibull log-likelihood (up to a
            (beta, eta)-independent constant, not needed for lrb).
            beta_, eta: broadcastable ndarrays (or plain floats).
            """
            log_eta = np.log(eta)
            f_sum = np.zeros_like(eta, dtype=np.float64)
            for lx in log_df:
                f_sum = f_sum + np.exp((lx - log_eta) * beta_)
            return (n_df * (np.log(beta_) - beta_ * log_eta)
                    + (beta_ - 1) * sum_log_df - f_sum)

        ll_func = ll_full_no_cens if self.ds is None else ll_full

        def lr_z(beta_, eta):
            """
            Likelihood-ratio test statistic shifted by the chi-squared
            critical value, so that lr_z == 0 marks the LRB contour.
            2 * log(exp(ll - ll_ref)) simplifies to 2 * (ll - ll_ref);
            skipping the exp/log round trip avoids needless overflow/
            underflow on large meshes and is also considerably faster.
            """
            return 2 * (ll_func(beta_, eta) - self.ll_ref) + self.chi2_val

        def crossing_idx(z, axis):
            """Row/col index pairs where z changes sign along `axis`."""
            diff = np.diff(np.sign(z), axis=axis)
            diff = np.ma.array(diff, mask=np.isnan(diff))
            idx = np.argwhere(diff)
            return idx[:, 0], idx[:, 1]

        def collapse_to_envelope(fixed_pts, scan_pts, group_idx):
            """
            Collapses groups (rows/columns, identified by group_idx)
            with more than two crossings down to just the two most
            extreme scan-axis values.

            Near a tangent point of the contour, z(scan_coord) can hug
            zero over a whole stretch; catastrophic cancellation in
            2 * (ll - ll_ref) then makes it flicker sign at the
            numerical-noise level instead of crossing cleanly, which
            without this collapse shows up as a "comb" of many points
            sharing the same fixed coordinate - rendered by
            contour_plot's convex hull as a spurious straight line
            segment. Collapsing to the two extremes is harmless for a
            genuine (convex) boundary too, since only the outermost
            points survive a convex hull anyway.
            """
            if len(group_idx) == 0:
                return fixed_pts, scan_pts, group_idx

            order = np.argsort(group_idx, kind='stable')
            group_sorted = group_idx[order]
            fixed_sorted = fixed_pts[order]
            scan_sorted = scan_pts[order]

            out_fixed, out_scan, out_group = [], [], []
            start = 0
            n = len(group_sorted)
            for end in range(1, n + 1):
                if end == n or group_sorted[end] != group_sorted[start]:
                    block_scan = scan_sorted[start:end]
                    if block_scan.size > 2:
                        lo, hi = block_scan.min(), block_scan.max()
                        out_fixed.extend([fixed_sorted[start], fixed_sorted[start]])
                        out_scan.extend([lo, hi])
                        out_group.extend([group_sorted[start], group_sorted[start]])
                    else:
                        out_fixed.extend(fixed_sorted[start:end])
                        out_scan.extend(scan_sorted[start:end])
                        out_group.extend(group_sorted[start:end])
                    start = end
            return np.array(out_fixed), np.array(out_scan), np.array(out_group, dtype=int)

        def bisect_refine(fixed_vals, lo, hi, is_beta_fixed, n_iter=40):
            """
            Vectorized bisection that pins every crossing down to near
            machine precision within its coarse mesh bracket [lo, hi],
            instead of just taking the bracket's midpoint.

            Without this, a near-vertical stretch of the contour (where
            many different fixed values share the same coarse bracket
            of the scanned coordinate) would have all of its solutions
            quantized to the same coordinate - visible as a spurious
            straight line segment once rendered. Refining every
            crossing individually resolves them to their own, slightly
            different position instead, tracing the true curve. All
            crossings for a direction are refined together in one
            vectorized pass (`fixed_vals`, `lo`, `hi` are arrays), so
            this stays cheap even for thousands of crossings.
            """
            if is_beta_fixed:
                f_lo = lr_z(fixed_vals, lo)
            else:
                f_lo = lr_z(lo, fixed_vals)
            for _ in range(n_iter):
                mid = (lo + hi) / 2
                f_mid = lr_z(fixed_vals, mid) if is_beta_fixed else lr_z(mid, fixed_vals)
                keep_lo = (np.sign(f_mid) == np.sign(f_lo)) | (f_mid == 0)
                lo = np.where(keep_lo, mid, lo)
                hi = np.where(keep_lo, hi, mid)
                f_lo = np.where(keep_lo, f_mid, f_lo)
            return (lo + hi) / 2

        def scan_axis(z, beta_range, eta_range, along):
            """
            Scans the mesh for solution pairs in one direction.

            along='eta': for every beta row, find the eta value(s) where
                z changes sign (the classic "fix beta, solve eta" scan).
            along='beta': for every eta column, find the beta value(s)
                where z changes sign ("fix eta, solve beta"), which is
                what recovers points near the left/right extremes of the
                contour that the eta-wise scan alone can miss.

            Returns beta_pts, eta_pts and the row/column index that each
            solution belongs to (used afterwards to spot rows/columns
            that only yielded a single solution).
            """
            if along == 'eta':
                row, col = crossing_idx(z, axis=1)
                beta_pts = beta_range[row]
                eta_pts = bisect_refine(beta_pts, eta_range[col], eta_range[col + 1],
                                         is_beta_fixed=True)
                beta_pts, eta_pts, group_idx = collapse_to_envelope(beta_pts, eta_pts, row)
                n_groups = len(beta_range)
            else:
                row, col = crossing_idx(z, axis=0)
                eta_pts = eta_range[col]
                beta_pts = bisect_refine(eta_pts, beta_range[row], beta_range[row + 1],
                                          is_beta_fixed=False)
                eta_pts, beta_pts, group_idx = collapse_to_envelope(eta_pts, beta_pts, col)
                n_groups = len(eta_range)
            return beta_pts, eta_pts, group_idx, n_groups

        def refine_lone_crossings(beta_pts, eta_pts, group_idx, n_groups,
                                   fixed_range, other_range, is_beta_fixed,
                                   n_local=4000, expand=3.0):
            """
            Workaround for rows/columns where only one crossing was
            found on the shared fine mesh (typically near-tangent
            points close to the extremes of the contour, where the
            missing second crossing lies outside the current mesh
            range). For every such row/column, redo a dense 1D scan
            over a widened range of the other coordinate to try to
            recover the missing solution. Only the two most extreme
            crossings found are kept (see collapse_to_envelope): a
            near-tangent point can make z flicker with many spurious
            noise-level sign changes rather than one clean crossing.
            """
            counts = np.bincount(group_idx, minlength=n_groups)
            lone = np.where(counts == 1)[0]
            if lone.size == 0:
                return np.array([]), np.array([])

            span = other_range.max() - other_range.min()
            lo = max(other_range.min() - expand * span, 1e-6)
            hi = other_range.max() + expand * span
            local_vals = np.linspace(lo, hi, n_local)

            extra_beta, extra_eta = [], []
            for i in lone:
                fixed_val = fixed_range[i]
                if is_beta_fixed:
                    b_local, e_local = np.full(n_local, fixed_val), local_vals
                else:
                    b_local, e_local = local_vals, np.full(n_local, fixed_val)

                z_local = lr_z(b_local, e_local)
                diff_local = np.diff(np.sign(z_local))
                cross = np.argwhere(np.ma.array(diff_local, mask=np.isnan(diff_local))).ravel()
                if cross.size >= 2:
                    lo_c, hi_c = cross.min(), cross.max()
                    for c in (lo_c, hi_c):
                        mid = (local_vals[c] + local_vals[c + 1]) / 2
                        if is_beta_fixed:
                            extra_beta.append(fixed_val)
                            extra_eta.append(mid)
                        else:
                            extra_beta.append(mid)
                            extra_eta.append(fixed_val)

            return np.array(extra_beta), np.array(extra_eta)

        def zerofinder(beta_range_init, eta_range_init, z):
            """
            Returns arrays with all solution pairs (beta, eta) for the
            LRB contour, found by scanning the mesh in both directions
            and locally refining rows/columns with only one crossing.
            """
            # Step 1: rough bounding box for the refined mesh, combining
            # both scan directions on the coarse mesh.
            row_e, col_e = crossing_idx(z, axis=1)
            row_b, col_b = crossing_idx(z, axis=0)

            beta_rough = np.concatenate([beta_range_init[row_e], beta_range_init[row_b]])
            eta_rough = np.concatenate([eta_range_init[col_e], eta_range_init[col_b]])

            if beta_rough.size == 0:
                raise RuntimeError('lrb(): no solutions found on the initial mesh; '
                                    'the parameter range may need to be widened.')

            # Step 2: one shared, finer mesh around the combined rough
            # bounding box, scanned in both directions. Sharing a single
            # mesh (instead of one per direction) halves the refine cost.
            # delta_* is a safety margin, e.g. for samples with a high
            # suspension ratio.
            delta_beta = max(beta_rough.max() - beta_rough.min(), 1e-6)
            delta_eta = max(eta_rough.max() - eta_rough.min(), 1e-6)
            beta_range = np.linspace(max(beta_rough.min() - delta_beta, 1e-6),
                                      beta_rough.max() + delta_beta, 1000)
            eta_range = np.linspace(max(eta_rough.min() - delta_eta, 1e-6),
                                     eta_rough.max() + delta_eta, 1000)

            bb, ee = np.meshgrid(beta_range, eta_range, indexing='ij')
            z_fine = lr_z(bb, ee)

            beta_eta_pts, eta_eta_pts, group_eta, n_beta = scan_axis(
                z_fine, beta_range, eta_range, along='eta')
            beta_beta_pts, eta_beta_pts, group_beta, n_eta = scan_axis(
                z_fine, beta_range, eta_range, along='beta')

            # Step 3: workaround for rows/columns with only one crossing.
            extra_beta_1, extra_eta_1 = refine_lone_crossings(
                beta_eta_pts, eta_eta_pts, group_eta, n_beta,
                fixed_range=beta_range, other_range=eta_range, is_beta_fixed=True)
            extra_beta_2, extra_eta_2 = refine_lone_crossings(
                beta_beta_pts, eta_beta_pts, group_beta, n_eta,
                fixed_range=eta_range, other_range=beta_range, is_beta_fixed=False)

            beta_pairs_lrb = np.concatenate(
                [beta_eta_pts, beta_beta_pts, extra_beta_1, extra_beta_2])
            eta_pairs_lrb = np.concatenate(
                [eta_eta_pts, eta_beta_pts, extra_eta_1, extra_eta_2])

            return (beta_pairs_lrb, eta_pairs_lrb)

        # 1.1 Calculate bounds for parameters using Fisher method

        # Check if parameters are bias-corrected
        if self.bcm is None:
            b = self.beta
            eta = self.eta
        else:
            if self.bcm == 'c4':
                b = self.beta_c4
                eta = self.eta
            elif self.bcm == 'hrbu':
                b = self.beta_hrbu
                eta = self.eta_hrbu
            elif self.bcm == 'np_bs':
                b = self.beta_np_bs
                eta = self.eta_np_bs
            elif self.bcm == 'p_bs':
                b = self.beta_p_bs
                eta = self.eta_p_bs
            else:
                raise ValueError('No valid bias-correction method is defined')


        # Hotfix to pass actual corrected beta and eta to the zerofinder
        self.sol_b = b
        self.sol_eta = eta

        # Compute elements F information matrix calculation
        if self.ds is None:
            element_1 = np.sum([(-1 / b ** 2) - (x / eta) ** b
                                * (np.log(x / eta) ** 2) for x in self.df])
            element_2 = np.sum([b / (eta ** 2) - ((x / eta) ** b)
                                * (b / (eta ** 2)) * (b + 1) for x in self.df])
            element_3 = np.sum([-1 / eta + ((x / eta) ** b)
                                * (1 / eta) * ((b * np.log(x / eta)) + 1) for x in self.df])
        else:  # does not work for interval censoring -> to do
            element_1 = (np.sum([(-1 / b ** 2) - (x / eta) ** b
                                * (np.log(x / eta) ** 2) for x in self.df])
                         + np.sum([-1 * (x / eta) ** b * np.log(x / eta) ** 2 for x in self.ds]))
            element_2 = (np.sum([b / (eta ** 2) - ((x / eta) ** b)
                                * (b / (eta ** 2)) * (b + 1) for x in self.df])
                         + np.sum([-1 * ((x / eta) ** b) * b / (eta ** 2)
                                   * (b + 1) for x in self.ds]))
            element_3 = (np.sum([-1 / eta + ((x / eta) ** b) * (1 / eta)
                                 * ((b * np.log(x / eta)) + 1) for x in self.df])
                         + np.sum([(x / eta) ** b / eta
                                   * (b * np.log(x / eta) + 1) for x in self.ds]))

        # f_inv equals cov_matrix, i.e. [var(beta), covar()],[covar(), var(eta)]
        self.f = np.matrix([[-element_1, -element_3], [-element_3, -element_2]])
        self.f_inv = np.linalg.inv(self.f)

        # K_a needed for bounds
        # Setting k_a_bounds manually is sufficient just to get a paramter range
        if self.bounds_type == '2s':
            self.k_a_bound = norm.ppf((1.0 - self.cl) / 2 + self.cl)
            # No need to adapt self.cl for 2-sided bounds
            self.cl_lrb = self.cl
        elif self.bounds_type == '1su':
            # 1-sided upper
            self.k_a_bound = norm.ppf(self.cl)
            self.cl_lrb = 2 * self.cl - 1
        elif self.bounds_type == '1sl':
            # 1-sided lower
            self.k_a_bound = norm.ppf(1.0 - self.cl)
            self.cl_lrb = 2 * self.cl - 1
        else:
            print('break')

        # Calculate inital parameter bounds using f_inv
        beta_lower = b / (np.exp(self.k_a_bound * np.sqrt(self.f_inv.item(0)) / b))
        beta_upper = b * (np.exp(self.k_a_bound * np.sqrt(self.f_inv.item(0)) / b))
        eta_lower = eta / (np.exp(self.k_a_bound * np.sqrt(self.f_inv.item(3)) / eta))
        eta_upper = eta * (np.exp(self.k_a_bound * np.sqrt(self.f_inv.item(3)) / eta))

        self.beta_f_range = [beta_lower, beta_upper]
        self.eta_f_range = [eta_lower, eta_upper]


        # 1.2 Create mesh with repect to the paramter range
        # Beta bounds are critical to the mesh resolution, hence two steps to produce beta range

        self.beta_range_init = np.arange(.2 * self.beta_f_range[0], 2 * self.beta_f_range[1], 0.02)
        self.eta_range_init = np.linspace(.5 * self.eta_f_range[0], 2 * self.eta_f_range[1], 1000)

        # Create mesh
        bb, ee = np.meshgrid(self.beta_range_init, self.eta_range_init, indexing='ij')

        # Precompute the reference log-likelihood and chi-squared
        # critical value once, both reused for every mesh/refinement
        # evaluation inside lr_z().
        self.chi2_val = chi2.ppf(self.cl_lrb, 1)
        self.ll_ref = ll_func(b, eta)

        # Ignore log(-inf) since this is not relevant
        with np.errstate(divide='ignore', invalid='ignore'):
            self.z = lr_z(bb, ee)
            self.beta_lrb, self.eta_lrb = zerofinder(self.beta_range_init,
                                                     self.eta_range_init,
                                                     self.z)

        # Calculate Solutions with Zerofinder
        self.bounds_lower, self.bounds_upper = t_bounds_from_pars(self.beta_lrb,
                                                                  self.eta_lrb,
                                                                  self.unrel)
    def plot(self):
        """
        Creates a probability plot for the distribution set via dist= in
        the constructor.
        """
        if self.dist == 'normal':
            return self._plot_normal()
        elif self.dist == 'lognormal':
            return self._plot_lognormal()
        elif self.dist == 'exponential':
            return self._plot_exponential()

        # Some needed functions:
        def weibull_prob_paper(x):
            """
            Needed to adjust figure to the Weibull probability plot.
            """
            x = np.asarray(x)

            # Prevent np.log(0) error raise
            x[x > .9999] = np.nan
            return np.log(-np.log(1 - x))

        # Just for y_tickslabel on the y-axis
        def weibull_ticks(y_i, _):
            """
            Adjusts the y-axis tick labels
            """
            return '{:.1f}'.format((100 * (1 - np.exp(-np.exp(y_i)))))

        def unrel_func(x_est, beta_, eta):
            if type(x_est) == list:
                x_est = np.asarray(x_est)
            y_est = (1 - np.exp(-(x_est / eta) ** beta_))
            y_est_lnln = weibull_prob_paper(y_est)

            return y_est_lnln

        def inverse_weibull(perc, beta, eta):
            """
            Computes time to failure data points.
            This function is being used to plot Weibull lines.

            Parameters
            ----------
            perc : float
                Percentage points fo which time to failure data should be computed.
            beta : float
                Weibull shape parameter.
            eta : float
                Weibull scale parameter.

            Returns
            -------
            float
                Time to failure data points.

            """
            return ((-np.log(1 -perc)) ** (1 / beta)) * eta

        # Generate Weibull Plot Figure
        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= self.y_min)
                                          & (PROBABILITY_PLOT_TICKS <= self.y_max)]
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.yticks(lny_ticks, color='black')
        ax.set_yticks([weibull_prob_paper(0.632)], minor=True)

        # Plots the horizontal dashed line for 63.2%
        plt.grid(True, which='minor', axis='y', linestyle='--')

        # X-Axis scaling
        if self.bcm == 'c4':
            if self.bounds is None:
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_c4)) * self.eta
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_c4)) * self.eta
            else:
                if self.bounds_type == '2s':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1su':
                    self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_c4)) * self.eta
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1sl':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_c4)) * self.eta
        elif self.bcm == 'hrbu':
            if self.bounds is None:
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_hrbu)) * self.eta_hrbu
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_hrbu)) * self.eta_hrbu
            else:
                if self.bounds_type == '2s':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1su':
                    self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_hrbu)) * self.eta_hrbu
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1sl':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_hrbu)) * self.eta_hrbu
        elif self.bcm == 'np_bs':
            if self.bounds is None:
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_np_bs)) * self.eta_np_bs
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_np_bs)) * self.eta_np_bs
            else:
                if self.bounds_type == '2s':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1su':
                    self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_np_bs)) * self.eta_np_bs
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1sl':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_np_bs)) * self.eta_np_bs
        elif self.bcm == 'p_bs':
            if self.bounds is None:
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_p_bs)) * self.eta_p_bs
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_p_bs)) * self.eta_p_bs
            else:
                if self.bounds_type == '2s':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1su':
                    self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta_p_bs)) * self.eta_p_bs
                    self.tmax_plot = max(self.bounds_upper)
                if self.bounds_type == '1sl':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta_p_bs)) * self.eta_p_bs
        elif self.bcm is None:
            if self.bounds is None:
                self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
                self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta
            else:
                if self.bounds_type == '2s':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = max(self.bounds_upper)
                elif self.bounds_type == '1su':
                    self.tmin_plot = ((-np.log(0.999)) ** (1 / self.beta)) * self.eta
                    self.tmax_plot = max(self.bounds_upper)
                elif self.bounds_type == '1sl':
                    self.tmin_plot = min(self.bounds_lower)
                    self.tmax_plot = ((-np.log(0.001)) ** (1 / self.beta)) * self.eta

        self.xplot = np.linspace(self.tmin_plot, self.tmax_plot, 100)
        left = (10 ** (np.ceil(np.log10(self.tmin_plot)) - 1))
        right = (10 ** (np.ceil(np.log10(self.tmax_plot))))
        plt.xlim(left, right)
        plt.tick_params(axis='x', colors='black')

        # Set labels and legends
        plt.title(self.plot_title, color='black', fontsize=self.plot_title_fontsize)
        plt.xlabel(f'{self.x_label}{" in "+self.unit if self.unit!="-" else ""}', color='black', fontsize=self.xy_fontsize)
        plt.xticks(fontsize=self.tick_fontsize)
        plt.ylabel(f'{self.y_label} in %', color='black', fontsize=self.xy_fontsize)
        plt.yticks(fontsize=self.tick_fontsize)

        # General style properties
        self.unrel = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01,
                               0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4,
                               0.5, 0.6, 0.632, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])

        # Plot legend
        if self.ds is None:
            susp_num = 0
        else:
            susp_num = len(self.ds)

        # Check if bias-corrections are applied
        if self.bcm == 'c4':
            # Plot corrected line
            xvals_c4 = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta_c4, self.eta_c4))
            plt.semilogx(xvals_c4, unrel_func(xvals_c4, self.beta_c4,self.eta_c4),
                         color=PREDICTR_FIT_COLOR, linestyle='-',
                         linewidth=1.5, zorder=2)

            # Plot biased line
            xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
            plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                         color='grey', linestyle='--',
                         linewidth=1.5, zorder=1)


            # Define title in legend
            leg_title = 'MLE C4'

            if self.bounds is not None:
                # Adapt bounds' name for legend, if bounds are applied
                if self.bounds == 'fb':
                    bounds_legend = 'Fisher bounds'
                elif self.bounds == 'lrb':
                    bounds_legend = 'LRB'
                elif self.bounds == 'npbb':
                    bounds_legend = 'Non-Par. Boostrap bounds'
                elif self.bounds == 'pbb':
                    bounds_legend = 'Par. Bootstrap bounds'

                # Plot legend
                if self.bounds_type == '2s':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_c4)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_c4),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1su':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_c4)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_c4),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\1su @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1sl':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_c4)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_c4),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
            else:
                if self.show_legend:
                    plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                  + susp_num,
                                                                  len(self.df),
                                                                  susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta_c4)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta),
                                '\nuncorrected MLE:\n'
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$'.format(self.eta)),
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=leg_title)
        elif self.bcm == 'hrbu':
            # Plot corrected line
            xvals_hrbu = list(inverse_weibull(np.array([0.001, 0.9999]),
                                              self.beta_hrbu, self.eta_hrbu))
            plt.semilogx(xvals_hrbu, unrel_func(xvals_hrbu, self.beta_hrbu,self.eta_hrbu),
                         color=PREDICTR_FIT_COLOR, linestyle='-',
                         linewidth=1.5, zorder=2)

            # Plot biased line
            xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
            plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                         color='grey', linestyle='--',
                         linewidth=1.5, zorder=1)

            # Define title in legend
            leg_title = 'MLE HRBU'

            if self.bounds is not None:
                # Adapt bounds' name for legend, if bounds are applied
                if self.bounds == 'fb':
                    bounds_legend = 'Fisher bounds'
                elif self.bounds == 'lrb':
                    bounds_legend = 'LRB'
                elif self.bounds == 'npbb':
                    bounds_legend = 'Non-Par. Boostrap bounds'
                elif self.bounds == 'pbb':
                    bounds_legend = 'Par. Bootstrap bounds'

                # Plot legend
                if self.bounds_type == '2s':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR,
                                 linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_hrbu)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_hrbu),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1su':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_hrbu)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_hrbu),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1su @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1sl':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_hrbu)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta_hrbu),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
            else:
                if self.show_legend:
                    plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                  + susp_num,
                                                                  len(self.df),
                                                                  susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta_hrbu)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta_hrbu),
                                '\nuncorrected MLE:\n'
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$'.format(self.eta)),
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=leg_title)
        elif self.bcm == 'np_bs':
            # Plot corrected line
            xvals_np_bs = list(inverse_weibull(np.array([0.001, 0.9999]),
                                              self.beta_np_bs, self.eta_np_bs))
            plt.semilogx(xvals_np_bs, unrel_func(xvals_np_bs, self.beta_np_bs,self.eta_np_bs),
                         color=PREDICTR_FIT_COLOR, linestyle='-',
                         linewidth=1.5, zorder=2)

            # Plot biased line
            xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
            plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                         color='grey', linestyle='--',
                         linewidth=1.5, zorder=1)

            # Define title in legend
            leg_title = 'MLE n.-p. Bootstrap'

            if self.bounds is not None:
                # Adapt bounds' name for legend, if bounds are applied
                if self.bounds == 'fb':
                    bounds_legend = 'Fisher bounds'
                elif self.bounds == 'lrb':
                    bounds_legend = 'LRB'
                elif self.bounds == 'npbb':
                    bounds_legend = 'Non-Par. Boostrap bounds'
                elif self.bounds == 'pbb':
                    bounds_legend = 'Par. Bootstrap bounds'

                # Plot legend
                if self.bounds_type == '2s':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_np_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_np_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1su':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_np_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_np_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\1su @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1sl':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_np_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_np_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
            else:
                if self.show_legend:
                    plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                  + susp_num,
                                                                  len(self.df),
                                                                  susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta_np_bs)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta_np_bs)
                                + '\nstatistic: {}'.format(self.est_type)
                                +'\nBS samples: {}'.format(self.bs_size),
                                '\nuncorrected MLE:\n' + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$'.format(self.eta)),
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=leg_title)
        elif self.bcm == 'p_bs':
            # Plot corrected line
            xvals_p_bs = list(inverse_weibull(np.array([0.001, 0.9999]),
                                              self.beta_p_bs, self.eta_p_bs))
            plt.semilogx(xvals_p_bs, unrel_func(xvals_p_bs, self.beta_p_bs,self.eta_p_bs),
                         color=PREDICTR_FIT_COLOR, linestyle='-',
                         linewidth=1.5, zorder=2)

            # Plot biased line
            xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
            plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                         color='grey', linestyle='--',
                         linewidth=1.5, zorder=1)

            # Define title in legend
            leg_title = 'MLE par. Bootstrap'

            if self.bounds is not None:
                # Adapt bounds' name for legend, if bounds are applied
                if self.bounds == 'fb':
                    bounds_legend = 'Fisher bounds'
                elif self.bounds == 'lrb':
                    bounds_legend = 'LRB'
                elif self.bounds == 'npbb':
                    bounds_legend = 'Non-Par. Boostrap bounds'
                elif self.bounds == 'pbb':
                    bounds_legend = 'Par. Bootstrap bounds'

                # Plot legend
                if self.bounds_type == '2s':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_p_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_p_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1su':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_p_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_p_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1su @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1sl':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta_p_bs)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta_p_bs)
                                    + '\nstatistic: {}'.format(self.est_type)
                                    +'\nBS samples: {}'.format(self.bs_size),
                                    '\nuncorrected MLE:\n'
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$'.format(self.eta),
                                    '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
            else:
                if self.show_legend:
                    plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                  + susp_num,
                                                                  len(self.df),
                                                                  susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta_p_bs)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta_p_bs)
                                + '\nstatistic: {}'.format(self.est_type)
                                +'\nBS samples: {}'.format(self.bs_size),
                                '\nuncorrected MLE:\n'
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$'.format(self.eta)),
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=leg_title)
        else:
            # Plot biased line
            xvals = list(inverse_weibull(np.array([0.001, 0.9999]), self.beta, self.eta))
            plt.semilogx(xvals, unrel_func(xvals, self.beta,self.eta),
                         color=PREDICTR_FIT_COLOR, linestyle='-',
                         linewidth=1.5)

            # Define title in legend
            leg_title = 'MLE'

            if self.bounds is not None:
                # Adapt bounds' name for legend, if bounds are applied
                if self.bounds == 'fb':
                    bounds_legend = 'Fisher bounds'
                elif self.bounds == 'lrb':
                    bounds_legend = 'LRB'
                elif self.bounds == 'npbb':
                    bounds_legend = 'Non-Par. Boostrap bounds'
                elif self.bounds == 'pbb':
                    bounds_legend = 'Par. Bootstrap bounds'

                # 2-sided bounds
                if self.bounds_type == '2s':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = PREDICTR_FIT_COLOR, label='_nolegend_')
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta),
                                    '\n{}:\n2s @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                elif self.bounds_type == '1su':
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta),
                                    '\n{}:\n1su @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
                # 1-sided lower bounds
                elif self.bounds_type == '1sl':
                    plt.semilogx(self.bounds_lower, weibull_prob_paper(self.unrel),
                                 color=PREDICTR_FIT_COLOR, linestyle='-', linewidth=1)
                    if self.show_legend:
                        plt.legend(('n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                      + susp_num,
                                                                      len(self.df),
                                                                      susp_num)
                                    + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                    + r'$\widehat\eta={:.3f}$ '.format(self.eta),
                                    '\n{}:\n1sl @{}%'.format((bounds_legend), self.cl * 100)),
                                   loc='lower left', bbox_to_anchor=(0.65, 0.0),
                                   fontsize=self.legend_fontsize, title=leg_title)
            else:
                if self.show_legend:
                    plt.legend(['n = {} (f: {} | s: {})\n'.format(len(self.df)
                                                                  + susp_num,
                                                                  len(self.df),
                                                                  susp_num)
                                + r'$\widehat\beta={:.3f}$ | '.format(self.beta)
                                + r'$\widehat\eta={:.3f}$ '.format(self.eta)],
                               loc='lower left', bbox_to_anchor=(0.65, 0.0),
                               fontsize=self.legend_fontsize, title=leg_title)
        if self.plot_ranks:
            # Plot median ranks
            if self.ds is None:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank()), marker='o',
                             markerfacecolor=PREDICTR_FIT_COLOR, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)
            else:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank_cens()), marker='o',
                             markerfacecolor=PREDICTR_FIT_COLOR, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder= 3)

        # Pin the y-axis to [y_min, y_max] as the very last step: fill_between/
        # fill_betweenx above auto-expand the view when they're added
        # (matplotlib's add_collection autoscales regardless of an earlier
        # plt.ylim() call), so setting this any earlier would get silently
        # overridden once the bounds shading is drawn.
        ax.set_ylim(bottom=weibull_prob_paper(self.y_min),
                    top=weibull_prob_paper(self.y_max))
        plt.tight_layout()

        # bbox_to_anchor=(0.65, 0.0) above can push the legend past the
        # right edge of the figure depending on bounds/bounds_type; widen
        # the figure so the legend stays fully visible instead of being
        # clipped. Measured via a throwaway Agg renderer so this works
        # reliably regardless of the active display backend (interactive
        # backends like macosx don't reliably populate a renderer via a
        # plain canvas.draw()). Must run AFTER tight_layout(): tight_layout
        # tightens the axes margins, shifting the (axes-relative) legend
        # anchor further right - measuring before it would just get
        # squeezed away again. Deliberately not calling tight_layout()
        # again afterwards, so the added width stays as real margin instead
        # of being re-absorbed.
        legend = plt.gca().get_legend()
        if legend is not None:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            original_canvas = fig.canvas
            FigureCanvasAgg(fig)
            fig.canvas.draw()
            legend_bbox = legend.get_window_extent(fig.canvas.get_renderer())
            fig.canvas = original_canvas
            fig_px_width = fig.get_size_inches()[0] * fig.dpi
            if legend_bbox.x1 > fig_px_width:
                extra_in = (legend_bbox.x1 - fig_px_width) / fig.dpi
                width_in, height_in = fig.get_size_inches()
                fig.set_size_inches(width_in + extra_in + 0.25, height_in)
        plt.grid(True, which='both')

        # Save plot
        if self.save:
            try:
                plt.savefig(self.save_path)
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, self.show, self.save)

    def _plot_normal(self):
        """
        Creates the Normal probability plot. Uses predictr's shared plot
        style/kwargs (fig_size, fonts, legend, save/show, ...) like the
        Weibull plot() above, but with axis scaling appropriate for the
        Normal distribution instead of Weibull's: the x-axis (time) stays
        linear rather than log-scaled, since the Normal quantile function
        t_p = mu + sigma*z_p is already linear - there's nothing to
        linearize by taking logs the way Weibull's plot needs to. The
        y-axis is transformed via the standard normal quantile function
        (probit) instead of Weibull's ln(-ln(1-F)).
        """
        def normal_ticks(y_i, _):
            """
            Converts a z-value on the transformed y-axis back to the
            unreliability percentage it represents, for the tick labels.
            """
            return '{:.1f}'.format(100 * norm.cdf(y_i))

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))

        # Y-axis: probit scale, same tick percentages as the Weibull plot
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= self.y_min)
                                          & (PROBABILITY_PLOT_TICKS <= self.y_max)]
        z_ticks = norm.ppf(y_ticks)
        plt.yticks(z_ticks, color='black')
        ax.set_yticks([0.0], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        # PREDICTR_PALETTE (the categorical 6-color palette) is for telling
        # multiple *different* Analysis objects apart, e.g. in PlotAll -
        # not for the elements of a single result. A single fit instead
        # uses predictr's own single-result convention, exactly like
        # Weibull's plot() above: PREDICTR_FIT_COLOR for the MLE line, the
        # bounds (thinner, plus a light fill) and the marker faces.
        x_data = np.array(self.df, dtype=float)
        susp_num = len(self.ds) if self.ds is not None else 0

        # Fitted MLE line. Evaluated across a fixed, wide percentile range
        # (matching Weibull's inverse_weibull(np.array([0.001, 0.9999]), ...)
        # convention) rather than padded around the data's x-range, so the
        # line always reaches the edges of whatever [y_min, y_max] window is
        # actually visible instead of stopping short of it. Only 2 points
        # are needed since the line is exactly straight on this axis.
        z_line = norm.ppf(np.array([0.0001, 0.9999]))
        x_line = self.mu + self.sigma * z_line
        plt.plot(x_line, z_line, color=PREDICTR_FIT_COLOR, linestyle='-',
                 linewidth=1.5, zorder=2)

        leg_title = 'MLE'
        leg_text = ('n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                        len(self.df), susp_num)
                    + r'$\widehat\mu={:.3f}$ | '.format(self.mu)
                    + r'$\widehat\sigma={:.3f}$'.format(self.sigma))
        legend_labels = (leg_text,)

        if self.bounds in ('fb', 'lrb') and (self.bounds_lower is not None
                                              or self.bounds_upper is not None):
            bounds_legend = 'Fisher bounds' if self.bounds == 'fb' else 'LRB'
            z_p = norm.ppf(self.unrel)
            if self.bounds_type == '2s':
                plt.plot(self.bounds_lower, z_p, color=PREDICTR_FIT_COLOR,
                         linestyle='-', linewidth=1)
                plt.plot(self.bounds_upper, z_p, color=PREDICTR_FIT_COLOR,
                         linestyle='-', linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=z_p, x1=self.bounds_lower, x2=self.bounds_upper,
                                   alpha=0.1, color=PREDICTR_FIT_COLOR, label='_nolegend_')
                bt_legend = '2s'
            elif self.bounds_type == '1su':
                plt.plot(self.bounds_upper, z_p, color=PREDICTR_FIT_COLOR,
                         linestyle='-', linewidth=1)
                bt_legend = '1su'
            elif self.bounds_type == '1sl':
                plt.plot(self.bounds_lower, z_p, color=PREDICTR_FIT_COLOR,
                         linestyle='-', linewidth=1)
                bt_legend = '1sl'
            legend_labels = (leg_text, '\n{}:\n{} @{}%'.format(
                bounds_legend, bt_legend, self.cl * 100))

        plt.xlabel(f'{self.x_label}{" in " + self.unit if self.unit != "-" else ""}',
                   color='black', fontsize=self.xy_fontsize)
        plt.ylabel(self.y_label + ' in %', color='black', fontsize=self.xy_fontsize)
        plt.title(self.plot_title, color='black', fontsize=self.plot_title_fontsize)
        plt.tick_params(labelsize=self.tick_fontsize)
        plt.grid(True, which='major')
        if self.show_legend:
            plt.legend(legend_labels, loc='lower left', bbox_to_anchor=(0.65, 0.0),
                       fontsize=self.legend_fontsize, title=leg_title)

        # Data points via median ranks (plotted after legend(), like
        # Weibull's plot_ranks block, so they're never picked up as an
        # extra legend entry) - identical machinery to the Weibull plot,
        # since plotting positions don't depend on the assumed
        # distribution, only the axis transform applied to them afterwards.
        if self.plot_ranks:
            ranks = np.array(self.median_rank() if self.ds is None
                              else self.median_rank_cens())
            y_data = norm.ppf(ranks)
            plt.plot(x_data, y_data, marker='o', markerfacecolor=PREDICTR_FIT_COLOR,
                     markeredgecolor='black', markersize=4, alpha=.5,
                     linestyle='None', zorder=3)

        # Pin the y-axis to [y_min, y_max] as the very last step - see the
        # identical comment in plot()/plot_mrr() for why fill_betweenx above
        # would otherwise silently re-expand the view.
        ax.set_ylim(bottom=norm.ppf(self.y_min), top=norm.ppf(self.y_max))

        plt.tight_layout()

        # Same legend-overflow-safety as plot()/plot_mrr() above: widen the
        # figure if the legend box would otherwise get clipped at the right
        # edge. See the comment on the first occurrence of this pattern
        # (plot()) for why it must run after tight_layout().
        legend = plt.gca().get_legend()
        if legend is not None:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            original_canvas = fig.canvas
            FigureCanvasAgg(fig)
            fig.canvas.draw()
            legend_bbox = legend.get_window_extent(fig.canvas.get_renderer())
            fig.canvas = original_canvas
            fig_px_width = fig.get_size_inches()[0] * fig.dpi
            if legend_bbox.x1 > fig_px_width:
                extra_in = (legend_bbox.x1 - fig_px_width) / fig.dpi
                width_in, height_in = fig.get_size_inches()
                fig.set_size_inches(width_in + extra_in + 0.25, height_in)

        if self.save:
            try:
                plt.savefig(self.save_path)
            except:
                raise ValueError('Path is faulty.')
        return _finish_plot(fig, self.show, self.save)

    def _plot_lognormal(self):
        """
        Creates the LogNormal probability plot. Combines both of the other
        two axis conventions at once: the y-axis is probit-transformed
        exactly like _plot_normal() (mu/sigma are the Normal parameters of
        ln(T), so the same z = Phi^-1(F) linearizes it), but the x-axis is
        log-scaled like Weibull's plot() - because it's ln(T), not T
        itself, that's linear in z. That combination is exactly why
        LogNormal needs its own probability paper rather than reusing
        either of the other two as-is.
        """
        def normal_ticks(y_i, _):
            return '{:.1f}'.format(100 * norm.cdf(y_i))

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))

        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= self.y_min)
                                          & (PROBABILITY_PLOT_TICKS <= self.y_max)]
        z_ticks = norm.ppf(y_ticks)
        plt.yticks(z_ticks, color='black')
        ax.set_yticks([0.0], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        x_data = np.array(self.df, dtype=float)
        susp_num = len(self.ds) if self.ds is not None else 0

        # Fitted MLE line, evaluated across a fixed wide percentile range
        # like _plot_normal() - see that method's comment for why.
        z_line = norm.ppf(np.array([0.0001, 0.9999]))
        x_line = np.exp(self.mu + self.sigma * z_line)
        plt.semilogx(x_line, z_line, color=PREDICTR_FIT_COLOR, linestyle='-',
                     linewidth=1.5, zorder=2)

        leg_title = 'MLE'
        leg_text = ('n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                        len(self.df), susp_num)
                    + r'$\widehat\mu={:.3f}$ | '.format(self.mu)
                    + r'$\widehat\sigma={:.3f}$'.format(self.sigma))
        legend_labels = (leg_text,)

        if self.bounds in ('fb', 'lrb') and (self.bounds_lower is not None
                                              or self.bounds_upper is not None):
            bounds_legend = 'Fisher bounds' if self.bounds == 'fb' else 'LRB'
            z_p = norm.ppf(self.unrel)
            if self.bounds_type == '2s':
                plt.semilogx(self.bounds_lower, z_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                plt.semilogx(self.bounds_upper, z_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=z_p, x1=self.bounds_lower, x2=self.bounds_upper,
                                   alpha=0.1, color=PREDICTR_FIT_COLOR, label='_nolegend_')
                bt_legend = '2s'
            elif self.bounds_type == '1su':
                plt.semilogx(self.bounds_upper, z_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                bt_legend = '1su'
            elif self.bounds_type == '1sl':
                plt.semilogx(self.bounds_lower, z_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                bt_legend = '1sl'
            legend_labels = (leg_text, '\n{}:\n{} @{}%'.format(
                bounds_legend, bt_legend, self.cl * 100))

        plt.xlabel(f'{self.x_label}{" in " + self.unit if self.unit != "-" else ""}',
                   color='black', fontsize=self.xy_fontsize)
        plt.ylabel(self.y_label + ' in %', color='black', fontsize=self.xy_fontsize)
        plt.title(self.plot_title, color='black', fontsize=self.plot_title_fontsize)
        plt.tick_params(labelsize=self.tick_fontsize)
        plt.grid(True, which='major')
        if self.show_legend:
            plt.legend(legend_labels, loc='lower left', bbox_to_anchor=(0.65, 0.0),
                       fontsize=self.legend_fontsize, title=leg_title)

        if self.plot_ranks:
            ranks = np.array(self.median_rank() if self.ds is None
                              else self.median_rank_cens())
            y_data = norm.ppf(ranks)
            plt.semilogx(x_data, y_data, marker='o', markerfacecolor=PREDICTR_FIT_COLOR,
                         markeredgecolor='black', markersize=4, alpha=.5,
                         linestyle='None', zorder=3)

        ax.set_ylim(bottom=norm.ppf(self.y_min), top=norm.ppf(self.y_max))

        # Pin the x-axis to a sensible data-driven range too (mirroring
        # Weibull's plot()): otherwise autoscale would size it to the MLE
        # line's own evaluation points (z=[norm.ppf(0.0001), ...]), which on
        # a log axis can extend towards 0 far past anything meaningful.
        if self.bounds in ('fb', 'lrb') and (self.bounds_lower is not None
                                              or self.bounds_upper is not None):
            if self.bounds_type == '2s':
                tmin_plot, tmax_plot = min(self.bounds_lower), max(self.bounds_upper)
            elif self.bounds_type == '1su':
                tmin_plot = np.exp(self.mu + self.sigma * norm.ppf(0.001))
                tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1sl':
                tmin_plot = min(self.bounds_lower)
                tmax_plot = np.exp(self.mu + self.sigma * norm.ppf(0.999))
        else:
            tmin_plot = np.exp(self.mu + self.sigma * norm.ppf(0.001))
            tmax_plot = np.exp(self.mu + self.sigma * norm.ppf(0.999))
        ax.set_xlim(10 ** (np.ceil(np.log10(tmin_plot)) - 1),
                    10 ** (np.ceil(np.log10(tmax_plot))))

        plt.tight_layout()

        legend = plt.gca().get_legend()
        if legend is not None:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            original_canvas = fig.canvas
            FigureCanvasAgg(fig)
            fig.canvas.draw()
            legend_bbox = legend.get_window_extent(fig.canvas.get_renderer())
            fig.canvas = original_canvas
            fig_px_width = fig.get_size_inches()[0] * fig.dpi
            if legend_bbox.x1 > fig_px_width:
                extra_in = (legend_bbox.x1 - fig_px_width) / fig.dpi
                width_in, height_in = fig.get_size_inches()
                fig.set_size_inches(width_in + extra_in + 0.25, height_in)

        if self.save:
            try:
                plt.savefig(self.save_path)
            except:
                raise ValueError('Path is faulty.')
        return _finish_plot(fig, self.show, self.save)

    def _plot_exponential(self):
        """
        Creates the Exponential probability plot on the same log-x /
        double-log-y paper as Weibull's plot() above, rather than the
        simpler "native" exponential paper (linear x, single-log
        y = -ln(1-F)) that's also a mathematically valid linearization.
        Both work because t_p = theta*(-ln(1-p)) is exactly Weibull's
        quantile function with beta=1, eta=theta: Weibull's
        ln(t_p) = ln(eta) + (1/beta)*ln(-ln(1-p)) becomes, at beta=1,
        ln(t_p) = ln(theta) + ln(-ln(1-p)) - a straight (unit-slope) line
        on Weibull's own paper too. Plotting it this way, treated as the
        beta=1 special case of Weibull, avoids the low-percentile
        compression the simpler paper has (its stretch factor 1/(1-F) is
        ~1 near F=0 but diverges towards F=1, whereas this double-log
        transform stretches both tails).
        """
        def weibull_prob_paper(x):
            x = np.asarray(x, dtype=float)
            x = np.where(x > .9999, np.nan, x)
            return np.log(-np.log(1 - x))

        def weibull_ticks(y_i, _):
            return '{:.1f}'.format(100 * (1 - np.exp(-np.exp(y_i))))

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))

        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= self.y_min)
                                          & (PROBABILITY_PLOT_TICKS <= self.y_max)]
        lny_ticks = weibull_prob_paper(y_ticks)
        plt.yticks(lny_ticks, color='black')
        ax.set_yticks([weibull_prob_paper(0.632)], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        x_data = np.array(self.df, dtype=float)
        susp_num = len(self.ds) if self.ds is not None else 0

        # Fitted MLE line: t_p = theta*(-ln(1-p)) is Weibull's quantile
        # function at beta=1, eta=theta - 2 points suffice since it's a
        # straight line on this paper, evaluated across a fixed wide
        # percentile range like the other _plot_*() methods.
        p_line = np.array([0.0001, 0.9999])
        x_line = self.theta * (-np.log(1 - p_line))
        plt.semilogx(x_line, weibull_prob_paper(p_line), color=PREDICTR_FIT_COLOR,
                     linestyle='-', linewidth=1.5, zorder=2)

        leg_title = 'MLE'
        leg_text = ('n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                        len(self.df), susp_num)
                    + r'$\widehat\theta={:.3f}$'.format(self.theta))
        legend_labels = (leg_text,)

        if self.bounds in ('fb', 'chi2') and (self.bounds_lower is not None
                                     or self.bounds_upper is not None):
            y_p = weibull_prob_paper(self.unrel)
            if self.bounds_type == '2s':
                plt.semilogx(self.bounds_lower, y_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                plt.semilogx(self.bounds_upper, y_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=y_p, x1=self.bounds_lower, x2=self.bounds_upper,
                                   alpha=0.1, color=PREDICTR_FIT_COLOR, label='_nolegend_')
                bt_legend = '2s'
            elif self.bounds_type == '1su':
                plt.semilogx(self.bounds_upper, y_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                bt_legend = '1su'
            elif self.bounds_type == '1sl':
                plt.semilogx(self.bounds_lower, y_p, color=PREDICTR_FIT_COLOR,
                             linestyle='-', linewidth=1)
                bt_legend = '1sl'
            # Exponential is the only dist with two distinct bounds
            # methods (see fisher_bounds()'s docstring), so the legend
            # names whichever one was actually used instead of always
            # saying "Fisher bounds" like the other three plots' legends.
            bounds_method_legend = ('Exact bounds (χ²)' if self.bounds == 'chi2'
                                     else 'Fisher bounds')
            legend_labels = (leg_text, '\n{}:\n{} @{}%'.format(
                bounds_method_legend, bt_legend, self.cl * 100))

        plt.xlabel(f'{self.x_label}{" in " + self.unit if self.unit != "-" else ""}',
                   color='black', fontsize=self.xy_fontsize)
        plt.ylabel(self.y_label + ' in %', color='black', fontsize=self.xy_fontsize)
        plt.title(self.plot_title, color='black', fontsize=self.plot_title_fontsize)
        plt.tick_params(labelsize=self.tick_fontsize)
        plt.grid(True, which='major')
        if self.show_legend:
            plt.legend(legend_labels, loc='lower left', bbox_to_anchor=(0.65, 0.0),
                       fontsize=self.legend_fontsize, title=leg_title)

        if self.plot_ranks:
            ranks = np.array(self.median_rank() if self.ds is None
                              else self.median_rank_cens())
            y_data = weibull_prob_paper(ranks)
            plt.semilogx(x_data, y_data, marker='o', markerfacecolor=PREDICTR_FIT_COLOR,
                         markeredgecolor='black', markersize=4, alpha=.5,
                         linestyle='None', zorder=3)

        ax.set_ylim(bottom=weibull_prob_paper(self.y_min), top=weibull_prob_paper(self.y_max))

        # Pin the x-axis to a sensible data-driven range too (mirroring
        # Weibull's plot()): otherwise autoscale would size it to the MLE
        # line's own evaluation points (p=[0.0001, ...]), which on a log
        # axis can extend towards 0 far past anything meaningful.
        if self.bounds in ('fb', 'chi2') and (self.bounds_lower is not None
                                     or self.bounds_upper is not None):
            if self.bounds_type == '2s':
                tmin_plot, tmax_plot = min(self.bounds_lower), max(self.bounds_upper)
            elif self.bounds_type == '1su':
                tmin_plot = self.theta * (-np.log(1 - 0.001))
                tmax_plot = max(self.bounds_upper)
            elif self.bounds_type == '1sl':
                tmin_plot = min(self.bounds_lower)
                tmax_plot = self.theta * (-np.log(1 - 0.999))
        else:
            tmin_plot = self.theta * (-np.log(1 - 0.001))
            tmax_plot = self.theta * (-np.log(1 - 0.999))
        ax.set_xlim(10 ** (np.ceil(np.log10(tmin_plot)) - 1),
                    10 ** (np.ceil(np.log10(tmax_plot))))

        plt.tight_layout()

        legend = plt.gca().get_legend()
        if legend is not None:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            original_canvas = fig.canvas
            FigureCanvasAgg(fig)
            fig.canvas.draw()
            legend_bbox = legend.get_window_extent(fig.canvas.get_renderer())
            fig.canvas = original_canvas
            fig_px_width = fig.get_size_inches()[0] * fig.dpi
            if legend_bbox.x1 > fig_px_width:
                extra_in = (legend_bbox.x1 - fig_px_width) / fig.dpi
                width_in, height_in = fig.get_size_inches()
                fig.set_size_inches(width_in + extra_in + 0.25, height_in)

        if self.save:
            try:
                plt.savefig(self.save_path)
            except:
                raise ValueError('Path is faulty.')
        return _finish_plot(fig, self.show, self.save)

    @classmethod
    def get_bx_percentile(cls, time, beta_, eta_):
        """
        Computes the unreliability at given input time.

        Parameters
        ----------
        time : float or list of floats
            Lifetime for which the percentiles are computed. If time is a list, the percentiles for
            each element of the list will be computed and returned.
        beta_ : float
            Weibull shape parameter.
        eta_ : float
            Weibull scale parameter.
        eta_ : float
            Weibull scale parameter.
        Returns
        -------
        unrel : float or list of floats
            Percentiles for the given BXlife.

        """

        # Weibull function
        def unrel_func(time, beta_, eta_):
            unrel = (1 - np.exp(-(time / eta_) ** beta_))
            return unrel

        # Check if bx is of type: list
        if isinstance(time, list):
            percentiles = [unrel_func(i, beta_, eta_) for i in time]
        else:
            percentiles = unrel_func(time, beta_, eta_)

        return percentiles

class PlotAll:
    """
    Plots cdfs, pdfs and Weibull plots for multiple instances
    """

    def __init__(self, objects=None, **kwargs):
        if objects is not None:
            self.objects = objects
            for key, val in objects.items():
                self.unit = getattr(val, 'unit')
                self.plot_style = getattr(val, 'plot_style')
            self.unrel = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01 , 0.02 , 0.03 , 0.05 ,
                                   0.07 , 0.1  , 0.2  , 0.3  , 0.4  , 0.5  , 0.6  , 0.632, 0.7,
                                   0.8  , 0.9  , 0.95 , 0.99 , 0.999])
            # if len(self.objects.keys()) > 6:
            #     raise ValueError('mult_weibull only support up to six instances being plotted.')

        # # Set colormap for Weibull plot
        # if 'set_cmap' in kwargs:
        #     self.color = iter(kwargs['set_cmap'])
        # else:
        #     self.color = iter(['royalblue', 'salmon', 'mediumseagreen',
        #                        'darkorange', 'peru', 'darkcyan'])
    def median_rank(self, cl=0.5):
        """
        Mediran ranks for uncensored data. Returns a list with
        median ranks.
        """
        ranks = []
        n = len(self.df)
        for i in range(1, n+1):
            ranks.append(beta.ppf(cl, i, n-i+1))
        return ranks

    def median_rank_cens(self, cl=0.5):
        """
        Returns adjusted median ranks as described in the
        Weibull Handbook. Returns a list with adjusted median ranks.
        """

        def bernard(adj_r, n, cl):
            """
            Returns Bernards Approximation for the adjusted ranks
            """
            #return (np.array(i) - 0.3) / (len(self.df+self.ds) + 0.4)
            return [beta.ppf(cl, i, n-i+1) for i in adj_r]

        n = len(self.df + self.ds)
        # Reverse ranks need to consider suspensions and their order
        all_ = self.df + self.ds
        rev_rank = []
        prev = 0
        for j in self.df:
            # Check if failure time is entered multiple times
            if self.df.count(j) > 1:
                # Ignore same elements after first time
                if prev == j:
                    pass
                else:
                    # Number of times element is in df
                    count_element = self.df.count(j)
                    # Loop through identical failure time
                    for i in range(count_element):
                        count = sum(map(lambda x : x < j, all_)) + i
                        rev_rank.append(len(all_) - count)
                prev = j
            else:
                count = sum(map(lambda x : x < j, all_))
                rev_rank.append(len(all_) - count)

        #Calculate adjusted rank
        adj_ranks = []
        prev_rank = 0
        for i in range(1, len(self.df)+1):
            adj_ranks.append((rev_rank[i-1] * prev_rank + n + 1) / (rev_rank[i-1] +1))
            prev_rank = adj_ranks[-1]
        self.adj_ranks = bernard(adj_ranks, n, cl)
        return self.adj_ranks



    def mult_weibull(self, x_label='Time To Failure', y_label='Unreliability',
                     plot_title='Weibull Probability Plot', xy_fontsize=12,
                     plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7),
                     x_bounds=None, plot_ranks=True, save=False, color=None, linestyle=None,
                     y_min=0.01, y_max=0.99, show=True,
                     **kwargs):
        """
        Plots multiple Analysis class objects in one figure

        Parameters
        ----------
        y_min : float, optional
            Lower y-axis limit (unreliability, as a fraction) shown on the
            plot. Must satisfy 0 < y_min < y_max < 1. The default is 0.01,
            matching Analysis's default - see Analysis(y_min=, y_max=).
        y_max : float, optional
            Upper y-axis limit (unreliability, as a fraction) shown on the
            plot. Must satisfy 0 < y_min < y_max < 1. The default is 0.99.
        """
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        non_weibull = [key for key, val in self.objects.items()
                       if getattr(val, 'dist', 'weibull') != 'weibull']
        if non_weibull:
            raise ValueError(f'mult_weibull() only supports Analysis objects '
                              f'with dist="weibull", but {non_weibull} do not.')

        def inverse_weibull(perc, beta, eta):
            return ((-np.log(1 -perc)) ** (1 / beta)) * eta

        def weibull_prob_paper(x):
            """
            Needed to adjust figure to the Weibull probability plot.
            """
            x = np.asarray(x)

            # Prevent np.log(0) error raise
            x[x > .9999] = np.nan
            return np.log(-np.log(1 - x))

        # Just for y_tickslabel on the y-axis
        def weibull_ticks(y_i, _):
            """
            Adjusts the y-axis tick labels
            """
            return '{:.1f}'.format((100 * (1 - np.exp(-np.exp(y_i)))))

        def unrel_func(x_est, beta_, eta):
            if type(x_est) == list:
                x_est = np.asarray(x_est)
            y_est = (1 - np.exp(-(x_est / eta) ** beta_))
            y_est_lnln = weibull_prob_paper(y_est)
            return y_est_lnln

        def median_rank(df, cl=0.5):
            """
            Mediran ranks for uncensored data. Returns a list with
            median ranks.
            """
            ranks = []
            n = len(df)
            for i in range(1, n+1):
                ranks.append(beta.ppf(cl, i, n-i+1))
            return ranks

        def median_rank_cens(df, ds, cl=0.5):
            """
            Returns adjusted median ranks as described in the
            Weibull Handbook. Returns a list with adjusted median ranks.
            """

            def bernard(adj_r, n, cl):
                """
                Returns Bernards Approximation for the adjusted ranks
                """
                #return (np.array(i) - 0.3) / (len(self.df+self.ds) + 0.4)
                return [beta.ppf(cl, i, n-i+1) for i in adj_r]

            n = len(df + ds)
            # Reverse ranks need to consider suspensions and their order
            all_ = df + ds
            rev_rank = []
            prev = 0
            for j in df:
                # Check if failure time is entered multiple times
                if df.count(j) > 1:
                    # Ignore same elements after first time
                    if prev == j:
                        pass
                    else:
                        # Number of times element is in df
                        count_element = df.count(j)
                        # Loop through identical failure time
                        for i in range(count_element):
                            count = sum(map(lambda x : x < j, all_)) + i
                            rev_rank.append(len(all_) - count)
                    prev = j
                else:
                    count = sum(map(lambda x : x < j, all_))
                    rev_rank.append(len(all_) - count)

            #Calculate adjusted rank
            adj_ranks = []
            prev_rank = 0
            for i in range(1, len(df)+1):
                adj_ranks.append((rev_rank[i-1] * prev_rank + n + 1) / (rev_rank[i-1] +1))
                prev_rank = adj_ranks[-1]
            adj_ranks = bernard(adj_ranks, n, cl)
            return adj_ranks

        # Set line color/linestyle. When neither is given, predictr's
        # categorical palette is used and, for more than 6 objects, the
        # linestyle is cycled alongside the color so each dataset stays
        # visually distinguishable instead of repeating a bare color (see
        # _categorical_style()).
        num_objects = len(self.objects)
        if color is None and linestyle is None:
            default_colors, default_linestyles = _categorical_style(num_objects)
            color = iter(default_colors)
            l_style = iter(default_linestyles)
        else:
            if color is not None:
                color = iter(color)
            else:
                default_colors, _ = _categorical_style(num_objects)
                color = iter(default_colors)

            # Check linestyle input
            if linestyle is not None:
                # Check if provided number of linstyle is in accordance with number of objects
                if num_objects != len(linestyle):
                    raise ValueError(f'Number of linestyles ({len(linestyle)}) is'\
                                     ' not in accordance with the number of'\
                                         f' objects ({num_objects}).')
                l_style = iter(linestyle)
            else:
                l_style = iter(num_objects * ['-'])

        # Get t_min and t_max to plot
        temp_list = []
        for key, val in self.objects.items():
            # Check for special case of percentil bounds, e.g. 'bbb'
            if (getattr(val, 'bounds')) != 'bbb':
                if (getattr(val, 'bounds_lower')) is not None:
                    temp_list.append(min(getattr(val, 'bounds_lower')))
                if (getattr(val, 'bounds_upper')) is not None:
                    temp_list.append(max(getattr(val, 'bounds_upper')))
            if ((getattr(val, 'bounds_lower')) is None
                and (getattr(val, 'bounds_upper')) is None):
                if getattr(val, 'beta_c4') is not None:
                    dat = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_c4'),
                                               getattr(val, 'eta_c4')))
                    temp_list.append(min(dat))
                    temp_list.append(max(dat))
                elif getattr(val, 'beta_hrbu') is not None:
                    dat = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_hrbu'),
                                               getattr(val, 'eta_hrbu')))
                    temp_list.append(min(dat))
                    temp_list.append(max(dat))
                elif getattr(val, 'beta_np_bs') is not None:
                    dat = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_np_bs'),
                                               getattr(val, 'eta_np_bs')))
                    temp_list.append(min(dat))
                    temp_list.append(max(dat))
                elif getattr(val, 'beta_p_bs') is not None:
                    dat = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_p_bs'),
                                               getattr(val, 'eta_p_bs')))
                    temp_list.append(min(dat))
                    temp_list.append(max(dat))
                else:
                    dat = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta'),
                                               getattr(val, 'eta')))
                    temp_list.append(min(dat))
                    temp_list.append(max(dat))
        
        # Check for custom x-axis limits 
        if x_bounds is None:
            x_axis_min = min(temp_list)
            x_axis_max = max(temp_list)
        else:
            # Check if x_bounds is of type list
            if type(x_bounds) != list:
                raise TypeError(f'x_bounds need to be of type list and not {type(x_bounds)}.') 
            x_axis_min = x_bounds[0]
            x_axis_max = x_bounds[1]

        # Generate Weibull Plot Figure
        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=fig_size, num=next(_figure_counter))

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= y_min)
                                          & (PROBABILITY_PLOT_TICKS <= y_max)]
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.yticks(lny_ticks, color='black')
        ax.set_yticks([weibull_prob_paper(0.632)], minor=True)

        # Plots the horizontal dashed line for 63.2%
        plt.grid(True, which='minor', axis='y', linestyle='--')
        #xplot = np.linspace(x_axis_min, x_axis_max, 2000)
        left = (10 ** (np.ceil(np.log10(x_axis_min)) - 1))
        right = (10 ** (np.ceil(np.log10(x_axis_max))))
        plt.xlim(left, right)
        plt.tick_params(axis='x', colors='black')

        # Set labels and legends
        plt.title(plot_title, color='black', fontsize=plot_title_fontsize)
        
        plt.xlabel(f'{x_label}{" in "+self.unit if self.unit!="-" else ""}', color='black', fontsize=xy_fontsize)
        plt.ylabel(f'{y_label} in %', color='black', fontsize=xy_fontsize)

        # Plot Weibull lines
        for key, val in self.objects.items():
            if getattr(val, 'bounds') is None:
                col = next(color)
                ls = next(l_style)
                if getattr(val, 'beta_c4') is not None:
                    xvals = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_c4'),
                                               getattr(val, 'eta_c4')))
                    xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xplot, unrel_func(xplot,
                                                   getattr(val, 'beta_c4'),
                                                   getattr(val, 'eta_c4')),
                                 color=col, linestyle=ls,
                                 linewidth=1.5,label=f'{key}')
                elif getattr(val, 'beta_hrbu') is not None:
                    xvals = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_hrbu'),
                                               getattr(val, 'eta_hrbu')))
                    xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xplot, unrel_func(xplot,
                                                   getattr(val, 'beta_hrbu'),
                                                   getattr(val, 'eta_hrbu')),
                                 color=col, linestyle=ls,
                                 linewidth=1.5, label=f'{key}')
                elif getattr(val, 'beta_np_bs') is not None:
                    xvals = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_np_bs'),
                                               getattr(val, 'eta_np_bs')))
                    xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xplot, unrel_func(xplot, getattr(val, 'beta_np_bs'),
                                                   getattr(val, 'eta_np_bs')),
                                 color=col, linestyle=ls,
                                 linewidth=1.5, label=f'{key}')
                elif getattr(val, 'beta_p_bs') is not None:
                    xvals = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_p_bs'),
                                               getattr(val, 'eta_p_bs')))
                    xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xplot, unrel_func(xplot, getattr(val, 'beta_p_bs'),
                                                   getattr(val, 'eta_p_bs')),
                                 color=col, linestyle=ls,
                                 linewidth=1.5, label=f'{key}')
                else:
                    xvals = list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta'),
                                               getattr(val, 'eta')))
                    xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xplot, unrel_func(xplot, getattr(val, 'beta'),
                                                   getattr(val, 'eta')),
                                 color=col, linestyle=ls,
                                 linewidth=1.5, label=f'{key}')
            if getattr(val, 'bounds_type') == '2s':
                if getattr(val, 'beta_c4') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_c4'),
                                               getattr(val, 'eta_c4'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=getattr(val, 'bounds_lower'),
                                     x2=getattr(val, 'bounds_upper'),
                                     alpha=0.1, color = col, label='_nolegend_')
                elif getattr(val, 'beta_hrbu') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_hrbu'),
                                               getattr(val, 'eta_hrbu'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=getattr(val, 'bounds_lower'),
                                     x2=getattr(val, 'bounds_upper'),
                                     alpha=0.1, color = col, label='_nolegend_')
                elif getattr(val, 'beta_np_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_np_bs'),
                                               getattr(val, 'eta_np_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=getattr(val, 'bounds_lower'),
                                     x2=getattr(val, 'bounds_upper'),
                                     alpha=0.1, color = col, label='_nolegend_')
                elif getattr(val, 'beta_p_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_p_bs'),
                                               getattr(val, 'eta_p_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=getattr(val, 'bounds_lower'),
                                     x2=getattr(val, 'bounds_upper'),
                                     alpha=0.1, color = col, label='_nolegend_')
                elif (getattr(val, 'bounds')) == 'bbb':
                    pass
                else:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta'),
                                               getattr(val, 'eta'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=getattr(val, 'bounds_lower'),
                                     x2=getattr(val, 'bounds_upper'),
                                     alpha=0.1, color = col, label='_nolegend_')
            if getattr(val, 'bounds_type') == '1su':
                if getattr(val, 'beta_c4') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_c4'),
                                               getattr(val, 'eta_c4'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_hrbu') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_hrbu'),
                                               getattr(val, 'eta_hrbu'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_np_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_np_bs'),
                                               getattr(val, 'eta_np_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_p_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_p_bs'),
                                               getattr(val, 'eta_p_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                else:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta'),
                                               getattr(val, 'eta'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_upper'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
            if getattr(val, 'bounds_type') == '1sl':
                if getattr(val, 'beta_c4') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_c4'),
                                               getattr(val, 'eta_c4'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_hrbu') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_hrbu'),
                                               getattr(val, 'eta_hrbu'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_np_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_np_bs'),
                                               getattr(val, 'eta_np_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                elif getattr(val, 'beta_p_bs') is not None:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta_p_bs'),
                                               getattr(val, 'eta_p_bs'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
                else:
                    col = next(color)
                    ls = next(l_style)
                    xvals = sorted(list(inverse_weibull(np.array([0.001, 0.9999]),
                                               getattr(val, 'beta'),
                                               getattr(val, 'eta'))))

                    #xplot = np.linspace(min(xvals), max(xvals), 100)
                    plt.semilogx(xvals, np.log(-np.log(1 - np.array([0.001, 0.9999]))),
                                 color=col, linestyle=ls, linewidth=1.5, zorder = 2,
                                 label=f'{key}')
                    plt.semilogx(getattr(val, 'bounds_lower'), weibull_prob_paper(self.unrel),
                                 color=col, linestyle=ls, linewidth=1, label='_nolegend_')
            

            # Plot discrete median ranks
            if plot_ranks:
                if getattr(val, 'ds') is None:
                    plt.semilogx(getattr(val, 'df'),
                                 weibull_prob_paper(median_rank(getattr(val, 'df'))),
                                 marker='o',
                                 markerfacecolor=col, markeredgecolor='black',
                                 markersize=4, alpha=.5, linestyle='None', zorder=3)
                else:
                    plt.semilogx(getattr(val, 'df'),
                                 weibull_prob_paper(median_rank_cens(getattr(val, 'df'), getattr(val, 'ds'))),
                                 marker='o',
                                 markerfacecolor=col, markeredgecolor='black',
                                 markersize=4, alpha=.5, linestyle='None', zorder= 3)

        # Pin the y-axis to [y_min, y_max] as the very last step - see the
        # identical comment in Analysis.plot()/plot_mrr() for why the
        # fill_betweenx calls above would otherwise silently re-expand it.
        ax.set_ylim(bottom=weibull_prob_paper(y_min), top=weibull_prob_paper(y_max))

        plt.tight_layout()
        plt.grid(True, which='both')
        plt.legend(fontsize= legend_fontsize)

        # Save plot
        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, show, save)

    @staticmethod
    def _mult_median_rank(df, cl=0.5):
        """
        Median ranks for uncensored data, as a standalone helper (mirrors
        Analysis.median_rank()) so mult_normal/mult_lognormal/
        mult_exponential don't depend on self.df/self.ds.
        """
        n = len(df)
        return [beta.ppf(cl, i, n - i + 1) for i in range(1, n + 1)]

    @staticmethod
    def _mult_median_rank_cens(df, ds, cl=0.5):
        """
        Adjusted median ranks for censored data, as a standalone helper
        (mirrors Analysis.median_rank_cens()).
        """
        def bernard(adj_r, n, cl):
            return [beta.ppf(cl, i, n - i + 1) for i in adj_r]

        n = len(df + ds)
        all_ = df + ds
        rev_rank = []
        prev = 0
        for j in df:
            if df.count(j) > 1:
                if prev == j:
                    pass
                else:
                    count_element = df.count(j)
                    for i in range(count_element):
                        count = sum(map(lambda x: x < j, all_)) + i
                        rev_rank.append(len(all_) - count)
                prev = j
            else:
                count = sum(map(lambda x: x < j, all_))
                rev_rank.append(len(all_) - count)

        adj_ranks = []
        prev_rank = 0
        for i in range(1, len(df) + 1):
            adj_ranks.append((rev_rank[i - 1] * prev_rank + n + 1) / (rev_rank[i - 1] + 1))
            prev_rank = adj_ranks[-1]
        return bernard(adj_ranks, n, cl)

    def _mult_style(self, color, linestyle):
        """
        Shared color/linestyle iterator setup used by mult_weibull() and
        its normal/lognormal/exponential counterparts below.
        """
        num_objects = len(self.objects)
        if color is None and linestyle is None:
            default_colors, default_linestyles = _categorical_style(num_objects)
            color_it = iter(default_colors)
            l_style = iter(default_linestyles)
        else:
            if color is not None:
                color_it = iter(color)
            else:
                default_colors, _ = _categorical_style(num_objects)
                color_it = iter(default_colors)

            if linestyle is not None:
                if num_objects != len(linestyle):
                    raise ValueError(f'Number of linestyles ({len(linestyle)}) is'
                                      ' not in accordance with the number of'
                                      f' objects ({num_objects}).')
                l_style = iter(linestyle)
            else:
                l_style = iter(num_objects * ['-'])
        return color_it, l_style

    def mult_normal(self, x_label='Time To Failure', y_label='Unreliability',
                     plot_title='Normal Probability Plot', xy_fontsize=12,
                     plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7),
                     x_bounds=None, plot_ranks=True, save=False, color=None,
                     linestyle=None, y_min=0.01, y_max=0.99, show=True,
                     **kwargs):
        """
        Plots multiple Analysis class objects with dist='normal' in one
        figure, analogous to mult_weibull(). The x-axis stays linear (the
        Normal quantile function is already linear in z), unlike
        mult_weibull()/mult_lognormal()/mult_exponential() which use a
        log x-axis - see Analysis._plot_normal()'s docstring.

        Parameters
        ----------
        y_min, y_max : float, optional
            Y-axis limits (unreliability, as a fraction). Must satisfy
            0 < y_min < y_max < 1. Defaults match Analysis's defaults.
        """
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        non_normal = [key for key, val in self.objects.items()
                      if getattr(val, 'dist', 'weibull') != 'normal']
        if non_normal:
            raise ValueError(f'mult_normal() only supports Analysis objects '
                              f'with dist="normal", but {non_normal} do not.')

        def normal_ticks(y_i, _):
            return '{:.1f}'.format(100 * norm.cdf(y_i))

        color, l_style = self._mult_style(color, linestyle)

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=fig_size, num=next(_figure_counter))

        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= y_min)
                                          & (PROBABILITY_PLOT_TICKS <= y_max)]
        z_ticks = norm.ppf(y_ticks)
        plt.yticks(z_ticks, color='black')
        ax.set_yticks([0.0], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')
        plt.tick_params(axis='x', colors='black')

        plt.title(plot_title, color='black', fontsize=plot_title_fontsize)
        plt.xlabel(f'{x_label}{" in "+self.unit if self.unit!="-" else ""}',
                   color='black', fontsize=xy_fontsize)
        plt.ylabel(f'{y_label} in %', color='black', fontsize=xy_fontsize)

        for key, val in self.objects.items():
            col = next(color)
            ls = next(l_style)
            z_line = norm.ppf(np.array([0.0001, 0.9999]))
            x_line = getattr(val, 'mu') + getattr(val, 'sigma') * z_line
            plt.plot(x_line, z_line, color=col, linestyle=ls, linewidth=1.5,
                     label=f'{key}')

            bounds_type = getattr(val, 'bounds_type')
            if bounds_type == '2s' and getattr(val, 'bounds_lower') is not None:
                z_p = norm.ppf(self.unrel)
                plt.plot(getattr(val, 'bounds_lower'), z_p, color=col,
                         linestyle=ls, linewidth=1, label='_nolegend_')
                plt.plot(getattr(val, 'bounds_upper'), z_p, color=col,
                         linestyle=ls, linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=z_p, x1=getattr(val, 'bounds_lower'),
                                   x2=getattr(val, 'bounds_upper'),
                                   alpha=0.1, color=col, label='_nolegend_')
            elif bounds_type == '1su' and getattr(val, 'bounds_upper') is not None:
                z_p = norm.ppf(self.unrel)
                plt.plot(getattr(val, 'bounds_upper'), z_p, color=col,
                         linestyle=ls, linewidth=1, label='_nolegend_')
            elif bounds_type == '1sl' and getattr(val, 'bounds_lower') is not None:
                z_p = norm.ppf(self.unrel)
                plt.plot(getattr(val, 'bounds_lower'), z_p, color=col,
                         linestyle=ls, linewidth=1, label='_nolegend_')

            if plot_ranks:
                if getattr(val, 'ds') is None:
                    ranks = self._mult_median_rank(getattr(val, 'df'))
                else:
                    ranks = self._mult_median_rank_cens(getattr(val, 'df'),
                                                         getattr(val, 'ds'))
                plt.plot(getattr(val, 'df'), norm.ppf(ranks), marker='o',
                         markerfacecolor=col, markeredgecolor='black',
                         markersize=4, alpha=.5, linestyle='None', zorder=3)

        if x_bounds is not None:
            if type(x_bounds) != list:
                raise TypeError(f'x_bounds need to be of type list and not {type(x_bounds)}.')
            plt.xlim(x_bounds[0], x_bounds[1])

        ax.set_ylim(bottom=norm.ppf(y_min), top=norm.ppf(y_max))

        plt.tight_layout()
        plt.grid(True, which='both')
        plt.legend(fontsize=legend_fontsize)

        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, show, save)

    def mult_lognormal(self, x_label='Time To Failure', y_label='Unreliability',
                        plot_title='LogNormal Probability Plot', xy_fontsize=12,
                        plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7),
                        x_bounds=None, plot_ranks=True, save=False, color=None,
                        linestyle=None, y_min=0.01, y_max=0.99, show=True,
                        **kwargs):
        """
        Plots multiple Analysis class objects with dist='lognormal' in one
        figure, analogous to mult_weibull(). Like mult_weibull(), the
        x-axis is log-scaled - see Analysis._plot_lognormal()'s docstring
        for why LogNormal needs this combination of a probit y-axis and a
        log x-axis.

        Parameters
        ----------
        y_min, y_max : float, optional
            Y-axis limits (unreliability, as a fraction). Must satisfy
            0 < y_min < y_max < 1. Defaults match Analysis's defaults.
        """
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        non_lognormal = [key for key, val in self.objects.items()
                         if getattr(val, 'dist', 'weibull') != 'lognormal']
        if non_lognormal:
            raise ValueError(f'mult_lognormal() only supports Analysis objects '
                              f'with dist="lognormal", but {non_lognormal} do not.')

        def normal_ticks(y_i, _):
            return '{:.1f}'.format(100 * norm.cdf(y_i))

        color, l_style = self._mult_style(color, linestyle)

        # Get x-axis bounds to plot
        temp_list = []
        for key, val in self.objects.items():
            bounds_type = getattr(val, 'bounds_type')
            mu, sigma = getattr(val, 'mu'), getattr(val, 'sigma')
            if bounds_type == '2s' and getattr(val, 'bounds_lower') is not None:
                temp_list.append(min(getattr(val, 'bounds_lower')))
                temp_list.append(max(getattr(val, 'bounds_upper')))
            elif bounds_type == '1su' and getattr(val, 'bounds_upper') is not None:
                temp_list.append(np.exp(mu + sigma * norm.ppf(0.001)))
                temp_list.append(max(getattr(val, 'bounds_upper')))
            elif bounds_type == '1sl' and getattr(val, 'bounds_lower') is not None:
                temp_list.append(min(getattr(val, 'bounds_lower')))
                temp_list.append(np.exp(mu + sigma * norm.ppf(0.999)))
            else:
                temp_list.append(np.exp(mu + sigma * norm.ppf(0.001)))
                temp_list.append(np.exp(mu + sigma * norm.ppf(0.999)))

        if x_bounds is None:
            x_axis_min = min(temp_list)
            x_axis_max = max(temp_list)
        else:
            if type(x_bounds) != list:
                raise TypeError(f'x_bounds need to be of type list and not {type(x_bounds)}.')
            x_axis_min, x_axis_max = x_bounds

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=fig_size, num=next(_figure_counter))

        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= y_min)
                                          & (PROBABILITY_PLOT_TICKS <= y_max)]
        z_ticks = norm.ppf(y_ticks)
        plt.yticks(z_ticks, color='black')
        ax.set_yticks([0.0], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        left = 10 ** (np.ceil(np.log10(x_axis_min)) - 1)
        right = 10 ** (np.ceil(np.log10(x_axis_max)))
        plt.xlim(left, right)
        plt.tick_params(axis='x', colors='black')

        plt.title(plot_title, color='black', fontsize=plot_title_fontsize)
        plt.xlabel(f'{x_label}{" in "+self.unit if self.unit!="-" else ""}',
                   color='black', fontsize=xy_fontsize)
        plt.ylabel(f'{y_label} in %', color='black', fontsize=xy_fontsize)

        for key, val in self.objects.items():
            col = next(color)
            ls = next(l_style)
            z_line = norm.ppf(np.array([0.0001, 0.9999]))
            x_line = np.exp(getattr(val, 'mu') + getattr(val, 'sigma') * z_line)
            plt.semilogx(x_line, z_line, color=col, linestyle=ls, linewidth=1.5,
                         label=f'{key}')

            bounds_type = getattr(val, 'bounds_type')
            if bounds_type == '2s' and getattr(val, 'bounds_lower') is not None:
                z_p = norm.ppf(self.unrel)
                plt.semilogx(getattr(val, 'bounds_lower'), z_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
                plt.semilogx(getattr(val, 'bounds_upper'), z_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=z_p, x1=getattr(val, 'bounds_lower'),
                                   x2=getattr(val, 'bounds_upper'),
                                   alpha=0.1, color=col, label='_nolegend_')
            elif bounds_type == '1su' and getattr(val, 'bounds_upper') is not None:
                z_p = norm.ppf(self.unrel)
                plt.semilogx(getattr(val, 'bounds_upper'), z_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
            elif bounds_type == '1sl' and getattr(val, 'bounds_lower') is not None:
                z_p = norm.ppf(self.unrel)
                plt.semilogx(getattr(val, 'bounds_lower'), z_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')

            if plot_ranks:
                if getattr(val, 'ds') is None:
                    ranks = self._mult_median_rank(getattr(val, 'df'))
                else:
                    ranks = self._mult_median_rank_cens(getattr(val, 'df'),
                                                         getattr(val, 'ds'))
                plt.semilogx(getattr(val, 'df'), norm.ppf(ranks), marker='o',
                             markerfacecolor=col, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)

        ax.set_ylim(bottom=norm.ppf(y_min), top=norm.ppf(y_max))

        plt.tight_layout()
        plt.grid(True, which='both')
        plt.legend(fontsize=legend_fontsize)

        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, show, save)

    def mult_exponential(self, x_label='Time To Failure', y_label='Unreliability',
                          plot_title='Exponential Probability Plot', xy_fontsize=12,
                          plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7),
                          x_bounds=None, plot_ranks=True, save=False, color=None,
                          linestyle=None, y_min=0.01, y_max=0.99, show=True,
                          **kwargs):
        """
        Plots multiple Analysis class objects with dist='exponential' in
        one figure, analogous to mult_weibull(). Uses the same Weibull
        probability paper as mult_weibull() (log x-axis, double-log
        y-axis) - see Analysis._plot_exponential()'s docstring for why
        Exponential is plotted as Weibull's beta=1 special case.

        Parameters
        ----------
        y_min, y_max : float, optional
            Y-axis limits (unreliability, as a fraction). Must satisfy
            0 < y_min < y_max < 1. Defaults match Analysis's defaults.
        """
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        non_exponential = [key for key, val in self.objects.items()
                           if getattr(val, 'dist', 'weibull') != 'exponential']
        if non_exponential:
            raise ValueError(f'mult_exponential() only supports Analysis objects '
                              f'with dist="exponential", but {non_exponential} do not.')

        def weibull_prob_paper(x):
            x = np.asarray(x, dtype=float)
            x = np.where(x > .9999, np.nan, x)
            return np.log(-np.log(1 - x))

        def weibull_ticks(y_i, _):
            return '{:.1f}'.format(100 * (1 - np.exp(-np.exp(y_i))))

        color, l_style = self._mult_style(color, linestyle)

        # Get x-axis bounds to plot
        temp_list = []
        for key, val in self.objects.items():
            bounds_type = getattr(val, 'bounds_type')
            theta = getattr(val, 'theta')
            if bounds_type == '2s' and getattr(val, 'bounds_lower') is not None:
                temp_list.append(min(getattr(val, 'bounds_lower')))
                temp_list.append(max(getattr(val, 'bounds_upper')))
            elif bounds_type == '1su' and getattr(val, 'bounds_upper') is not None:
                temp_list.append(theta * (-np.log(1 - 0.001)))
                temp_list.append(max(getattr(val, 'bounds_upper')))
            elif bounds_type == '1sl' and getattr(val, 'bounds_lower') is not None:
                temp_list.append(min(getattr(val, 'bounds_lower')))
                temp_list.append(theta * (-np.log(1 - 0.999)))
            else:
                temp_list.append(theta * (-np.log(1 - 0.001)))
                temp_list.append(theta * (-np.log(1 - 0.999)))

        if x_bounds is None:
            x_axis_min = min(temp_list)
            x_axis_max = max(temp_list)
        else:
            if type(x_bounds) != list:
                raise TypeError(f'x_bounds need to be of type list and not {type(x_bounds)}.')
            x_axis_min, x_axis_max = x_bounds

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=fig_size, num=next(_figure_counter))

        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= y_min)
                                          & (PROBABILITY_PLOT_TICKS <= y_max)]
        lny_ticks = weibull_prob_paper(y_ticks)
        plt.yticks(lny_ticks, color='black')
        ax.set_yticks([weibull_prob_paper(0.632)], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        left = 10 ** (np.ceil(np.log10(x_axis_min)) - 1)
        right = 10 ** (np.ceil(np.log10(x_axis_max)))
        plt.xlim(left, right)
        plt.tick_params(axis='x', colors='black')

        plt.title(plot_title, color='black', fontsize=plot_title_fontsize)
        plt.xlabel(f'{x_label}{" in "+self.unit if self.unit!="-" else ""}',
                   color='black', fontsize=xy_fontsize)
        plt.ylabel(f'{y_label} in %', color='black', fontsize=xy_fontsize)

        for key, val in self.objects.items():
            col = next(color)
            ls = next(l_style)
            p_line = np.array([0.0001, 0.9999])
            x_line = getattr(val, 'theta') * (-np.log(1 - p_line))
            plt.semilogx(x_line, weibull_prob_paper(p_line), color=col,
                         linestyle=ls, linewidth=1.5, label=f'{key}')

            bounds_type = getattr(val, 'bounds_type')
            if bounds_type == '2s' and getattr(val, 'bounds_lower') is not None:
                y_p = weibull_prob_paper(self.unrel)
                plt.semilogx(getattr(val, 'bounds_lower'), y_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
                plt.semilogx(getattr(val, 'bounds_upper'), y_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=y_p, x1=getattr(val, 'bounds_lower'),
                                   x2=getattr(val, 'bounds_upper'),
                                   alpha=0.1, color=col, label='_nolegend_')
            elif bounds_type == '1su' and getattr(val, 'bounds_upper') is not None:
                y_p = weibull_prob_paper(self.unrel)
                plt.semilogx(getattr(val, 'bounds_upper'), y_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')
            elif bounds_type == '1sl' and getattr(val, 'bounds_lower') is not None:
                y_p = weibull_prob_paper(self.unrel)
                plt.semilogx(getattr(val, 'bounds_lower'), y_p, color=col,
                             linestyle=ls, linewidth=1, label='_nolegend_')

            if plot_ranks:
                if getattr(val, 'ds') is None:
                    ranks = self._mult_median_rank(getattr(val, 'df'))
                else:
                    ranks = self._mult_median_rank_cens(getattr(val, 'df'),
                                                         getattr(val, 'ds'))
                plt.semilogx(getattr(val, 'df'), weibull_prob_paper(np.array(ranks)),
                             marker='o', markerfacecolor=col, markeredgecolor='black',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)

        ax.set_ylim(bottom=weibull_prob_paper(y_min), top=weibull_prob_paper(y_max))

        plt.tight_layout()
        plt.grid(True, which='both')
        plt.legend(fontsize=legend_fontsize)

        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, show, save)

    def compare(self, df, ds=None, bounds=None, bounds_type='2s', cl=0.9,
                x_label='Time to Failure', y_label='Unreliability',
                fig_size=(7.7, 7), y_min=0.01, y_max=0.99, plot_ranks=False,
                criteria='aic', plot_pdf=True, pdf_xy_fontsize=12,
                pdf_tick_fontsize=10, pdf_legend_fontsize=9,
                pdf_plot_title_fontsize=14, show=True, save=False,
                plot_style='predictr', **kwargs):
        """
        Fits every distribution predictr supports (see
        Analysis.SUPPORTED_DISTRIBUTIONS) to the given data and plots one
        probability-plot panel per distribution, each on its own native
        paper (unlike mult_weibull(), which overlays several fits of the
        *same* dist on one shared paper), arranged in a grid and ranked by
        AIC (best fit first). Unlike mult_weibull()/contour_plot(), this is
        for comparing *different* distributions fit to one dataset, not the
        same distribution fit to different datasets.

        Parameters
        ----------
        df : list of floats
            Failure data every distribution is fit to.
        ds : list of floats, optional
            Suspensions (right-censored data), if any. The default is None.
        bounds : {None, 'fb'}, optional
            Confidence bound method used for every fit. Only bounds='fb' is
            implemented for dist='normal'/'lognormal'/'exponential' -
            asymptotic Fisher-information bounds for all three there,
            including 'exponential' (see _fisher_bounds_exponential()'s
            docstring) - so that's the only bound type compare() can show
            consistently across all panels. bounds='chi2' is *not* accepted
            here even though Analysis(dist='exponential') supports it too
            (an exact chi-square alternative - see
            _exact_bounds_exponential()'s docstring): it only applies to
            that one distribution, and compare() needs a single bounds=
            value that works across every panel. The default is None (no
            bounds).
        bounds_type : {'2s', '1sl', '1su'}, optional
            Passed through to every fit. The default is '2s'.
        cl : float, optional
            Confidence level for the bounds. The default is 0.9.
        x_label : string, optional
            Label for each panel's x-axis. The default is 'Time to Failure'.
        y_label : string, optional
            Label for each panel's y-axis. The default is 'Unreliability'.
        fig_size : tuple of floats, optional
            Sets width and height in inches: (width, height). The default
            is (7.7, 7), tuned for the common 4-distribution case.
        y_min : float, optional
            Lower y-axis limit (unreliability, as a fraction) shown on
            every panel. Must satisfy 0 < y_min < y_max < 1. The default
            is 0.01.
        y_max : float, optional
            Upper y-axis limit shown on every panel. The default is 0.99.
        plot_ranks : boolean, optional
            If True, each panel also plots the data's median-rank points.
            Median ranks are a plotting-position convenience (where to draw
            each failure on the probability paper) and play no part in the
            MLE fit or the ranking compare() is actually about, so - unlike
            Analysis(plot_ranks=), which defaults to True - this defaults to
            False here. The default is False.
        criteria : {'aic', 'ad'}, optional
            Which goodness-of-fit measure ranks and labels the panels.
            'aic' (the default) ranks by AIC, comparable across the
            different distributions by construction. 'ad' ranks by an
            "adjusted" Anderson-Darling statistic instead (see
            _anderson_darling()'s docstring) - lower is a closer fit - and
            the subtitle notes that, unlike AIC, comparing AD values
            *across* distributions isn't fully rigorous.
        plot_pdf : boolean, optional
            If True, also produces a second, separate figure overlaying
            every fitted distribution's probability density function (PDF)
            in one linear-scale plot - unlike the probability-plot grid
            above, which puts each distribution on its own native paper and
            can't show them on shared axes. Every curve uses the same
            predictr palette color/order as the ranked panels; the raw data
            is marked along the x-axis ('x' for failures, 'o' for
            suspensions, both black). This second figure follows show/save
            the same way the main one does, saved to a '_pdf' suffixed
            filename when save=True. The default is True.
        pdf_xy_fontsize, pdf_tick_fontsize, pdf_legend_fontsize,
        pdf_plot_title_fontsize : float, optional
            Font sizes for the PDF figure (only used when plot_pdf=True).
            Unlike the probability-plot grid's panel titles/labels/ticks/
            legend - deliberately shrunk relative to fig_size via a `scale`
            factor since several panels share that one figure (see
            _set_log_minor_ticks()'s docstring) - the PDF figure is its own
            single, full-size figure, so these default to the same fixed
            sizes any other standalone predictr plot uses (e.g.
            Analysis.plot()'s xy_fontsize/tick_fontsize/legend_fontsize/
            plot_title_fontsize defaults), not a fraction of them.
        show : boolean, optional
            If True, the figure is displayed via plt.show(). The default
            is True.
        save : boolean, optional
            If True, the plot is saved according to the path. The default
            is False.
        plot_style : string, optional
            Passed through to every fit. The default is 'predictr'.
        **kwargs :
            path: string
                Path defines the directory and format of the figure E.g.
                r'var/user/.../test.pdf'
        """
        if not (0 < y_min < y_max < 1):
            raise ValueError('y_min and y_max must satisfy 0 < y_min < '
                              'y_max < 1.')
        if bounds == 'chi2':
            raise ValueError("compare() fits every distribution to the same "
                              "data, so it needs one bounds= value that "
                              "works for all of them - bounds='chi2' (the "
                              "exact chi-square bounds; see "
                              "_exact_bounds_exponential()'s docstring) only "
                              "applies to dist='exponential'. Use "
                              "bounds='fb' instead: it's valid for every "
                              "distribution, including 'exponential' (there "
                              "via the asymptotic Fisher-information bounds "
                              "in _fisher_bounds_exponential(), not the "
                              "exact chi-square ones).")
        if bounds not in (None, 'fb'):
            raise ValueError("compare() only supports bounds=None or "
                              "bounds='fb': that's the only bound type "
                              "implemented for dist='normal'/'lognormal'/"
                              "'exponential', so it's the only one every "
                              "panel can show consistently.")
        if criteria not in ('aic', 'ad'):
            raise ValueError(f'"{criteria}" is not a supported criteria. '
                              f'Supported: \'aic\', \'ad\'.')

        # Weibull/LogNormal/Exponential all raise on negative df/ds (see
        # Analysis.__init__) - only Normal is defined on the whole real
        # line - so compare() fits Normal alone in that case instead of
        # letting the other three raise, and notes why in the subtitle.
        has_negative = any(x < 0 for x in list(df) + list(ds or []))
        dists_to_fit = ['normal'] if has_negative else sorted(Analysis.SUPPORTED_DISTRIBUTIONS)

        objects = {}
        for dist in dists_to_fit:
            fit = Analysis(df=df, ds=ds, dist=dist, bounds=bounds,
                           bounds_type=bounds_type, cl=cl, show=False,
                           plot_style=plot_style)
            fit.mle()
            objects[dist] = fit

        def weibull_paper(x):
            x = np.asarray(x, dtype=float)
            x = np.where(x > .9999, np.nan, x)
            return np.log(-np.log(1 - x))

        def weibull_ticks(y_i, _):
            return '{:.0f}'.format(100 * (1 - np.exp(-np.exp(y_i))))

        def normal_ticks(y_i, _):
            return '{:.0f}'.format(100 * norm.cdf(y_i))

        y_ticks_frac = PROBABILITY_PLOT_TICKS[(PROBABILITY_PLOT_TICKS >= y_min)
                                               & (PROBABILITY_PLOT_TICKS <= y_max)]

        if criteria == 'ad':
            criteria_values = {key: _anderson_darling(obj) for key, obj in objects.items()}
        else:
            criteria_values = {key: obj.aic for key, obj in objects.items()}

        ranked = sorted(objects.items(), key=lambda kv: criteria_values[kv[0]])
        best_value = criteria_values[ranked[0][0]]
        n = len(ranked)
        ncols = 2 if n > 1 else 1
        nrows = int(np.ceil(n / ncols))

        # Every fontsize below is scaled relative to the (11, 10) size this
        # layout was originally tuned at, so the title/subtitle/panel-title/
        # label/legend/tick text - and the fixed-fraction y-positions near
        # the end of this method - all shrink or grow together with
        # fig_size instead of the text staying a fixed point size while the
        # panels around it shrink, which is what caused labels to overlap
        # when fig_size was first made configurable.
        scale = fig_size[1] / 10.0

        _apply_plot_style(plot_style)
        fig, axes = plt.subplots(nrows, ncols, figsize=fig_size, num=next(_figure_counter))
        axes = np.atleast_1d(axes).flatten()
        for extra_ax in axes[n:]:
            extra_ax.set_visible(False)

        for panel_i, (key, obj) in enumerate(ranked):
            ax = axes[panel_i]
            rank = panel_i + 1
            delta = criteria_values[key] - best_value
            dist = obj.dist

            if dist in ('weibull', 'exponential'):
                ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
                y_ticks = weibull_paper(y_ticks_frac)
                ax.set_yticks(y_ticks)
                ax.set_yticks([weibull_paper(0.632)], minor=True)
                ax.set_xscale('log')
                p_line = np.array([0.001, 0.999])
                if dist == 'weibull':
                    x_line = obj.eta * (-np.log(1 - p_line)) ** (1 / obj.beta)
                else:
                    x_line = obj.theta * (-np.log(1 - p_line))
                y_line = weibull_paper(p_line)
            else:
                ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
                y_ticks = norm.ppf(y_ticks_frac)
                ax.set_yticks(y_ticks)
                ax.set_yticks([0.0], minor=True)
                z_line = norm.ppf(np.array([0.001, 0.999]))
                if dist == 'lognormal':
                    ax.set_xscale('log')
                    x_line = np.exp(obj.mu + obj.sigma * z_line)
                else:
                    x_line = obj.mu + obj.sigma * z_line
                y_line = z_line

            ax.grid(True, which='minor', axis='y', linestyle='--')

            # bounds='fb' confidence bounds: asymptotic Fisher-information
            # bounds for every distribution shown here, including
            # Exponential (see _fisher_bounds_exponential()'s docstring -
            # compare() never uses Exponential's exact chi-square
            # alternative, bounds='chi2', since it can't be applied to the
            # other three panels too). Same fit colour + light fill
            # convention as the single-distribution plots. All
            # labeled '_nolegend_' so they don't get auto-picked-up as
            # legend handles below - only the (unlabeled) MLE line becomes
            # the legend's icon.
            if obj.bounds_lower is not None or obj.bounds_upper is not None:
                z_p = (weibull_paper(obj.unrel) if dist in ('weibull', 'exponential')
                       else norm.ppf(obj.unrel))
                if obj.bounds_lower is not None:
                    ax.plot(obj.bounds_lower, z_p, color=PREDICTR_FIT_COLOR,
                            linewidth=1, label='_nolegend_')
                if obj.bounds_upper is not None:
                    ax.plot(obj.bounds_upper, z_p, color=PREDICTR_FIT_COLOR,
                            linewidth=1, label='_nolegend_')
                if obj.bounds_lower is not None and obj.bounds_upper is not None:
                    ax.fill_betweenx(z_p, obj.bounds_lower, obj.bounds_upper,
                                      alpha=0.1, color=PREDICTR_FIT_COLOR, label='_nolegend_')

            if plot_ranks:
                ranks = np.array(obj.median_rank() if obj.ds is None
                                  else obj.median_rank_cens())
                x_data = np.array(obj.df, dtype=float)
                y_data = (weibull_paper(ranks) if dist in ('weibull', 'exponential')
                          else norm.ppf(ranks))
                ax.plot(x_data, y_data, marker='o', markerfacecolor=PREDICTR_FIT_COLOR,
                        markeredgecolor='black', markersize=5, alpha=.6,
                        linestyle='None', zorder=3, label='_nolegend_')
            ax.plot(x_line, y_line, color=PREDICTR_FIT_COLOR, linewidth=1.5, zorder=2)

            if ax.get_xscale() == 'log':
                _set_log_minor_ticks(ax)

            ax.set_ylim(y_ticks[0], y_ticks[-1])
            # Deliberate deviation from the single-plot title convention
            # (plain, non-bold): a comparison grid needs a stronger title
            # hierarchy (suptitle > panel title > axis labels) than one
            # standalone plot does.
            ax.set_title(f'{rank}.  {dist.capitalize()}', color='black',
                         fontsize=15 * scale, fontweight='bold')
            # One rung down from xy_fontsize's default (12) on predictr's
            # existing size ladder (14 title / 12 labels / 10 ticks /
            # 9 legend) - lands on 10.
            ax.set_xlabel(x_label, color='black', fontsize=10 * scale)
            ax.set_ylabel(y_label + ' in %', color='black', fontsize=10 * scale)

            if criteria == 'ad':
                stat_text = (f'log-likelihood = {obj.loglik:.2f}\n'
                             f'AD = {criteria_values[key]:.3f}'
                             + ('  (best)' if rank == 1 else f'   ΔAD = +{delta:.3f}'))
            else:
                stat_text = (f'log-likelihood = {obj.loglik:.2f}\n'
                             f'AIC = {criteria_values[key]:.2f}'
                             + ('  (best)' if rank == 1 else f'   ΔAIC = +{delta:.2f}'))
            ax.legend((stat_text,), loc='lower right', fontsize=9 * scale,
                      title='Goodness of fit', title_fontsize=9 * scale)

            ax.grid(True, which='major')
            ax.tick_params(which='both', labelsize=10 * scale)
            # Large/small axis values (e.g. df in the hundred-thousands, as
            # with dist='normal' on data that isn't log-scaled) make
            # matplotlib draw a "1eN" offset-text next to the tick labels.
            # That Text object isn't covered by tick_params(labelsize=...)
            # above, so without this it stays at its unscaled default size
            # and looks oversized/misplaced once the panel has been shrunk
            # by `scale`.
            ax.xaxis.get_offset_text().set_fontsize(10 * scale)
            ax.yaxis.get_offset_text().set_fontsize(10 * scale)

        susp_num = len(ds) if ds is not None else 0
        fig.suptitle('Distribution comparison', color='black', fontsize=18 * scale,
                     fontweight='bold', y=0.953)
        if has_negative:
            subtitle = ('Only Normal shown - Weibull/LogNormal/Exponential '
                         'require non-negative data  (n = {} (f: {} | s: {}))').format(
                             len(df) + susp_num, len(df), susp_num)
        elif criteria == 'ad':
            subtitle = ('ranked by Anderson-Darling (adjusted) - comparisons '
                         'across distributions are limited  '
                         '(n = {} (f: {} | s: {}))').format(
                             len(df) + susp_num, len(df), susp_num)
        else:
            subtitle = 'ranked by AIC  (n = {} (f: {} | s: {}))'.format(
                len(df) + susp_num, len(df), susp_num)
        fig.text(0.5, 0.908, subtitle,
                 ha='center', color='black', fontsize=12 * scale, fontweight='normal')
        # tight_layout's own computed top margin is tighter than any rect
        # ceiling passed to it once enough vertical room is available, so
        # rect no longer controls the subtitle-to-panel-title gap past
        # that point - only the two fig.text/suptitle y positions above do
        # (verified empirically via the Figure's rendered text bounding
        # boxes, not just visually). The y=0.953/0.908/rect=0.93 fractions
        # are tied to fig_size via `scale` above (all fontsizes scale with
        # it too), so they hold regardless of fig_size instead of only
        # being valid at the one size they were originally tuned at.
        plt.tight_layout(rect=[0, 0, 1, 0.93])

        if save:
            try:
                fig.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')
        result = _finish_plot(fig, show, save)

        if plot_pdf:
            # Second, separate figure: every fitted PDF overlaid on one
            # shared linear-scale axis (the probability-plot grid above
            # can't do this - each distribution needs its own native paper
            # there). Colors follow objects' own insertion order (fixed,
            # alphabetical - see dists_to_fit above), not `ranked`'s
            # criteria-dependent order, so a given distribution keeps the
            # same color across different compare() calls/criteria.
            _apply_plot_style(plot_style)
            fig_pdf = plt.figure(figsize=fig_size, num=next(_figure_counter))
            ax_pdf = fig_pdf.gca()

            def _quantile(o, p):
                if o.dist == 'weibull':
                    return o.eta * (-np.log(1 - p)) ** (1 / o.beta)
                elif o.dist == 'exponential':
                    return o.theta * (-np.log(1 - p))
                elif o.dist == 'normal':
                    return o.mu + o.sigma * norm.ppf(p)
                else:
                    return np.exp(o.mu + o.sigma * norm.ppf(p))

            def _pdf(o, x):
                if o.dist == 'weibull':
                    return weibull_min.pdf(x, o.beta, scale=o.eta)
                elif o.dist == 'exponential':
                    return expon.pdf(x, scale=o.theta)
                elif o.dist == 'normal':
                    return norm.pdf(x, o.mu, o.sigma)
                else:
                    y = np.zeros_like(x)
                    positive = x > 0
                    y[positive] = norm.pdf(np.log(x[positive]), o.mu, o.sigma) / x[positive]
                    return y

            # x-range covers every fitted PDF's own 1%-99% quantile range
            # (tighter than the 0.1%-99.9% tmin_plot/tmax_plot convention
            # used elsewhere in this file, deliberately: on a linear axis, a
            # heavy-tailed candidate's far-out 99.9th percentile compresses
            # every curve's actual peak - the visually interesting part for
            # a shape-vs-shape comparison - into a sliver of the plot) plus
            # the raw data itself, so no rug-plot marker on an outlier gets
            # clipped even if it sits beyond every fitted curve's own 99%.
            x_bounds = ([_quantile(o, p) for o in objects.values() for p in (0.01, 0.99)]
                        + list(df) + list(ds or []))
            x_pdf = np.linspace(min(x_bounds), max(x_bounds), 500)

            colors, _ = _categorical_style(len(objects))
            for (dist_key, obj), color in zip(objects.items(), colors):
                ax_pdf.plot(x_pdf, _pdf(obj, x_pdf), color=color, linewidth=1.5,
                            label=obj.dist.capitalize())

            # Raw data as a rug plot along the x-axis - black regardless of
            # which curve's color it's nearest, since it's the same
            # observed data for every distribution, not tied to one of them.
            ax_pdf.plot(df, [0] * len(df), marker='x', color='black',
                        linestyle='None', markersize=7, label='Failures')
            if ds:
                ax_pdf.plot(ds, [0] * len(ds), marker='o', markerfacecolor='none',
                            markeredgecolor='black', linestyle='None',
                            markersize=7, label='Suspensions')

            ax_pdf.set_xlabel(x_label, color='black', fontsize=pdf_xy_fontsize)
            ax_pdf.set_ylabel('Probability Density Function', color='black',
                               fontsize=pdf_xy_fontsize)
            ax_pdf.set_title('Probability Density Function Comparison', color='black',
                              fontsize=pdf_plot_title_fontsize, fontweight='bold')
            ax_pdf.tick_params(labelsize=pdf_tick_fontsize)
            ax_pdf.grid(True)
            ax_pdf.legend(fontsize=pdf_legend_fontsize)
            fig_pdf.tight_layout()

            if save:
                try:
                    path = kwargs['path']
                    if '.' in path:
                        base, ext = path.rsplit('.', 1)
                        pdf_path = f'{base}_pdf.{ext}'
                    else:
                        pdf_path = f'{path}_pdf'
                    fig_pdf.savefig(pdf_path)
                except:
                    raise ValueError('Path is faulty.')
            _finish_plot(fig_pdf, show, save)

        return result

    def contour_plot(self, show=True, style='hull', show_weibull=False, show_legend=True, color=None, x_label=None,
                     y_label=None, plot_title='Contour Plot', xy_fontsize=12,
                     plot_title_fontsize=14, legend_fontsize=9, fig_size=(6.4, 4.8), save=False,
                     scale_mode='auto', log_ratio_threshold=10, cl_set=None,
                     curve_fill=True, fill_alpha=0.25, **kwargs):
        """
        Plots the likelihood-ratio-bounds contour (bounds='lrb') for one or
        more already-fitted Analysis objects, in whichever of the two
        parameter spaces its dist uses: Weibull's (beta, eta), or Normal/
        LogNormal's (mu, sigma) - see Analysis._lrb_param_names().

        Every object plotted together must have the exact same dist. This
        is a hard requirement, not just a default: Weibull's beta (a
        unitless shape parameter) and Normal/LogNormal's mu (a location
        parameter in the data's own units, or in ln(data)'s units for
        LogNormal) aren't the same kind of quantity, and even Normal vs.
        LogNormal - despite both using a (mu, sigma) pair - live on
        different scales (LogNormal's mu/sigma are of ln(data), Normal's
        are of the data itself; e.g. a LogNormal mu of ~4.6 and a Normal mu
        of ~100 can be the *same* underlying data). Overlaying curves from
        different dists on one shared pair of axes would visually compare
        numbers that don't mean the same thing, so it's rejected outright
        (ValueError) rather than silently drawn - there's no established
        reliability-engineering practice (see e.g. Meeker & Escobar) for
        putting them on one plot; every reference treats each fitted
        model's contour separately. Cross-distribution comparison belongs
        in likelihood/AIC space (see PlotAll.compare()), not parameter
        space.

        Every plotted object must also already have its LRB contour
        computed (i.e. .mle() was called with bounds='lrb', or .lrb() was
        called directly) - contour_plot() checks this explicitly per
        object (see Analysis._lrb_arrays()) and raises a clear error
        rather than crashing on a missing value.

        Parameters
        ----------
        x_label, y_label : str, optional
            Axis labels. Default (None) picks the label pair matching the
            objects' shared dist (see Analysis._lrb_axis_labels()); only
            override if you want something other than that default.
        scale_mode : {'auto', 'linear', 'log'}, optional
            Scaling of the spread-parameter (y) axis - eta for Weibull,
            sigma for Normal/LogNormal. 'auto' (default) inspects that
            range across all objects and switches to a logarithmic scale
            when it spans more than `log_ratio_threshold`x in magnitude.
            This avoids small-scale ellipses being crushed into a flat
            line next to large-scale ones when several objects with very
            different magnitudes are plotted together. 'linear' and 'log'
            force the respective scale regardless of the data.
        log_ratio_threshold : float, optional
            Ratio of max/min of the spread parameter across all objects
            above which 'auto' switches to a logarithmic scale.
            Default = 10.
        cl_set : list of float, optional
            Confidence levels to draw per dataset, e.g. [0.95, 0.9, 0.8].
            If None or empty (default), each object's own `cl` attribute is
            used, i.e. exactly one curve per object, unchanged from before.
            When given, `cl_set` is sorted from largest to smallest and, for
            every object, its contour is recomputed and drawn once per
            confidence level - via a temporary deep copy of the object, so
            the original object's own `cl` and LRB arrays are left
            untouched. All curves of a dataset share the same base hue and
            are progressively shaded darker from the largest to the smallest
            confidence level (see the `shade_of` helper below), so that the
            different confidence levels of the same dataset stay clearly
            distinguishable without ever washing out towards black or grey.
        curve_fill : bool, optional
            If True, the area enclosed by each curve is filled with that
            curve's own (possibly shaded) color at `fill_alpha` opacity.
            Default = False.
        fill_alpha : float, optional
            Opacity used for the curve fill when `curve_fill` is True.
            Default = 0.25.

        """
        if scale_mode not in ('auto', 'linear', 'log'):
            raise ValueError("scale_mode must be one of 'auto', 'linear', 'log'.")

        if cl_set:
            if any(not (0 < c < 1) for c in cl_set):
                raise ValueError('All values in cl_set must be between 0 and 1.')
            cl_set = sorted(cl_set, reverse=True)

        # All objects must share one dist - see docstring for why mixing
        # parameter spaces on one plot isn't meaningful.
        dists = {getattr(val, 'dist') for val in self.objects.values()}
        if len(dists) > 1:
            raise ValueError(f'contour_plot() cannot mix dists {sorted(dists)} '
                              f'on one plot - Weibull\'s (beta, eta), Normal\'s '
                              f'(mu, sigma), and LogNormal\'s (mu, sigma) are not '
                              f'the same kind of quantity or on the same scale '
                              f'(see contour_plot()\'s docstring). Plot each dist '
                              f'separately.')
        if not dists:
            raise ValueError('contour_plot() has no objects to plot.')
        shared_dist = next(iter(dists))
        x_name, y_name = Analysis._lrb_param_names(shared_dist)

        # Configure plot
        _apply_plot_style(self.plot_style)
        fig, ax = plt.subplots(figsize=fig_size, num=next(_figure_counter))
        ax.set_title(plot_title, fontsize=plot_title_fontsize)


        # Set color/linestyle per dataset. When color isn't given explicitly,
        # predictr's categorical palette is used and, for more than 6
        # datasets, the linestyle is cycled alongside the color (see
        # _categorical_style()) so identity stays distinguishable instead of
        # silently repeating a bare color.
        num_objects = len(self.objects)
        if color is not None:
            color = iter(color)
            base_linestyles = iter(['-'] * num_objects)
        else:
            default_colors, default_linestyles = _categorical_style(num_objects)
            color = iter(default_colors)
            base_linestyles = iter(default_linestyles)

        def shade_of(c, t, value_drop=0.5, saturation_gain=0.2):
            """
            Returns a shade of color `c` for t in [0, 1]: t=0 is the color
            itself, t=1 is the darkest shade. Works in HSV space and only
            ramps value down (by up to `value_drop`, i.e. never below
            (1 - value_drop) of the original brightness) while nudging
            saturation up (by up to `saturation_gain`), keeping the hue
            fixed. This keeps a family of shades of the same base color
            (e.g. several cl_set levels of one dataset) clearly
            distinguishable and visibly the same hue - unlike naively
            multiplying all RGB channels by a shrinking factor, which
            desaturates towards black/grey and erases the hue after a few
            steps.
            """
            r, g, b = mpl.colors.to_rgb(c)
            h, s, v = colorsys.rgb_to_hsv(r, g, b)
            v = v * (1 - value_drop * t)
            s = min(1.0, s + saturation_gain * t)
            return colorsys.hsv_to_rgb(h, s, v)

        # Collect the (cl, x_lrb, y_lrb) curve(s) to draw per object - x/y
        # being (beta, eta) for Weibull or (mu, sigma) for Normal/LogNormal
        # (see _lrb_param_names()). If cl_set is given, the LRB region is
        # recomputed for every requested confidence level on a deep copy of
        # the object, leaving the original object untouched; otherwise the
        # object's own cl and already-computed LRB arrays (checked via
        # _lrb_arrays() - raises a clear error if .lrb() was never called)
        # are used, as before.
        curves = {}
        for key, val in self.objects.items():
            if cl_set:
                object_curves = []
                for cl_value in cl_set:
                    val_cl = copy.deepcopy(val)
                    val_cl.cl = cl_value
                    val_cl.lrb()
                    x_lrb, y_lrb = Analysis._lrb_arrays(val_cl)
                    object_curves.append((cl_value, x_lrb, y_lrb))
                curves[key] = object_curves
            else:
                x_lrb, y_lrb = Analysis._lrb_arrays(val)
                curves[key] = [(getattr(val, 'cl'), x_lrb, y_lrb)]

        # Analyze the spread-parameter (eta/sigma) data across all objects/
        # confidence levels to decide on the y-axis scale
        all_eta = np.concatenate([np.asarray(eta) for object_curves in curves.values()
                                  for _, _, eta in object_curves])
        if scale_mode == 'log':
            use_log = True
        elif scale_mode == 'linear':
            use_log = False
        else:
            use_log = (all_eta.max() / all_eta.min()) > log_ratio_threshold
        if use_log and np.any(all_eta <= 0):
            raise ValueError(f'Logarithmic scale requires strictly positive '
                              f'{y_name} values.')

        def label_anchors(object_curves):
            """
            Picks, for one object's family of curves, which of the 4 cardinal
            directions (top/bottom/left/right of each curve) best separates
            the inline confidence-level labels from each other, and returns
            one anchor point + (ha, va) alignment per curve in that
            direction.

            Nested cl_set curves of the same dataset can be strongly
            elongated (e.g. a wide beta range but a narrow eta range, or vice
            versa), so always anchoring at the topmost point can bunch all
            labels together on one side while leaving them far apart on
            another. Comparing the (scale-appropriate) spread of the topmost/
            bottommost/leftmost/rightmost points across the curve family and
            using whichever direction spreads them out the most keeps the
            labels legible regardless of the region's shape.
            """
            directions = {
                'top': (lambda b, e: np.argmax(e), 'y', 'center', 'bottom'),
                'bottom': (lambda b, e: np.argmin(e), 'y', 'center', 'top'),
                'right': (lambda b, e: np.argmax(b), 'x', 'left', 'center'),
                'left': (lambda b, e: np.argmin(b), 'x', 'right', 'center'),
            }
            best_name, best_pts, best_spread = None, None, -1
            for name, (pick_idx, coord, ha, va) in directions.items():
                pts = []
                for _, beta, eta in object_curves:
                    beta_arr, eta_arr = np.asarray(beta), np.asarray(eta)
                    idx = pick_idx(beta_arr, eta_arr)
                    pts.append((beta_arr[idx], eta_arr[idx]))
                if coord == 'y':
                    vals = [np.log10(p[1]) for p in pts] if use_log else [p[1] for p in pts]
                else:
                    vals = [p[0] for p in pts]
                spread = (max(vals) - min(vals)) if len(vals) > 1 else 0
                if spread > best_spread:
                    best_name, best_pts, best_spread = name, pts, spread
            ha, va = directions[best_name][2], directions[best_name][3]
            return best_pts, ha, va

        # Curves only use shading up to CURVE_T_MAX, reserving headroom up to
        # t=1.0 for the point estimate marker (see below), so the marker is
        # always one visible tick darker than the dataset's darkest curve,
        # regardless of how many cl_set levels that dataset has.
        CURVE_T_MAX = 0.8

        for key, val in self.objects.items():
            base_color = next(color)
            base_linestyle = next(base_linestyles)
            n_curves = len(curves[key])
            # t=0 for the (only) curve when cl_set is not used, so the
            # single-curve case keeps using the plain base color as before.
            t_values = np.linspace(0, CURVE_T_MAX, n_curves) if n_curves > 1 else [0.0]
            anchor_pts, anchor_ha, anchor_va = label_anchors(curves[key])

            for idx, (conf_level, beta, eta) in enumerate(curves[key]):
                curve_color = shade_of(base_color, t_values[idx])
                # One legend entry per dataset (its name from the objects
                # dict, as before), attached to its first/lightest curve;
                # the confidence level itself is now labeled inline on the
                # curve (see below), not in the legend anymore.
                label = key if idx == 0 else None

                if curve_fill or style == 'hull':
                    # Create Convex Hull around points
                    points = np.column_stack((beta, eta))
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices]
                    hull_points = np.vstack([hull_points, hull_points[0]])

                    if curve_fill:
                        ax.fill(hull_points[:, 0], hull_points[:, 1],
                                color=curve_color, alpha=fill_alpha, zorder=1)

                if style == 'scatter':
                    ax.scatter(beta, eta, s=3, c=[curve_color], label=label, zorder=3)
                elif style == 'hull':
                    ax.plot(hull_points[:, 0], hull_points[:, 1], c=curve_color,
                             linewidth=1.5, markersize=4, label=label, zorder=3,
                             linestyle=base_linestyle)

                # Inline confidence-level label (analogous to
                # matplotlib.axes.Axes.clabel on a regular contour plot),
                # placed at whichever side (top/bottom/left/right) best
                # separates this dataset's nested curves (see label_anchors).
                # A small pill in the axes background color sits behind the
                # text - the same trick clabel uses to break the line - so
                # the label stays readable without covering more of the
                # curve/grid than necessary.
                anchor_x, anchor_y = anchor_pts[idx]
                ax.text(anchor_x, anchor_y, f'{conf_level*100:.0f}%',
                        color=curve_color, fontsize=xy_fontsize * 0.6,
                        ha=anchor_ha, va=anchor_va, zorder=4,
                        bbox=dict(boxstyle='round,pad=0.15', facecolor=ax.get_facecolor(),
                                  edgecolor='none', alpha=0.85))

            # Point estimate marker (bias-corrected if a bcm was used - only
            # applies to Weibull, see _lrb_point()), always one tick darker
            # (t=1.0) than the dataset's darkest curve (which tops out at
            # CURVE_T_MAX < 1.0)
            point_x, point_y = Analysis._lrb_point(val)
            ax.scatter(point_x, point_y, s=40,
                        c=[shade_of(base_color, 1.0)], marker='o', zorder=5)

        default_x_label, default_y_label = Analysis._lrb_axis_labels(shared_dist)
        if x_label is None:
            x_label = default_x_label
        ax.set_xlabel(x_label, fontsize=xy_fontsize)
        if y_label is None:
            # default_y_label is a full '$...$' mathtext span; splice the
            # log(...) wrapper *inside* those delimiters (not around them)
            # so the result stays one single math span - r'$\log($' here
            # would leave a stray '$' right after the opening paren,
            # splitting this into two separate (and here, one malformed)
            # mathtext spans instead of one.
            y_label = (r'$\log(' + default_y_label[1:-1] + r')$') if use_log else default_y_label
        ax.set_ylabel(y_label, fontsize=xy_fontsize)
        if use_log:
            ax.set_yscale('log')
        ax.grid(True, which='both')
        fig.tight_layout()

        if show_legend:
            ax.legend(fontsize=legend_fontsize)

        # Save plot
        if save:
            try:
                fig.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')
        if show_weibull==True:
            if show or save:
                plt.close(fig)
            return fig
        return _finish_plot(fig, show, save)

    def weibull_pdf(self, beta=None, eta=None, linestyle=['-', '--', ':', '-.'], labels = None,
                    x_label = None, y_label=None, xy_fontsize=12, tick_fontsize=10,
                    legend_fontsize=9,
                    plot_title='Weibull PDF', plot_title_fontsize=14, x_bounds=None,
                    fig_size=None, color=None, save=False, show=True, plot_style='predictr', **kwargs):
        """
        Parameters
        ----------
        beta : list of floats
            Weibull shape parameter.
        eta : list of floats
            Weibull scale parameter
        linestyle : list of strings, optional
            Defines the linestyle(s) in the plot.
        labels : list of strings, optional
            List containing the labels for the plot legend.
        x_label : string, optional
            Label for the x-axis. The default is None.
        y_label : string, optional
            Label for the y-axis. The default is None.
        xy_fontsize : float, optional
            Fontsize for the axes labels. The default is 12.
        tick_fontsize : float, optional
            Fontsize for the tick labels. The default is 10.
        legend_fontsize : float, optional
            Fontsize for the legend. The default is 9.
        plot_title : string, optional
            Title for the plot. The default is 'Weibull PDF'.
        plot_title_fontsize : float, optional
            Fontsize of the plot title. The default is 14.
        x_bounds : list of floats,
            Sets x-axis boundaries: [start, end, steps]
        fig_size : tuple of floats, optional
            Sets width and height in inches: (width, height)
        color : list of strings, optional
            List containing the colormap for the plotted lines. Length of list must be equal to
            the beta and eta length of lists.
        save : boolean, optional
            If True, the plot is saved according to the path. The default is False.
        plot_style : TYPE, optional
            DESCRIPTION. The default is 'predictr'.
        **kwargs :
            path: string
                Path defines the directory and format of the figure E.g. r'var/user/.../test.pdf'
        """

        def wei_pdf(x, beta, eta):
            """
            Weibull probability density function.
            """
            return beta / eta * (x / eta) ** (beta - 1) * np.exp(-1 * ((x / eta) ** beta))

        # Check needed input data
        if beta is None or eta is None:
            raise ValueError('Beta and eta must be specified.')

        if x_bounds is None:
            raise ValueError('X axis bounds are not defined. \
                             Use x_bounds argument for this purpose.')

        # Check if number of linestyles and object count match eachother
        if len(linestyle) != len(beta) or len(linestyle) != len(eta):
            raise ValueError('Number of linestyles must match the list length of beta and eta.')

        # Set line color. linestyle is required to already match len(beta),
        # so distinguishing series beyond the 6-color palette is left to the
        # caller's own linestyle list here; itertools.cycle just keeps the
        # color from running out (instead of raising StopIteration) if more
        # than 6 curves are plotted.
        if color is not None:
            color = iter(color)
        else:
            color = itertools.cycle(PREDICTR_PALETTE)
        # Set x-axis
        xvals = np.linspace(x_bounds[0], x_bounds[1], x_bounds[2])

        # Configure plot
        _apply_plot_style(plot_style)

        # Always start a fresh figure/axes - without this, if fig_size is
        # left at its default (None), the plot would silently land on
        # whatever axes happen to still be active (e.g. a previous
        # mle()/mrr() Weibull probability plot), inheriting its log x-scale
        # and probability-paper y-axis formatting instead of a clean linear
        # PDF plot.
        if fig_size is None:
            fig_size = (6.4, 4.8)
        fig = plt.figure(figsize=fig_size, num=next(_figure_counter))

        # Set title
        plt.title(plot_title, fontsize=plot_title_fontsize)

        # Set x and y axis labels and fontsizes
        if x_label is not None:
            plt.xlabel(x_label, fontsize=xy_fontsize)
        else:
            plt.xticks(fontsize=tick_fontsize)

        if y_label is not None:
            plt.ylabel(y_label, fontsize=xy_fontsize)
        else:
            plt.yticks(fontsize=tick_fontsize)

        # Check if multiple lines need to be plotted
        if type(beta)==list and type(eta) == list:
            if labels is not None:
                for i, j, lab, line in zip(beta, eta, labels, linestyle):
                    plt.plot(xvals, wei_pdf(xvals, i, j),
                             linestyle=line, label=lab, color=next(color), linewidth=1.5)

                # Set legend

                plt.legend(fontsize=legend_fontsize)
            else:
                for i, j, line in zip(beta, eta, linestyle):
                    plt.plot(xvals, wei_pdf(xvals, i, j),
                             linestyle=line, color=next(color), linewidth=1.5)
        plt.tight_layout()

        # Save plot
        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return _finish_plot(fig, show, save)

    def simple_weibull(self, beta, eta, unit='-', x_label = 'Time to Failure',
                       y_label = 'Unreliability', xy_fontsize=12, tick_fontsize=10,
                       plot_title_fontsize=14, plot_title='Weibull Probability Plot',
                       fig_size=(6, 7), show_legend=True, legend_fontsize=9,
                       save=False, df=None, ds=None, **kwargs):
        """
        beta : float
            Weibull shape parameter.
        eta : float
            Weibull scale parameter
        x_label : string, optional
            Label for the x-axis. The default is None.
        y_label : string, optional
            Label for the y-axis. The default is None.
        xy_fontsize : float, optional
            Fontsize for the axes labels. The default is 12.
        tick_fontsize : float, optional
            Fontsize for the tick labels. The default is 10.
        legend_fontsize : float, optional
            Fontsize for the legend. The default is 9.
        plot_title : string, optional
            Title for the plot. The default is 'Weibull PDF'.
        plot_title_fontsize : float, optional
            Fontsize of the plot title. The default is 14.
        size : tuple of floats, optional
            Sets width and height in inches: (width, height)
        save : boolean, optional
            If True, the plot is saved according to the path. The default is False.
        plot_style : TYPE, optional
            DESCRIPTION. The default is 'predictr' (inherited from Analysis).
        unit : TYPE, optional
            DESCRIPTION. The default is '-'.
        show_legend : TYPE, optional
            DESCRIPTION. The default is True.
        df : list of floats, optional
            Contains the failures. If None, there will be no median ranks in the plot
        ds : list of floats, optional
            Contains suspensions. The default is None.
        **kwargs :
            path: raw-string
                Path defines the directory and format of the figure E.g. r'var/user/.../test.pdf'

        """

        # Create dummy object
        if df is None:
            df = []
        if ds is None:
            ds=[]
        x = Analysis(df=df, ds=ds, show_legend=show_legend)

        # Set attributes from input
        setattr(x, 'beta', beta)
        setattr(x, 'eta', eta)
        setattr(x, 'x_label', x_label)
        setattr(x, 'y_label', y_label)
        setattr(x, 'xy_fontsize', xy_fontsize)
        setattr(x, 'tick_fontsize', tick_fontsize)
        setattr(x, 'plot_title', plot_title)
        setattr(x, 'plot_title_fontsize', plot_title_fontsize)
        setattr(x, 'fig_size', fig_size)
        setattr(x, 'legend_fontsize', legend_fontsize)
        setattr(x, 'unit', unit)

        # Plot object
        fig = x.plot()

        # Save plot
        if save:
            try:
                fig.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        return fig

class Regression:
    """
    Lifetime (survival) regression with covariates, for uncensored and
    right-censored data. Two models, selected via model=:

      'weibull_aft' - a parametric Weibull accelerated-failure-time model.
          In log-location-scale form (see Meeker & Escobar, "Statistical
          Methods for Reliability Data", ch. 4; Lawless, "Statistical
          Models and Methods for Lifetime Data", 2nd ed., ch. 6):

              ln T_i = x_i' theta + sigma * W_i

          with W_i standard smallest-extreme-value distributed (density
          f_W(w) = exp(w - e^w), survival S_W(w) = exp(-e^w)). sigma is
          the scale; the equivalent Weibull shape is beta = 1 / sigma.
          exp(theta_j) is the time ratio (acceleration factor) for a
          one-unit increase in covariate j.

      'cox_ph' - a semiparametric proportional-hazards model
              h(t | x) = h0(t) * exp(x' beta)
          with an unspecified baseline hazard h0 (no intercept). beta is
          estimated from the partial likelihood (Cox, 1972); exp(beta_j)
          is the hazard ratio. Tied event times are handled by the Efron
          (default) or Breslow approximation (see ties= and
          Kalbfleisch & Prentice, "The Statistical Analysis of Failure
          Time Data", 2nd ed., ch. 4). The cumulative baseline hazard is
          the Breslow estimator; S(t | x) = S0(t) ** exp(x' beta).

    predictr conventions carry over from Analysis: the failures/suspensions
    data split (df/ds), cl/bounds/bounds_type for confidence bounds, and
    the shared plot kwargs (plot_style, fig_size, fonts, legend, show/save).

    Data can be supplied either the predictr way - failures in df,
    right-censored observations in ds, with covariate rows aligned in
    x_df/x_ds - or as one DataFrame via data=/duration_col/event_col
    /covariate_cols (event_col: 1 = failure, 0 = censored). String/category
    covariate columns are one-hot encoded (first level dropped).

    fit() estimates the model and fills the result attributes; summary()
    prints the full report and returns it as a DataFrame.
    """

    SUPPORTED_MODELS = {'weibull_aft', 'cox_ph'}
    SUPPORTED_TIES = {'efron', 'breslow'}
    SUPPORTED_BOUNDS = {None, 'fb', 'lrb', 'npbb', 'pbb'}
    # |beta| past this in the Cox fit is treated as "ran off to infinity".
    _COX_DIVERGENCE = 50.0

    # Boltzmann constant in eV/K - Arrhenius uses 1/(k_B T) so the fitted
    # coefficient is the activation energy Ea in electron-volts.
    _K_B_EV = 8.617333262e-5
    # Named life-stress (aging) laws. 'suffixes' name the design columns
    # built from one raw stress column ('<raw><suffix>'); 'kelvin' marks
    # the temperature laws that need stress_units handling.
    _STRESS_LAWS = {
        'arrhenius':     dict(suffixes=('_invkT',),          kelvin=True),
        'eyring':        dict(suffixes=('_invkT', '_lnT'),    kelvin=True),
        'inverse_power': dict(suffixes=('_ln',),             kelvin=False),
        'coffin_manson': dict(suffixes=('_ln',),             kelvin=False),
        'exponential':   dict(suffixes=('_lin',),            kelvin=False),
    }

    # Styling of the "last observed time" marker that plot_survival draws
    # when the time grid runs past the data: a dotted vertical at that time,
    # tagged with a small vertical 't_max' label reading bottom-to-top,
    # sitting in the gap near the top of the axes. label_y is a 0..1
    # fraction of the survival axis. Everything right of the line is
    # extrapolation.
    _EXTRAP_STYLE = dict(
        line_color='#8b8b8b', line_style=':', line_width=1.0,
        label=r'$t_\mathrm{max}$', label_color='#5a5a5a', label_y=0.95,
        label_fs=None,
    )

    def __init__(self, df=None, ds=None, x_df=None, x_ds=None,
                 data=None, duration_col=None, event_col=None,
                 covariate_cols=None, feature_names=None,
                 stress_model=None, stress_units=None,
                 model='weibull_aft', ties='efron', fit_intercept=True,
                 standardize=True, cl=0.9, bounds=None, bounds_type='2s',
                 max_iter=100, tol=1e-8, n_boot=1000, strata=None, entry_col=None,
                 show=False, save=False, plot_style='predictr', unit='-',
                 x_label='Time to Failure', y_label='Unreliability',
                 xy_fontsize=12, tick_fontsize=10, legend_fontsize=9,
                 plot_title=None, plot_title_fontsize=14,
                 fig_size=(9, 6), show_legend=True, **kwargs):
        """
        Parameters
        ----------
        df, ds : list of floats, optional
            Failures (df) and right-censored observations (ds), predictr
            style. Use together with x_df/x_ds. Mutually exclusive with
            data=.
        x_df, x_ds : DataFrame | dict of columns | 2D array-like, optional
            Covariate rows aligned row-for-row with df and ds. Same columns
            (order and meaning) in both. x_df is required with df; x_ds is
            required whenever ds is given.
        data : DataFrame, optional
            One row per unit. Mutually exclusive with df/ds. Requires
            duration_col and event_col.
        duration_col, event_col : str, optional
            Column names in data for the observed time and the event
            indicator (1 = failure, 0 = censored). Any non-zero event value
            is treated as a failure.
        covariate_cols : list of str, optional
            Columns in data to use as covariates. Defaults to every column
            except duration_col/event_col (and strata/entry_col once those
            are supported).
        feature_names : list of str, optional
            Names for the covariate columns when they cannot be inferred
            (x_df given as a bare array). Ignored when one-hot encoding
            changes the column count.
        stress_model : dict or str, optional
            {raw_stress_column: law_name} to fit a named life-stress
            (aging) relationship. A bare law name (str) is accepted as a
            shorthand when there is exactly one covariate column, and
            applies to that column. Each raw stress column is replaced by
            its physically motivated transform before fitting, so the
            coefficient carries the physical meaning:
              'arrhenius'     - temperature, term 1/(k_B T); coef = Ea [eV]
              'eyring'        - temperature, terms 1/(k_B T) and ln T
              'inverse_power' - positive stress, term ln(S); exponent n = -coef
              'coffin_manson' - thermal-cycling range, term ln(range)
              'exponential'   - stress used as-is (log-linear life)
            Temperature-humidity (Peck) = 'arrhenius' on T plus
            'inverse_power' on RH. predict_*/plot_survival then accept the
            operating point in raw units. The physical-parameter table
            (stress_params), acceleration_factor, plot_stress_life and
            check_shape are Weibull-AFT only; with model='cox_ph' the
            transformed columns still fit and predict_hazard_ratio accepts
            raw units, but the coefficients are hazard-scale, not
            life-stress-law parameters.
        stress_units : dict, optional
            {raw_stress_column: 'C' | 'K'} for the temperature laws.
            Default 'C' (Celsius, converted to Kelvin internally).
        model : {'weibull_aft', 'cox_ph'}, optional
            The default is 'weibull_aft'.
        ties : {'efron', 'breslow'}, optional
            Tie-handling for model='cox_ph'. The default is 'efron' (closer
            to the exact partial likelihood; 'breslow' is faster and biases
            coefficients toward zero as ties accumulate).
        fit_intercept : bool, optional
            model='weibull_aft' only. The default is True.
        standardize : bool, optional
            Center/scale covariates internally for numerical conditioning;
            coefficients and bounds are reported on the original scale. The
            default is True.
        cl : float, optional
            Confidence level for the bounds. The default is 0.9.
        bounds : {None, 'fb', 'lrb', 'npbb', 'pbb'}, optional
            None - no confidence bounds (the default; matches Analysis).
            'fb' - Wald bounds from the observed information (cheap).
            'lrb' - profile-likelihood bounds (more accurate with few
            events, slower). 'npbb' - non-parametric bootstrap percentile
            bounds (resample units, refit n_boot times). 'pbb' - parametric
            bootstrap (simulate the response from the fitted model).
        bounds_type : {'2s', '1sl', '1su'}, optional
            Two-sided, one-sided lower or one-sided upper. The default is
            '2s'. Ignored by 'npbb'/'pbb' (always two-sided percentiles).
        max_iter, tol : int, float, optional
            Newton-Raphson iteration cap and gradient tolerance.
        n_boot : int, optional
            Number of resamples for bounds='npbb'/'pbb'. The default is
            1000.
        strata, entry_col : optional
            Reserved for stratification and left-truncation. Not supported
            yet - passing a non-None value raises.
        show, save, plot_style, unit, x_label, y_label, xy_fontsize,
        tick_fontsize, legend_fontsize, plot_title, plot_title_fontsize,
        fig_size, show_legend : see Analysis - identical meaning. fig_size
            defaults to landscape (9, 6) here.
        **kwargs :
            path : str
                Output path/format when save=True.
        """
        if model not in self.SUPPORTED_MODELS:
            raise ValueError(f'"{model}" is not a supported model. '
                             f'Supported: {sorted(self.SUPPORTED_MODELS)}')
        if ties not in self.SUPPORTED_TIES:
            raise ValueError(f'"{ties}" is not a supported tie-handling. '
                             f'Supported: {sorted(self.SUPPORTED_TIES)}')
        if bounds not in self.SUPPORTED_BOUNDS:
            raise ValueError(f'"{bounds}" is not supported by Regression. '
                             f'Supported: {sorted(b for b in self.SUPPORTED_BOUNDS if b)} '
                             f'or None.')
        if bounds is not None and bounds_type not in ('2s', '1sl', '1su'):
            raise ValueError("bounds_type must be '2s', '1sl' or '1su'.")
        if strata is not None or entry_col is not None:
            raise ValueError('strata= and entry_col= are reserved for a '
                             'future version and not supported yet.')

        got_split = df is not None or ds is not None
        got_frame = data is not None
        if got_split == got_frame:
            raise ValueError('Provide the data either as df/ds (+ x_df/x_ds) '
                             'or as data= (+ duration_col/event_col), not '
                             'both and not neither.')

        self.model = model
        self.ties = ties
        self.fit_intercept = fit_intercept if model == 'weibull_aft' else False
        self.standardize = standardize
        self.cl = cl
        self.bounds = bounds
        self.bounds_type = bounds_type if bounds is not None else None
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.n_boot = int(n_boot)

        # Plot attributes - same names/defaults as Analysis.
        self.show = show
        self.save = save
        self.plot_style = plot_style
        self.unit = unit
        self.x_label, self.y_label = x_label, y_label
        self.xy_fontsize, self.tick_fontsize = xy_fontsize, tick_fontsize
        self.legend_fontsize = legend_fontsize
        self.plot_title_fontsize = plot_title_fontsize
        self.fig_size, self.show_legend = fig_size, show_legend
        if plot_title is None:
            plot_title = {'weibull_aft': 'Weibull AFT — Coefficients',
                          'cox_ph': 'Cox PH — Coefficients'}[model]
        self.plot_title = plot_title
        if self.save:
            try:
                self.save_path = kwargs['path']
            except KeyError:
                raise ValueError('Path is not defined.')

        # --- ingest: build one internal (t, event, covariate frame) ---
        if got_frame:
            if not isinstance(data, pd.DataFrame):
                raise ValueError('data= must be a pandas DataFrame.')
            if duration_col is None or event_col is None:
                raise ValueError('data= requires duration_col and event_col.')
            frame = data.reset_index(drop=True)
            t = np.asarray(frame[duration_col], dtype=float)
            event = (np.asarray(frame[event_col]) != 0).astype(int)
            if covariate_cols is None:
                covariate_cols = [c for c in frame.columns
                                  if c not in (duration_col, event_col)]
            x_all = frame[list(covariate_cols)].reset_index(drop=True)
        else:
            df = list(df) if df is not None else []
            ds = list(ds) if ds is not None else []
            if not df:
                raise ValueError('At least one failure (df) is required.')
            fdf = self._as_frame(x_df, len(df), 'x')
            sdf = self._as_frame(x_ds, len(ds), 'x')
            if fdf is None:
                raise ValueError('x_df is required when df is given.')
            if ds and sdf is None:
                raise ValueError('x_ds is required when ds is given.')
            x_all = pd.concat([fdf] + ([sdf] if sdf is not None else []),
                              ignore_index=True)
            t = np.asarray(df + ds, dtype=float)
            event = np.array([1] * len(df) + [0] * len(ds), dtype=int)

        # --- named life-stress (aging) laws: transform the raw columns ---
        self._stress_spec = {}          # raw column -> law name
        self._stress_units = {}         # raw column -> 'C' | 'K'
        self._stress_terms = {}         # raw column -> [design column names]
        if stress_model:
            # Shorthand: a bare law name is allowed when there is exactly
            # one covariate column - the law then applies to that column.
            if isinstance(stress_model, str):
                if x_all.shape[1] != 1:
                    raise ValueError(
                        "stress_model as a bare law name is only allowed "
                        "with a single covariate column; pass a dict "
                        "{column: law_name} otherwise.")
                stress_model = {x_all.columns[0]: stress_model}
            if not isinstance(stress_model, dict):
                raise ValueError('stress_model must be a dict '
                                 '{column: law_name} or a bare law name.')
            su = stress_units or {}
            if not isinstance(su, dict):
                raise ValueError('stress_units must be a dict {column: '
                                 "'C'|'K'}.")
            for raw, law in stress_model.items():
                if law not in self._STRESS_LAWS:
                    raise ValueError(
                        f"stress_model['{raw}'] = '{law}' is not supported. "
                        f"Supported: {sorted(self._STRESS_LAWS)}")
                if raw not in x_all.columns:
                    raise ValueError(f"stress column '{raw}' not found in the "
                                     f"data.")
                unit = str(su.get(raw, 'C')).upper()
                if unit not in ('C', 'K'):
                    raise ValueError(f"stress_units['{raw}'] must be 'C' or "
                                     f"'K'.")
                self._stress_spec[raw] = law
                self._stress_units[raw] = unit
            x_all = self._stress_transform_frame(x_all, record=True)

        # raw covariate frame (pre one-hot), kept so the descriptive
        # estimators can split on an original column via by=.
        self._x_frame = x_all.reset_index(drop=True).copy()

        x_all = pd.get_dummies(x_all, drop_first=True)
        x_all = x_all.apply(pd.to_numeric, errors='coerce')
        if (feature_names is not None and not self._stress_spec
                and len(feature_names) == x_all.shape[1]):
            x_all.columns = list(feature_names)
        self.feature_names = list(x_all.columns)
        x_raw = x_all.to_numpy(dtype=float)

        if x_raw.shape[1] == 0:
            raise ValueError('No covariates found.')
        if not np.all(np.isfinite(x_raw)):
            raise ValueError('Covariates contain NaN/inf (or non-numeric '
                             'values that could not be encoded).')
        if not np.all(np.isfinite(t)) or np.any(t <= 0):
            raise ValueError('Durations must be finite and strictly positive.')
        if set(np.unique(event)) - {0, 1}:
            raise ValueError('Event indicator must be 0/1.')
        if event.sum() < 1:
            raise ValueError('At least one failure (event == 1) is required.')
        const_cols = [self.feature_names[i]
                      for i in np.where(x_raw.var(axis=0) == 0)[0]]
        if const_cols:
            stress_terms = {c for cols in self._stress_terms.values()
                            for c in cols}
            if set(const_cols) & stress_terms:
                bad = sorted(set(const_cols) & stress_terms)
                raise ValueError(
                    f'Life-stress term(s) {bad} are constant - the test has '
                    f'only one level of that stress. At least two distinct '
                    f'stress levels are needed to fit the aging law.')
            raise ValueError(f'Constant covariate column(s): {const_cols}. '
                             f'Remove them (the intercept already carries a '
                             f'constant term).')

        self._t = t
        self._event = event.astype(float)
        self._x_raw = x_raw
        self.n = int(len(t))
        self.n_events = int(event.sum())
        self.censoring = 'right-censored' if self.n_events < self.n else None

        # Result attributes - filled by fit().
        self.coef = None
        self.se_coef = None
        self.z_values = None
        self.p_values = None
        self.ci_lower = None
        self.ci_upper = None
        self.cov = None
        self._coef_cov = None
        self.params_ = None
        self.summary_ = None
        self.stress_params = None
        self.loglik = None
        self.loglik_null = None
        self.aic = None
        self.lr_stat = None
        self.lr_pvalue = None
        # weibull_aft only
        self.intercept = None
        self.sigma = None
        self.se_sigma = None
        self.beta = None
        # cox_ph only
        self.hazard_ratio = None
        self.hazard_ratio_ci = None
        self.baseline_cumhaz = None
        self.baseline_surv = None
        self.concordance = None
        # lazily filled caches for the survival-function bands
        self._cox_band = None
        self._boot_cache = {}

    # ------------------------------------------------------------------
    # named life-stress (aging) laws
    # ------------------------------------------------------------------
    def _stress_terms_for(self, law, raw, s, unit):
        """{design_column: values} for one named life-stress law applied to
        the raw stress values s (unit 'C' or 'K' for the thermal laws)."""
        s = np.asarray(s, dtype=float)
        kB = self._K_B_EV
        if self._STRESS_LAWS[law]['kelvin']:
            TK = s if unit == 'K' else s + 273.15
            if not np.all(np.isfinite(TK)) or np.any(TK <= 0):
                raise ValueError(f"stress '{raw}': temperature must be a "
                                 f"finite value > 0 K (stress_units="
                                 f"'{unit}').")
            out = {f'{raw}_invkT': 1.0 / (kB * TK)}
            if law == 'eyring':
                out[f'{raw}_lnT'] = np.log(TK)
            return out
        if law in ('inverse_power', 'coffin_manson'):
            if not np.all(np.isfinite(s)) or np.any(s <= 0):
                raise ValueError(f"stress '{raw}': values for '{law}' must be "
                                 f"finite and > 0.")
            return {f'{raw}_ln': np.log(s)}
        if not np.all(np.isfinite(s)):
            raise ValueError(f"stress '{raw}': values must be finite.")
        return {f'{raw}_lin': s}

    def _stress_transform_frame(self, frame, record=False):
        """Replace every raw stress column present in `frame` by its
        life-stress design columns. With record=True also stores the
        column names in self._stress_terms (used once, at ingest)."""
        out = frame.copy()
        for raw, law in self._stress_spec.items():
            if raw not in out.columns:
                continue
            s = pd.to_numeric(out.pop(raw), errors='coerce').to_numpy()
            terms = self._stress_terms_for(
                law, raw, s, self._stress_units.get(raw, 'C'))
            if record:
                self._stress_terms[raw] = list(terms)
            for name, vals in terms.items():
                out[name] = vals
        return out

    @staticmethod
    def _coerce_pred_frame(X):
        """dict / list-of-dicts / DataFrame -> DataFrame, accepting a plain
        dict of scalars as a single row."""
        if isinstance(X, dict):
            scalarish = all(np.ndim(v) == 0 for v in X.values())
            return pd.DataFrame([X] if scalarish else X).reset_index(drop=True)
        return pd.DataFrame(X).reset_index(drop=True)

    # ------------------------------------------------------------------
    # data helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _as_frame(x, n_expected, name_hint):
        """Coerce a covariate block (DataFrame / dict / 2D array-like) to a
        DataFrame and check its row count."""
        if x is None:
            return None
        if isinstance(x, pd.DataFrame):
            frame = x.reset_index(drop=True).copy()
        elif isinstance(x, dict):
            frame = pd.DataFrame(x)
        else:
            arr = np.asarray(x)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            if arr.ndim != 2:
                raise ValueError('Covariates must be 2D (rows = units).')
            frame = pd.DataFrame(
                arr, columns=[f'{name_hint}{i}' for i in range(arr.shape[1])])
        if len(frame) != n_expected:
            raise ValueError(f'Covariate block has {len(frame)} rows but '
                             f'{n_expected} were expected.')
        return frame

    def _prep_X(self, X):
        """Turn user-supplied covariate values for prediction into an
        (m, n_features) array in self.feature_names order. With a
        stress_model the raw stress columns are accepted and transformed
        exactly as at fit time."""
        if isinstance(X, (pd.DataFrame, dict)):
            frame = self._coerce_pred_frame(X)
            if self._stress_spec:
                have = [r for r in self._stress_spec if r in frame.columns]
                if have:
                    frame = self._stress_transform_frame(frame)
                missing = [r for r in self._stress_spec
                           if r not in have
                           and not set(self._stress_terms[r])
                           .issubset(frame.columns)]
                if missing:
                    raise ValueError(f'prediction input is missing stress '
                                     f'column(s) {missing} (give them in raw '
                                     f'units).')
            frame = pd.get_dummies(frame, drop_first=True)
            frame = frame.reindex(columns=self.feature_names, fill_value=0)
            return frame.to_numpy(dtype=float)
        arr = np.asarray(X, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.shape[1] != len(self.feature_names):
            raise ValueError(f'Expected {len(self.feature_names)} covariate '
                             f'columns, got {arr.shape[1]}.')
        return arr

    # ------------------------------------------------------------------
    # Weibull AFT: log-likelihood, score, observed information
    # ------------------------------------------------------------------
    # Parameters are packed as [coef (k), logsig]. With design row z_i
    # (including the intercept column when fit_intercept), y_i = ln t_i,
    # sigma = exp(logsig), standardized residual r_i = (y_i - z_i'coef)/sigma
    # and event indicator d_i:
    #
    #   ln L = sum_i [ d_i (r_i - logsig - y_i) - exp(r_i) ]
    #
    # (the -y_i term is parameter-free but kept so the reported loglik is a
    # genuine data log-likelihood). With dr_i/dcoef = -z_i/sigma and
    # dr_i/dlogsig = -r_i:
    #
    #   g_coef   = -(1/sigma) * sum_i (d_i - exp r_i) z_i
    #   g_logsig = -sum_i [ (d_i - exp r_i) r_i + d_i ]
    #
    #   H_coefcoef   = -(1/sigma^2) sum_i exp(r_i) z_i z_i'
    #   H_coeflogsig =  (1/sigma)   sum_i [ d_i - exp(r_i)(1 + r_i) ] z_i
    #   H_logsig^2   = -sum_i [ exp(r_i) r_i (r_i + 1) - r_i d_i ]
    #
    # observed information = -H.
    # exp(r) is the smallest-extreme-value hazard term; a standardized
    # log-time residual r anywhere near this cap only happens at parameter
    # values the line search / profile search is about to reject, so
    # capping it keeps the objective and its derivatives finite (no
    # overflow to +/-inf) without affecting any step that is actually
    # accepted - at the optimum r is O(1).
    _R_CAP = 30.0

    @staticmethod
    def _aft_negll(params, Z, y, d):
        k = Z.shape[1]
        coef, logsig = params[:k], params[k]
        # clamp only the runaway excursions an optimiser may probe; at any
        # real optimum logsig is O(1) so this bound never binds.
        sig = np.exp(np.clip(logsig, -Regression._R_CAP, Regression._R_CAP))
        r = np.minimum((y - Z @ coef) / sig, Regression._R_CAP)
        return -float(np.sum(d * (r - logsig - y) - np.exp(r)))

    @staticmethod
    def _aft_neggrad(params, Z, y, d):
        k = Z.shape[1]
        coef, logsig = params[:k], params[k]
        # clamp only the runaway excursions an optimiser may probe; at any
        # real optimum logsig is O(1) so this bound never binds.
        sig = np.exp(np.clip(logsig, -Regression._R_CAP, Regression._R_CAP))
        r = np.minimum((y - Z @ coef) / sig, Regression._R_CAP)
        er = np.exp(r)
        g_coef = -(1.0 / sig) * (Z.T @ (d - er))
        g_logsig = -float(np.sum((d - er) * r + d))
        return -np.concatenate([g_coef, [g_logsig]])

    @staticmethod
    def _aft_info(params, Z, y, d):
        k = Z.shape[1]
        coef, logsig = params[:k], params[k]
        # clamp only the runaway excursions an optimiser may probe; at any
        # real optimum logsig is O(1) so this bound never binds.
        sig = np.exp(np.clip(logsig, -Regression._R_CAP, Regression._R_CAP))
        r = np.minimum((y - Z @ coef) / sig, Regression._R_CAP)
        er = np.exp(r)
        h_cc = -(1.0 / sig ** 2) * (Z.T @ (er[:, None] * Z))
        h_cl = (1.0 / sig) * (Z.T @ (d - er * (1.0 + r)))
        h_ll = -float(np.sum(er * r * (r + 1.0) - r * d))
        H = np.empty((k + 1, k + 1))
        H[:k, :k] = h_cc
        H[:k, k] = h_cl
        H[k, :k] = h_cl
        H[k, k] = h_ll
        return -H

    def _aft_fit(self, Z, y, d, start=None):
        """Maximize the AFT log-likelihood. Without a start: OLS over the
        events for the location, residual spread for logsig, a bounded
        quasi-Newton climb, then exact Newton polish. With a start (used
        for the original-scale re-fit after standardized fitting): Newton
        with backtracking only."""
        k = Z.shape[1]
        if start is None:
            ev = d.astype(bool)
            coef0, *_ = np.linalg.lstsq(Z[ev], y[ev], rcond=None)
            resid = y[ev] - Z[ev] @ coef0
            s0 = np.std(resid) if resid.size > 1 and np.std(resid) > 0 else 1.0
            p = np.concatenate([coef0, [np.log(s0)]])
            res = optimize.minimize(
                self._aft_negll, p, args=(Z, y, d), jac=self._aft_neggrad,
                method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-12,
                                            'maxiter': self.max_iter * 20})
            p = np.asarray(res.x, dtype=float)
        else:
            p = np.asarray(start, dtype=float)

        ll = -self._aft_negll(p, Z, y, d)
        for _ in range(self.max_iter):
            grad = -self._aft_neggrad(p, Z, y, d)
            info = self._aft_info(p, Z, y, d)
            try:
                step = np.linalg.solve(info, grad)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(info, grad, rcond=None)[0]
            alpha = 1.0
            for _bt in range(30):
                cand = p + alpha * step
                ll_new = -self._aft_negll(cand, Z, y, d)
                if np.isfinite(ll_new) and ll_new >= ll - 1e-12:
                    break
                alpha *= 0.5
            else:
                break
            p, ll = cand, ll_new
            if np.max(np.abs(grad)) < self.tol:
                break
        info = self._aft_info(p, Z, y, d)
        try:
            cov = np.linalg.inv(info)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(info)
        return p, cov, ll

    # ------------------------------------------------------------------
    # Cox PH: partial log-likelihood, score, observed information
    # ------------------------------------------------------------------
    # With eta_i = x_i'beta, w_i = exp(eta_i), ordered distinct event times
    # tau_k, death set D_k (size d_k) and risk set R_k = {i : t_i >= tau_k},
    # and risk-set moments
    #   S0 = sum_{R_k} w_j,  S1 = sum_{R_k} w_j x_j,  S2 = sum_{R_k} w_j x_j x_j'
    # Breslow:
    #   ln L = sum_k [ sum_{D_k} eta_i - d_k ln S0 ]
    #   U    = sum_k [ sum_{D_k} x_i - d_k S1/S0 ]
    #   I    = sum_k d_k [ S2/S0 - (S1/S0)(S1/S0)' ]
    # Efron replaces, for l = 0..d_k-1, S0/S1/S2 by
    #   S0 - (l/d_k) SD0, ...  with SD* the same moments over D_k only,
    # accounting for the tied deaths gradually leaving the risk set.
    # Iterating distinct times from largest to smallest makes S0/S1/S2 plain
    # running sums. I is a sum of covariance matrices, hence PSD.
    @staticmethod
    def _cox_stats(coef, X, t, d, ties):
        n, p = X.shape
        eta = X @ coef
        # Every risk-set term is a ratio of exp(eta) values, so a constant
        # shift of eta cancels exactly: subtract max(eta) before exp() to
        # keep w in (0, 1] and never overflow, whatever the coef scale.
        eta = eta - np.max(eta)
        w = np.exp(eta)
        order = np.argsort(t, kind='mergesort')
        ts, ds_, Xs, ws, etas = t[order], d[order], X[order], w[order], eta[order]
        ll = 0.0
        grad = np.zeros(p)
        info = np.zeros((p, p))
        rs0 = 0.0
        rs1 = np.zeros(p)
        rs2 = np.zeros((p, p))
        uniq = np.unique(ts)
        for tau in uniq[::-1]:
            lo = np.searchsorted(ts, tau, side='left')
            hi = np.searchsorted(ts, tau, side='right')
            Xg, wg = Xs[lo:hi], ws[lo:hi]
            rs0 += wg.sum()
            rs1 += wg @ Xg
            rs2 += (Xg * wg[:, None]).T @ Xg
            dmask = ds_[lo:hi].astype(bool)
            dk = int(dmask.sum())
            if dk == 0:
                continue
            Xd, wd = Xg[dmask], wg[dmask]
            grad += Xd.sum(axis=0)
            ll += float(etas[lo:hi][dmask].sum())
            if ties == 'breslow':
                ll -= dk * np.log(rs0)
                grad -= dk * (rs1 / rs0)
                info += dk * (rs2 / rs0 - np.outer(rs1, rs1) / rs0 ** 2)
            else:
                dw0 = wd.sum()
                dw1 = wd @ Xd
                dw2 = (Xd * wd[:, None]).T @ Xd
                for l in range(dk):
                    f = l / dk
                    den = rs0 - f * dw0
                    num1 = rs1 - f * dw1
                    num2 = rs2 - f * dw2
                    ll -= np.log(den)
                    grad -= num1 / den
                    info += num2 / den - np.outer(num1, num1) / den ** 2
        return ll, grad, info

    def _cox_fit(self, X, t, d, start=None):
        p = X.shape[1]
        coef = np.zeros(p) if start is None else np.asarray(start, dtype=float)
        ll, grad, info = self._cox_stats(coef, X, t, d, self.ties)
        for _ in range(self.max_iter):
            try:
                step = np.linalg.solve(info, grad)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(info, grad, rcond=None)[0]
            alpha = 1.0
            for _bt in range(30):
                cand = coef + alpha * step
                ll_new, grad_new, info_new = self._cox_stats(
                    cand, X, t, d, self.ties)
                if np.isfinite(ll_new) and ll_new >= ll - 1e-12:
                    break
                alpha *= 0.5
            else:
                break
            coef, ll, grad, info = cand, ll_new, grad_new, info_new
            if np.max(np.abs(step)) < self.tol:
                break
        # A partial likelihood that stays monotone in some direction of
        # beta (a covariate, or a combination, that perfectly orders the
        # event times - "monotone/complete separation") has its maximum
        # only at infinity: the iteration then just drifts outward with a
        # near-singular information matrix. Flag it instead of returning a
        # meaningless huge coefficient.
        if (not np.all(np.isfinite(coef))
                or np.max(np.abs(coef)) > self._COX_DIVERGENCE
                or not np.all(np.isfinite(info))
                or np.linalg.cond(info) > 1e12):
            raise ValueError(
                'Cox partial likelihood did not converge - a covariate (or '
                'linear combination) almost perfectly separates the event '
                'ordering, so its coefficient is unbounded. Drop or combine '
                'that covariate, or collect more data.')
        try:
            cov = np.linalg.inv(info)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(info)
        return coef, cov, ll

    def _breslow_baseline(self, eta, t, d):
        """Breslow cumulative baseline hazard at the distinct event times:
        H0(tau_k) = sum_{m<=k} d_m / sum_{j: t_j >= tau_m} exp(eta_j)."""
        order = np.argsort(t, kind='mergesort')
        ts, ds_, we = t[order], d[order].astype(bool), np.exp(eta[order])
        ev_times = np.unique(ts[ds_])
        H0 = np.empty(len(ev_times))
        cum = 0.0
        for i, tau in enumerate(ev_times):
            lo = np.searchsorted(ts, tau, side='left')
            risk = we[lo:].sum()
            deaths = int(np.count_nonzero((ts == tau) & ds_))
            cum += deaths / risk
            H0[i] = cum
        return ev_times, H0

    @staticmethod
    def _concordance(risk, t, d):
        """Harrell's concordance index. A higher risk score is expected to
        go with a shorter lifetime; ties in the score count as half."""
        num = den = 0.0
        for i in np.where(d == 1)[0]:
            later = t > t[i]
            m = int(np.count_nonzero(later))
            if m == 0:
                continue
            den += m
            num += np.count_nonzero(risk[i] > risk[later])
            num += 0.5 * np.count_nonzero(risk[i] == risk[later])
        return num / den if den > 0 else np.nan

    # ------------------------------------------------------------------
    # profile-likelihood bounds (bounds='lrb')
    # ------------------------------------------------------------------
    def _profile_ci(self, idx, params, ll_hat, restricted_maxll, se):
        """Invert the likelihood-ratio test for parameter idx: the interval
        is the set of v with 2 (ll_hat - max_{rest} ll(v)) <= chi2(cl, 1)."""
        thresh = chi2.ppf(self.cl if self.bounds_type == '2s' else 2 * self.cl - 1, 1)
        v_hat = params[idx]
        span = se if np.isfinite(se) and se > 0 else max(abs(v_hat), 1.0)

        def g(v):
            return 2.0 * (ll_hat - restricted_maxll(idx, v, params)) - thresh

        roots = []
        for sign in (-1.0, 1.0):
            step, v_out, it = span, v_hat + sign * span, 0
            while g(v_out) < 0 and it < 60:
                step *= 1.6
                v_out += sign * step
                it += 1
            try:
                roots.append(optimize.brentq(g, min(v_hat, v_out),
                                             max(v_hat, v_out), xtol=1e-7,
                                             maxiter=200))
            except (ValueError, RuntimeError):
                roots.append(np.nan)
        return roots[0], roots[1]

    def _aft_restricted_maxll(self, idx, v, params):
        free = np.ones(len(params), dtype=bool)
        free[idx] = False
        x0 = params[free]

        def neg(z):
            full = params.copy()
            full[free] = z
            full[idx] = v
            return self._aft_negll(full, self._Z, self._y, self._event)

        def negjac(z):
            full = params.copy()
            full[free] = z
            full[idx] = v
            return self._aft_neggrad(full, self._Z, self._y, self._event)[free]

        res = optimize.minimize(neg, x0, jac=negjac, method='L-BFGS-B',
                                options={'ftol': 1e-14, 'gtol': 1e-10})
        return -res.fun

    def _cox_restricted_maxll(self, idx, v, params):
        """Maximize the Cox partial log-likelihood over every coefficient
        except idx (held at v), by Newton steps on the free sub-vector
        using the analytic score/information blocks."""
        free = np.ones(len(params), dtype=bool)
        free[idx] = False
        b = params.copy()
        b[idx] = v
        ll = self._cox_stats(b, self._x_raw, self._t, self._event, self.ties)[0]
        for _ in range(self.max_iter):
            _, g, I = self._cox_stats(b, self._x_raw, self._t,
                                      self._event, self.ties)
            gf = g[free]
            if np.max(np.abs(gf)) < self.tol:
                break
            If = I[np.ix_(free, free)]
            try:
                step = np.linalg.solve(If, gf)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(If, gf, rcond=None)[0]
            alpha = 1.0
            for _bt in range(30):
                cand = b.copy()
                cand[free] = b[free] + alpha * step
                ll_new = self._cox_stats(cand, self._x_raw, self._t,
                                         self._event, self.ties)[0]
                if np.isfinite(ll_new) and ll_new >= ll - 1e-12:
                    break
                alpha *= 0.5
            else:
                break
            b, ll = cand, ll_new
        return ll

    def _wald_ci(self, coef, se):
        coef = np.asarray(coef, dtype=float)
        se = np.asarray(se, dtype=float)
        if self.bounds_type == '2s':
            k = norm.ppf(1 - (1 - self.cl) / 2)
            return coef - k * se, coef + k * se
        k = norm.ppf(self.cl)
        if self.bounds_type == '1sl':
            return coef - k * se, np.full_like(coef, np.inf)
        return np.full_like(coef, -np.inf), coef + k * se

    # ------------------------------------------------------------------
    # fit
    # ------------------------------------------------------------------
    def fit(self):
        """Estimate the model and fill the result attributes. Returns self,
        so Regression(...).fit().summary() works."""
        x_raw = self._x_raw
        n, p = x_raw.shape
        d = self._event
        self._cox_band = None
        self._boot_cache = {}
        if self.standardize:
            self._x_mean = x_raw.mean(axis=0)
            self._x_std = np.where(x_raw.std(axis=0, ddof=0) == 0, 1.0,
                                   x_raw.std(axis=0, ddof=0))
        else:
            self._x_mean = np.zeros(p)
            self._x_std = np.ones(p)
        xs = (x_raw - self._x_mean) / self._x_std

        if self.model == 'weibull_aft':
            self._fit_weibull_aft(x_raw, xs, d)
        else:
            self._fit_cox_ph(x_raw, xs, d)

        self._build_stress_params()
        self._build_summary()
        if self.show or self.save:
            self.plot()
        return self

    def _fit_weibull_aft(self, x_raw, xs, d):
        n, p = x_raw.shape
        y = np.log(self._t)
        Zs = np.column_stack([np.ones(n), xs]) if self.fit_intercept else xs.copy()
        Zo = np.column_stack([np.ones(n), x_raw]) if self.fit_intercept else x_raw.copy()
        k = Zs.shape[1]

        params_s, _, _ = self._aft_fit(Zs, y, d)
        # map the standardized solution to the original scale, then polish
        b = params_s[:k]
        if self.fit_intercept:
            b_o = np.empty(k)
            b_o[0] = b[0] - np.sum(b[1:] * self._x_mean / self._x_std)
            b_o[1:] = b[1:] / self._x_std
        else:
            b_o = b / self._x_std
        start = np.concatenate([b_o, [params_s[k]]])
        params, cov, ll = self._aft_fit(Zo, y, d, start=start)

        self._Z, self._y = Zo, y  # for profile bounds
        self._aft_params = params

        se_all = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
        logsig = params[k]
        self.sigma = float(np.exp(logsig))
        self.se_sigma = float(self.sigma * se_all[k])
        self.beta = float(1.0 / self.sigma)

        names = list(self.feature_names)
        if self.fit_intercept:
            self.intercept = float(params[0])
            coef = params[1:k]
            se = se_all[1:k]
            idx_map = list(range(1, k))
            names = names + ['Intercept']
            coef_row = np.append(coef, params[0])
            se_row = np.append(se, se_all[0])
            idx_row = idx_map + [0]
        else:
            self.intercept = 0.0
            coef = params[:k]
            se = se_all[:k]
            coef_row, se_row, idx_row = coef, se, list(range(k))

        self.coef = np.asarray(coef, dtype=float)
        self.se_coef = np.asarray(se, dtype=float)
        self.cov = cov
        s0 = 1 if self.fit_intercept else 0
        self._coef_cov = cov[s0:s0 + len(coef), s0:s0 + len(coef)]
        self.params_ = pd.Series(coef_row, index=names)

        if self.bounds == 'lrb':
            lo_row, hi_row = self._profile_bounds(
                idx_row, params, ll, self._aft_restricted_maxll, se_row)
        elif self.bounds == 'fb':
            lo_row, hi_row = self._wald_ci(coef_row, se_row)
        elif self.bounds in ('npbb', 'pbb'):
            rows = self._boot_coef_rows(self.bounds)
            a = (1.0 - self.cl) / 2.0
            lo_row = np.quantile(rows, a, axis=0)
            hi_row = np.quantile(rows, 1.0 - a, axis=0)
        else:
            lo_row = hi_row = np.full(len(coef_row), np.nan)

        self._row_names = names
        self._coef_row = coef_row
        self._se_row = se_row
        self._lo_row = lo_row
        self._hi_row = hi_row
        self.ci_lower = np.asarray(lo_row[:len(self.coef)], dtype=float)
        self.ci_upper = np.asarray(hi_row[:len(self.coef)], dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.z_values = self.coef / self.se_coef
        self.p_values = 2.0 * norm.sf(np.abs(self.z_values))

        # Harrell's concordance: a longer predicted lifetime (larger x'coef)
        # is a lower risk, so rank by the negative linear predictor.
        self.concordance = float(self._concordance(-(x_raw @ coef),
                                                   self._t, d))

        # null model: intercept + scale only
        _, _, ll0 = self._aft_fit(np.ones((n, 1)), y, d)
        self.loglik = float(ll)
        self.loglik_null = float(ll0)
        n_par = p + (1 if self.fit_intercept else 0) + 1
        self.aic = float(2 * n_par - 2 * ll)
        self.lr_stat = float(2 * (ll - ll0))
        self.lr_pvalue = float(chi2.sf(self.lr_stat, p))

    def _fit_cox_ph(self, x_raw, xs, d):
        n, p = x_raw.shape
        coef_s, _, _ = self._cox_fit(xs, self._t, d)
        coef, cov, ll = self._cox_fit(x_raw, self._t, d,
                                      start=coef_s / self._x_std)
        self._cox_params = coef
        se = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))

        self.coef = np.asarray(coef, dtype=float)
        self.se_coef = np.asarray(se, dtype=float)
        self.cov = cov
        self._coef_cov = cov
        self.params_ = pd.Series(self.coef, index=list(self.feature_names))
        self.intercept = None

        # Baseline is taken at the mean covariate vector (not at x = 0,
        # which is usually far outside the data): baseline_cumhaz /
        # baseline_surv then describe an "average" unit, and predictions
        # multiply by exp((x - xbar)'beta). Computed before the bounds
        # block so parametric-bootstrap bounds can simulate from it.
        eta_c = (x_raw - self._x_mean) @ coef
        self.baseline_cumhaz = self._breslow_baseline(eta_c, self._t, d)
        bt, bH = self.baseline_cumhaz
        self.baseline_surv = (bt, np.exp(-bH))

        if self.bounds == 'lrb':
            lo, hi = self._profile_bounds(list(range(p)), coef, ll,
                                          self._cox_restricted_maxll, se)
        elif self.bounds == 'fb':
            lo, hi = self._wald_ci(coef, se)
        elif self.bounds in ('npbb', 'pbb'):
            rows = self._boot_coef_rows(self.bounds)
            a = (1.0 - self.cl) / 2.0
            lo = np.quantile(rows, a, axis=0)
            hi = np.quantile(rows, 1.0 - a, axis=0)
        else:
            lo = hi = np.full(p, np.nan)
        self.ci_lower, self.ci_upper = np.asarray(lo), np.asarray(hi)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.z_values = self.coef / self.se_coef
        self.p_values = 2.0 * norm.sf(np.abs(self.z_values))
        self._row_names = list(self.feature_names)
        self._coef_row, self._se_row = self.coef, self.se_coef
        self._lo_row, self._hi_row = self.ci_lower, self.ci_upper

        self.hazard_ratio = np.exp(self.coef)
        self.hazard_ratio_ci = (np.exp(self.ci_lower), np.exp(self.ci_upper))
        self.concordance = float(self._concordance(x_raw @ coef, self._t, d))

        ll0 = self._cox_stats(np.zeros(p), x_raw, self._t, d, self.ties)[0]
        self.loglik = float(ll)
        self.loglik_null = float(ll0)
        self.aic = float(2 * p - 2 * ll)
        self.lr_stat = float(2 * (ll - ll0))
        self.lr_pvalue = float(chi2.sf(self.lr_stat, p))

    def _profile_bounds(self, idx_list, params, ll_hat, restricted_maxll, se):
        lo = np.empty(len(idx_list))
        hi = np.empty(len(idx_list))
        for j, idx in enumerate(idx_list):
            try:
                a, b = self._profile_ci(idx, params, ll_hat,
                                        restricted_maxll, se[j])
            except Exception:
                a, b = np.nan, np.nan
            if not np.isfinite(a) or not np.isfinite(b):
                # fall back to the Wald interval for this parameter
                wl, wh = self._wald_ci(np.array([params[idx]]),
                                       np.array([se[j]]))
                a = a if np.isfinite(a) else float(wl[0])
                b = b if np.isfinite(b) else float(wh[0])
            lo[j], hi[j] = a, b
        return lo, hi

    # ==================================================================
    # confidence bands for the predicted survival function
    #
    # The band is built on the complementary-log-log scale
    #   eta(t | x) = ln(-ln S(t | x))
    # (for Weibull AFT: eta = (ln t - x'theta) / sigma; for Cox:
    #  eta = ln H0_mean(t) + (x - xbar)'beta), because eta is unbounded so
    # the back-transformed band always stays inside (0, 1) and its
    # small-sample coverage is better than a band drawn straight on S.
    #
    # bounds='fb'  -> delta method: Var(eta_hat) = g' V g with g the
    #                 gradient of eta w.r.t. the parameters and V their
    #                 covariance (Cox additionally carries the Breslow
    #                 baseline variance and its coupling to beta_hat -
    #                 Klein & Moeschberger, "Survival Analysis", ch. 8.6;
    #                 Kalbfleisch & Prentice, ch. 4.3).
    # bounds='lrb' -> Weibull AFT only: genuine profile-likelihood band
    #                 (pin eta(t) at a trial value, maximise over the
    #                 remaining parameters, invert the likelihood-ratio
    #                 statistic). For Cox there is no profile-likelihood
    #                 band for the Breslow survival curve, so 'lrb' uses
    #                 the same delta band as 'fb' there.
    # simultaneous  -> band valid over the whole observed time range at
    #                 once. The pointwise z-quantile is replaced by the
    #                 (1-alpha) quantile of sup_t |Z(t)| for the mean-zero
    #                 Gaussian process Z of eta_hat(t), obtained by Monte
    #                 Carlo from its estimated correlation matrix (an
    #                 equal-precision-type band; same construction for AFT
    #                 and Cox). The edges are then tightened to be
    #                 monotone in t without losing coverage, since S is.
    # ==================================================================
    def _resolve_band_opts(self, cl, bounds):
        cl = self.cl if cl is None else float(cl)
        if not 0.0 < cl < 1.0:
            raise ValueError('cl must be in (0, 1).')
        method = bounds if bounds is not None else (self.bounds or 'fb')
        if method not in ('fb', 'lrb', 'npbb', 'pbb'):
            raise ValueError("the survival band needs bounds 'fb', 'lrb', "
                             "'npbb' or 'pbb'.")
        return cl, method

    # ---- Weibull AFT -------------------------------------------------
    def _aft_eta_grad(self, xrow, times):
        """clog-log transform eta(t) = (ln t - lp) / sigma for one covariate
        row and its gradient w.r.t. the packed parameter vector
        [intercept?, coef_1..p, logsig] (the order of self.cov)."""
        xrow = np.asarray(xrow, dtype=float).ravel()
        lp = (self.intercept if self.fit_intercept else 0.0) + xrow @ self.coef
        eta = (np.log(times) - lp) / self.sigma
        p = len(self.coef)
        G = np.empty((len(times), self.cov.shape[0]))
        c = 0
        if self.fit_intercept:
            G[:, 0] = -1.0 / self.sigma
            c = 1
        G[:, c:c + p] = -xrow[None, :] / self.sigma
        G[:, c + p] = -eta                       # d eta / d logsig
        return eta, G

    def _aft_delta_band(self, xrow, times, cl):
        eta, G = self._aft_eta_grad(xrow, times)
        var = np.einsum('ti,ij,tj->t', G, self.cov, G)
        se = np.sqrt(np.clip(var, 0.0, np.inf))
        z = norm.ppf(1.0 - (1.0 - cl) / 2.0)
        # far into extrapolation eta +/- z*se can exceed the exp() range;
        # exp(-inf) = 0 is the intended clamped bound, so ignore the overflow.
        with np.errstate(over='ignore'):
            S = np.exp(-np.exp(eta))
            lower = np.exp(-np.exp(eta + z * se))
            upper = np.exp(-np.exp(eta - z * se))
        return S, lower, upper, eta, se, G

    def _aft_profile_band(self, xrow, times, cl):
        """Profile-likelihood pointwise band for the AFT survival curve.
        Needs an intercept: it is the parameter substituted out to pin
        eta = (ln t - lp) / sigma at each trial value. Evaluated on a
        coarse log-spaced sub-grid and linearly interpolated in ln t (eta
        is close to linear there)."""
        if not self.fit_intercept:
            raise ValueError("bounds='lrb' survival band needs "
                             "fit_intercept=True; use bounds='fb'.")
        xrow = np.asarray(xrow, dtype=float).ravel()
        y, d, Z = self._y, self._event, self._Z
        p = len(self.coef)
        thr = chi2.ppf(cl, 1)
        ll_max = -self._aft_negll(self._aft_params, Z, y, d)
        free0 = np.concatenate([self._aft_params[1:1 + p],
                                [self._aft_params[-1]]])

        def restr(free, eta0, tv):
            coef_f, logsig_f = free[:p], free[p]
            sig_f = np.exp(logsig_f)
            b0 = np.log(tv) - xrow @ coef_f - sig_f * eta0
            full = np.concatenate([[b0], coef_f, [logsig_f]])
            neg = self._aft_negll(full, Z, y, d)
            gf = self._aft_neggrad(full, Z, y, d)
            g_b0 = gf[0]
            return neg, np.concatenate([gf[1:1 + p] - g_b0 * xrow,
                                        [gf[p + 1] - g_b0 * sig_f * eta0]])

        def prof_ll(eta0, tv):
            r = optimize.minimize(restr, free0, args=(eta0, tv), jac=True,
                                  method='L-BFGS-B',
                                  options={'ftol': 1e-13, 'gtol': 1e-9})
            return -r.fun

        m = int(min(len(times), 24))
        gt = np.geomspace(max(times.min(), times.max() * 1e-6),
                          times.max(), m)
        lp = self.intercept + xrow @ self.coef
        e_lo = np.empty(m)
        e_hi = np.empty(m)
        for i, tv in enumerate(gt):
            e_hat = (np.log(tv) - lp) / self.sigma

            def g(e0):
                return 2.0 * (ll_max - prof_ll(e0, tv)) - thr

            for sgn, store in ((-1.0, e_lo), (1.0, e_hi)):
                step, v, it = 1.0, e_hat + sgn, 0
                while g(v) < 0 and it < 40:
                    step *= 1.6
                    v += sgn * step
                    it += 1
                try:
                    store[i] = optimize.brentq(g, min(e_hat, v), max(e_hat, v),
                                               xtol=1e-5, maxiter=100)
                except (ValueError, RuntimeError):
                    store[i] = e_hat + sgn * 6.0
        lt, lg = np.log(times), np.log(gt)
        eta_lo = np.interp(lt, lg, e_lo)
        eta_hi = np.interp(lt, lg, e_hi)
        eta_hat = (lt - lp) / self.sigma
        with np.errstate(over='ignore'):          # exp(-inf) = 0 in the tails
            return (np.exp(-np.exp(eta_hat)), np.exp(-np.exp(eta_hi)),
                    np.exp(-np.exp(eta_lo)), eta_hat, eta_lo, eta_hi)

    # ---- Cox PH ----------------------------------------------------------
    def _cox_band_pieces(self):
        """Cumulative sums over the distinct event times needed for the
        Breslow-survival variance: S0_i (risk-set sum of exp of centred
        covariates), the risk-set mean centred covariate Ebar_i, and the
        running sums  sum d_i / S0_i^2 ,  sum d_i / S0_i ,
        sum d_i Ebar_i / S0_i ."""
        if self._cox_band is not None:
            return self._cox_band
        X = self._x_raw - self._x_mean
        w = np.exp(X @ self.coef)
        t, dmask = self._t, self._event.astype(bool)
        order = np.argsort(t, kind='mergesort')
        ts, Xs, ws, dm = t[order], X[order], w[order], dmask[order]
        ev = np.unique(ts[dm])
        D, p = len(ev), X.shape[1]
        S0 = np.empty(D)
        Ebar = np.empty((D, p))
        dk = np.empty(D)
        for i, tau in enumerate(ev):
            lo = np.searchsorted(ts, tau, side='left')
            S0[i] = ws[lo:].sum()
            Ebar[i] = (ws[lo:] @ Xs[lo:]) / S0[i]
            dk[i] = np.count_nonzero((ts == tau) & dm)
        self._cox_band = {
            'ev': ev,
            'cum_q1': np.cumsum(dk / S0 ** 2),
            'cum_a': np.cumsum(dk / S0),
            'cum_b': np.cumsum((dk / S0)[:, None] * Ebar, axis=0),
        }
        return self._cox_band

    def _cox_delta_band(self, xrow, times, cl):
        pc = self._cox_band_pieces()
        xc = np.asarray(xrow, dtype=float).ravel() - self._x_mean
        theta = float(np.exp(xc @ self.coef))
        z = norm.ppf(1.0 - (1.0 - cl) / 2.0)
        kk = np.searchsorted(pc['ev'], times, side='right')
        has = kk > 0
        idx = np.clip(kk - 1, 0, len(pc['ev']) - 1)
        H = pc['cum_a'][idx] * theta
        Q1 = pc['cum_q1'][idx]
        Q3 = xc[None, :] * pc['cum_a'][idx][:, None] - pc['cum_b'][idx]
        varH = theta ** 2 * (Q1 + np.einsum('ti,ij,tj->t', Q3, self.cov, Q3))
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            se_eta = np.sqrt(np.clip(varH, 0.0, np.inf)) / H
            S = np.where(has, np.exp(-H), 1.0)
            lo = np.where(has, np.exp(-H * np.exp(z * se_eta)), 1.0)
            hi = np.where(has, np.exp(-H * np.exp(-z * se_eta)), 1.0)
        eta = np.where(has & (H > 0), np.log(np.where(H > 0, H, 1.0)), -np.inf)
        return S, lo, hi, eta, np.where(has, se_eta, np.nan)

    # ---- simultaneous critical value -------------------------------
    @staticmethod
    def _gauss_sup_crit(corr, cl, n_sim=4000):
        """(1-alpha) quantile of sup_t |Z(t)| for a mean-zero Gaussian
        process with the given correlation matrix - the simultaneous
        multiplier for the parametric (AFT) band."""
        vals, vecs = np.linalg.eigh(corr)
        A = vecs * np.sqrt(np.clip(vals, 0.0, None))
        rng = np.random.default_rng(12345)
        draws = A @ rng.standard_normal((corr.shape[0], n_sim))
        return float(np.quantile(np.max(np.abs(draws), axis=0), cl))

    @staticmethod
    def _enforce_monotone_band(lo, hi):
        """Tighten a simultaneous band's edges to be non-increasing in t
        without losing simultaneous coverage: the true S is non-increasing,
        so an upper edge can be replaced by its running minimum from the
        left and a lower edge by its running maximum from the right. NaN
        entries (outside the band range) are left as they are."""
        lo = np.asarray(lo, dtype=float).copy()
        hi = np.asarray(hi, dtype=float).copy()
        fh = np.where(np.isfinite(hi))[0]
        if fh.size:
            hi[fh] = np.minimum.accumulate(hi[fh])
        fl = np.where(np.isfinite(lo))[0]
        if fl.size:
            lo[fl] = np.maximum.accumulate(lo[fl][::-1])[::-1]
        return lo, hi

    def _cox_eta_cov(self, xrow, times):
        """Covariance matrix of eta(t | x) = ln H(t | x) over `times` for a
        Cox profile. Baseline part: cov = Q1(min(s, t)) (independent
        increments); coefficient part: Q3(s)' cov(beta) Q3(t); both divided
        by H0(s) H0(t) so theta cancels. Same pieces the pointwise band
        uses, so the diagonal reproduces its se(eta)^2."""
        pc = self._cox_band_pieces()
        ev = pc['ev']
        xc = np.asarray(xrow, dtype=float).ravel() - self._x_mean
        kk = np.searchsorted(ev, times, side='right')
        idx = np.clip(kk - 1, 0, len(ev) - 1)
        a = pc['cum_a'][idx]                       # H0(t)
        Q1 = pc['cum_q1'][idx]
        Q3 = xc[None, :] * a[:, None] - pc['cum_b'][idx]
        with np.errstate(divide='ignore', invalid='ignore'):
            cov_eta = (np.minimum.outer(Q1, Q1)
                       + Q3 @ self.cov @ Q3.T) / np.outer(a, a)
        bad = kk == 0
        cov_eta[~np.isfinite(cov_eta)] = 0.0
        cov_eta[bad, :] = 0.0
        cov_eta[:, bad] = 0.0
        return cov_eta

    def _survival_band(self, xrow, times, cl, method, simultaneous):
        """Returns (S, lower, upper, lower_sim, upper_sim) for one profile,
        each an array over `times`. *_sim is None unless simultaneous."""
        t_lo, t_hi = float(self._t.min()), float(self._t.max())
        in_range = (times >= t_lo) & (times <= t_hi)
        sim_lo = sim_hi = None
        if method in ('npbb', 'pbb'):
            # bootstrap band: pointwise percentiles of the resampled
            # survival curves. No simultaneous variant.
            S, lower, upper = self._boot_band(xrow, times, cl, method)
            return S, lower, upper, None, None
        if self.model == 'weibull_aft':
            if method == 'lrb':
                S, lower, upper, eta, e_lo, e_hi = self._aft_profile_band(
                    xrow, times, cl)
            else:
                S, lower, upper, eta, se, G = self._aft_delta_band(
                    xrow, times, cl)
            if simultaneous:
                _, _, _, eta_d, se_d, G = self._aft_delta_band(xrow, times, cl)
                sel = in_range & np.isfinite(se_d) & (se_d > 0)
                if sel.sum() >= 2:
                    Gs = G[sel]
                    cov_eta = Gs @ self.cov @ Gs.T
                    dd = np.sqrt(np.clip(np.diag(cov_eta), 1e-300, None))
                    corr = cov_eta / np.outer(dd, dd)
                    k = self._gauss_sup_crit(corr, cl)
                    with np.errstate(over='ignore'):
                        sim_lo = np.exp(-np.exp(eta_d + k * se_d))
                        sim_hi = np.exp(-np.exp(eta_d - k * se_d))
                    sim_lo, sim_hi = self._enforce_monotone_band(sim_lo, sim_hi)
        else:
            S, lower, upper, eta, se = self._cox_delta_band(xrow, times, cl)
            if simultaneous:
                # Simultaneous band as the supremum of the estimated
                # Gaussian process for eta(t | x) = ln H(t | x), same
                # construction as the AFT branch. (A Hall-Wellner Brownian-
                # bridge multiplier is not used here: for a covariate-
                # adjusted curve the pointwise variance of eta is not
                # monotone in t, so its bridge time-transform is invalid.)
                ev = self._cox_band_pieces()['ev']
                band_ok = (np.isfinite(eta) & np.isfinite(se) & (se > 0)
                           & (times >= ev[0]) & (times <= ev[-1]))
                if band_ok.sum() >= 2:
                    cov_eta = self._cox_eta_cov(xrow, times)
                    sel = band_ok & (np.diag(cov_eta) > 0)
                    Cs = cov_eta[np.ix_(sel, sel)]
                    dd = np.sqrt(np.clip(np.diag(Cs), 1e-300, None))
                    k = self._gauss_sup_crit(Cs / np.outer(dd, dd), cl)
                    with np.errstate(over='ignore'):
                        sim_lo = np.where(band_ok,
                                          np.exp(-np.exp(eta + k * se)), np.nan)
                        sim_hi = np.where(band_ok,
                                          np.exp(-np.exp(eta - k * se)), np.nan)
                    sim_lo, sim_hi = self._enforce_monotone_band(sim_lo, sim_hi)
        return S, lower, upper, sim_lo, sim_hi

    # ==================================================================
    # resampling: bootstrap confidence bounds (bounds='npbb' / 'pbb')
    # and Monte-Carlo power analysis
    # ==================================================================
    @staticmethod
    def _cross_time(ev, curve, target):
        """Earliest event time at which the right-continuous step curve is
        at or below `target` (the standard quantile convention for a step
        survival estimator - it matches where the plotted step drops
        through the line). inf if the curve never gets that low."""
        curve = np.asarray(curve, dtype=float)
        if curve.size == 0 or curve[-1] > target:
            return np.inf
        j = int(np.searchsorted(-curve, -target))
        return float(ev[min(j, len(ev) - 1)])

    def _point_surv(self, xrow, times):
        """Fitted-model S(t | xrow) on a time grid (one covariate row)."""
        xr = np.asarray(xrow, dtype=float).ravel()
        if self.model == 'weibull_aft':
            lp = (self.intercept if self.fit_intercept else 0.0) + xr @ self.coef
            return np.exp(-np.exp(np.clip((np.log(times) - lp) / self.sigma,
                                          -self._R_CAP, self._R_CAP)))
        ev, bs = self.baseline_surv
        pos = np.searchsorted(ev, times, side='right') - 1
        S0 = np.where(pos >= 0, bs[np.clip(pos, 0, len(bs) - 1)], 1.0)
        return S0 ** float(np.exp((xr - self._x_mean) @ self.coef))

    def _point_bq(self, xrow, bq):
        """Fitted-model B_q life for one covariate row (time at which
        S = 1 - bq)."""
        xr = np.asarray(xrow, dtype=float).ravel()
        if self.model == 'weibull_aft':
            wq = float(np.log(-np.log(1.0 - bq)))
            return float(np.exp((self.intercept if self.fit_intercept else 0.0)
                                + xr @ self.coef + self.sigma * wq))
        ev, bs = self.baseline_surv
        S = bs ** float(np.exp((xr - self._x_mean) @ self.coef))
        return self._cross_time(ev, S, 1.0 - bq)

    def _sim_response(self, x_raw, censored, rng, truth):
        """Simulate (time, event) for the given covariate rows from a
        ground-truth model. `censored` is the pool of observed censoring
        times drawn from with replacement (empty -> no censoring)."""
        n = len(x_raw)
        if self.model == 'weibull_aft':
            lp = truth['intercept'] + x_raw @ truth['coef']
            w = np.clip(np.log(-np.log(rng.uniform(size=n))),
                        -self._R_CAP, self._R_CAP)
            t_ev = np.exp(lp + truth['sigma'] * w)
        else:
            theta = np.exp((x_raw - truth['x_mean']) @ truth['coef'])
            ev, H0 = truth['ev'], truth['cumhaz']
            need = rng.exponential(size=n) / theta
            k = np.searchsorted(H0, need, side='left')
            t_ev = np.where(k < len(ev), ev[np.clip(k, 0, len(ev) - 1)], np.inf)
        event = np.ones(n)
        if len(censored):
            c = rng.choice(np.asarray(censored, dtype=float), size=n,
                           replace=True)
            event = (t_ev <= c).astype(float)
            t = np.minimum(t_ev, c)
        else:
            t = t_ev.copy()
        bad = ~np.isfinite(t)
        if bad.any():
            fill = float(np.max(t[~bad])) if (~bad).any() else 1.0
            t = np.where(bad, fill, t)
            event = np.where(bad, 0.0, event)
        return np.clip(t, 1e-12, None), event

    def _refit_light(self, x_raw, t, d):
        """Refit the model on one resample; returns a compact bundle.
        Raises on a degenerate resample so the caller can skip it."""
        d = np.asarray(d, dtype=float)
        if np.any(x_raw.var(axis=0) == 0) or int(d.sum()) < 2:
            raise ValueError('degenerate resample')
        n = len(t)
        xs = (x_raw - self._x_mean) / self._x_std
        if self.model == 'weibull_aft':
            y = np.log(t)
            Zs = (np.column_stack([np.ones(n), xs]) if self.fit_intercept
                  else xs.copy())
            Zo = (np.column_stack([np.ones(n), x_raw]) if self.fit_intercept
                  else x_raw.copy())
            k = Zs.shape[1]
            ps, _, _ = self._aft_fit(Zs, y, d)
            b = ps[:k]
            if self.fit_intercept:
                b_o = np.empty(k)
                b_o[0] = b[0] - np.sum(b[1:] * self._x_mean / self._x_std)
                b_o[1:] = b[1:] / self._x_std
            else:
                b_o = b / self._x_std
            params, _, _ = self._aft_fit(
                Zo, y, d, start=np.concatenate([b_o, [ps[k]]]))
            if not np.all(np.isfinite(params)):
                raise ValueError('non-finite refit')
            p = len(self.coef)
            return {
                'kind': 'aft',
                'intercept': float(params[0]) if self.fit_intercept else 0.0,
                'coef': np.asarray(params[1:1 + p] if self.fit_intercept
                                   else params[:p], dtype=float),
                'sigma': float(np.exp(params[k])),
            }
        cs, _, _ = self._cox_fit(xs, t, d)
        coef, _, _ = self._cox_fit(x_raw, t, d, start=cs / self._x_std)
        xm = x_raw.mean(axis=0)
        ev, H0 = self._breslow_baseline((x_raw - xm) @ coef, t, d)
        return {'kind': 'cox', 'coef': np.asarray(coef, dtype=float),
                'x_mean': xm, 'ev': ev, 'bs': np.exp(-H0)}

    def _bootstrap(self, kind):
        """List of `n_boot` light refit bundles. kind='npbb' resamples
        whole units with replacement; kind='pbb' keeps the covariates and
        simulates the response from the fitted model. Cached per
        (kind, n_boot)."""
        key = (kind, self.n_boot)
        if key in self._boot_cache:
            return self._boot_cache[key]
        rng = np.random.default_rng(0)
        n = self.n
        censored = self._t[self._event == 0.0]
        truth = None
        if kind == 'pbb':
            if self.model == 'weibull_aft':
                truth = {
                    'intercept': self.intercept if self.fit_intercept else 0.0,
                    'coef': np.asarray(self.coef, dtype=float),
                    'sigma': float(self.sigma)}
            else:
                truth = {'coef': np.asarray(self.coef, dtype=float),
                         'x_mean': self._x_mean,
                         'ev': self.baseline_cumhaz[0],
                         'cumhaz': self.baseline_cumhaz[1]}
        out = []
        attempts = 0
        cap = self.n_boot * 3 + 20
        while len(out) < self.n_boot and attempts < cap:
            attempts += 1
            if kind == 'npbb':
                ix = rng.integers(0, n, size=n)
                xb, tb, db = self._x_raw[ix], self._t[ix], self._event[ix]
            else:
                xb = self._x_raw
                tb, db = self._sim_response(xb, censored, rng, truth)
            try:
                out.append(self._refit_light(xb, tb, db))
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
        need = min(self.n_boot, max(50, self.n_boot // 5))
        if len(out) < need:
            raise ValueError(
                f"bounds='{kind}': only {len(out)} of {self.n_boot} resamples "
                f"produced a usable fit (near-separation or too few events). "
                f"Use more data, fewer covariates, or bounds='fb'/'lrb'.")
        self._boot_cache[key] = out
        return out

    def _boot_coef_rows(self, kind):
        """(n_boot, len(coef_row)) array of bootstrap coefficient draws,
        columns in the order of self._coef_row (AFT: coefs then
        Intercept)."""
        rows = []
        for bd in self._bootstrap(kind):
            if bd['kind'] == 'aft':
                rows.append(np.append(bd['coef'], bd['intercept'])
                            if self.fit_intercept else bd['coef'])
            else:
                rows.append(bd['coef'])
        return np.asarray(rows, dtype=float)

    def _boot_surv(self, bd, xrow, times):
        xr = np.asarray(xrow, dtype=float).ravel()
        if bd['kind'] == 'aft':
            lp = bd['intercept'] + xr @ bd['coef']
            return np.exp(-np.exp(np.clip((np.log(times) - lp) / bd['sigma'],
                                          -self._R_CAP, self._R_CAP)))
        ev, bs = bd['ev'], bd['bs']
        pos = np.searchsorted(ev, times, side='right') - 1
        S0 = np.where(pos >= 0, bs[np.clip(pos, 0, len(bs) - 1)], 1.0)
        return S0 ** float(np.exp((xr - bd['x_mean']) @ bd['coef']))

    def _boot_band(self, xrow, times, cl, kind):
        curves = np.vstack([self._boot_surv(bd, xrow, times)
                            for bd in self._bootstrap(kind)])
        a = (1.0 - cl) / 2.0
        lower = np.quantile(curves, a, axis=0)
        upper = np.quantile(curves, 1.0 - a, axis=0)
        return self._point_surv(xrow, times), lower, upper

    def _boot_bq(self, xrow, bq, cl, kind):
        """(lower, point, upper) B_q life from the bootstrap draws."""
        xr = np.asarray(xrow, dtype=float).ravel()
        wq = float(np.log(-np.log(1.0 - bq)))
        target = 1.0 - bq
        vals = []
        for bd in self._bootstrap(kind):
            if bd['kind'] == 'aft':
                vals.append(float(np.exp(bd['intercept'] + xr @ bd['coef']
                                         + bd['sigma'] * wq)))
            else:
                S = bd['bs'] ** float(np.exp((xr - bd['x_mean']) @ bd['coef']))
                vals.append(self._cross_time(bd['ev'], S, target))
        q = np.array([v if np.isfinite(v) else np.nan for v in vals])
        a = (1.0 - cl) / 2.0
        if np.all(np.isnan(q)):
            return np.inf, np.inf, np.inf
        return (float(np.nanquantile(q, a)), self._point_bq(xr, bq),
                float(np.nanquantile(q, 1.0 - a)))

    # ---- drop-one likelihood-ratio tests (explicit data args) ------
    @staticmethod
    def _aft_drop1(Z, y, d, full, ll_full, idx):
        free = np.ones(full.size, dtype=bool)
        free[idx] = False

        def neg(z):
            v = full.copy()
            v[free] = z
            v[idx] = 0.0
            return Regression._aft_negll(v, Z, y, d)

        def jac(z):
            v = full.copy()
            v[free] = z
            v[idx] = 0.0
            return Regression._aft_neggrad(v, Z, y, d)[free]

        r = optimize.minimize(neg, full[free], jac=jac, method='L-BFGS-B',
                              options={'ftol': 1e-13, 'gtol': 1e-10})
        return max(0.0, 2.0 * (ll_full - (-r.fun)))

    def _cox_drop1(self, X, t, d, full, ll_full, idx):
        free = np.ones(full.size, dtype=bool)
        free[idx] = False
        b = full.copy()
        b[idx] = 0.0
        ll = self._cox_stats(b, X, t, d, self.ties)[0]
        for _ in range(self.max_iter):
            _, g, I = self._cox_stats(b, X, t, d, self.ties)
            gf = g[free]
            if np.max(np.abs(gf)) < self.tol:
                break
            If = I[np.ix_(free, free)]
            try:
                step = np.linalg.solve(If, gf)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(If, gf, rcond=None)[0]
            alpha = 1.0
            for _bt in range(30):
                cand = b.copy()
                cand[free] = b[free] + alpha * step
                ll_new = self._cox_stats(cand, X, t, d, self.ties)[0]
                if np.isfinite(ll_new) and ll_new >= ll - 1e-12:
                    break
                alpha *= 0.5
            else:
                break
            b, ll = cand, ll_new
        return max(0.0, 2.0 * (ll_full - ll))

    def _lrt_terms(self, x_raw, t, d):
        """(lr_per_covariate, lr_model) on one simulated data set: a
        drop-one LR statistic for every covariate, plus the overall model
        LR vs. the covariate-free null."""
        n, p = x_raw.shape
        d = np.asarray(d, dtype=float)
        if np.any(x_raw.var(axis=0) == 0) or int(d.sum()) < max(2, p + 1):
            raise ValueError('degenerate simulation')
        std = np.where(x_raw.std(axis=0, ddof=0) == 0, 1.0,
                       x_raw.std(axis=0, ddof=0))
        if self.model == 'weibull_aft':
            y = np.log(t)
            Z = (np.column_stack([np.ones(n), x_raw]) if self.fit_intercept
                 else x_raw.copy())
            full, _, ll_full = self._aft_fit(Z, y, d)
            if not np.all(np.isfinite(full)):
                raise ValueError('non-finite fit')
            c0 = 1 if self.fit_intercept else 0
            lr = np.array([self._aft_drop1(Z, y, d, full, ll_full, c0 + j)
                           for j in range(p)])
            Z0 = np.ones((n, 1)) if self.fit_intercept else np.zeros((n, 0))
            _, _, ll0 = self._aft_fit(Z0, y, d)
            return lr, max(0.0, 2.0 * (ll_full - ll0))
        xs = (x_raw - x_raw.mean(axis=0)) / std
        cs, _, _ = self._cox_fit(xs, t, d)
        full, _, ll_full = self._cox_fit(x_raw, t, d, start=cs / std)
        lr = np.array([self._cox_drop1(x_raw, t, d, full, ll_full, j)
                       for j in range(p)])
        ll0 = self._cox_stats(np.zeros(p), x_raw, t, d, self.ties)[0]
        return lr, max(0.0, 2.0 * (ll_full - ll0))

    # ---- public API ------------------------------------------------
    def power_analysis(self, n=None, n_sim=500, alpha=None, coef=None,
                       sigma=None, seed=0, print_report=True):
        """Monte-Carlo power for each covariate: the fraction of simulated
        data sets in which a drop-one likelihood-ratio test rejects
        H0: coef = 0 at level `alpha`.

        The fitted model is the ground truth. Pass `coef` (length p) and,
        for AFT, `sigma` to explore a different effect size / spread. The
        covariate design is drawn with replacement from the observed rows
        to size `n` (default: the fitted sample size); right-censoring
        times, if the data are censored, are resampled from the observed
        ones. Returns a pandas Series of powers indexed by covariate name,
        with an extra '(model)' entry for the overall LR test.
        """
        self._check_fitted()
        p = len(self.feature_names)
        n = int(self.n if n is None else n)
        if n < p + 2:
            raise ValueError(f'n must be at least {p + 2} for {p} covariates.')
        alpha = float((1.0 - self.cl) if alpha is None else alpha)
        truth_coef = np.asarray(self.coef if coef is None else coef,
                                dtype=float).ravel()
        if truth_coef.shape != (p,):
            raise ValueError(f'coef must have length {p}.')
        rng = np.random.default_rng(seed)
        censored = self._t[self._event == 0.0]
        if self.model == 'weibull_aft':
            truth = {
                'intercept': self.intercept if self.fit_intercept else 0.0,
                'coef': truth_coef,
                'sigma': float(self.sigma if sigma is None else sigma)}
        else:
            truth = {'coef': truth_coef, 'x_mean': self._x_mean,
                     'ev': self.baseline_cumhaz[0],
                     'cumhaz': self.baseline_cumhaz[1]}
        hits = np.zeros(p)
        model_hits = 0.0
        ok = 0
        for _ in range(int(n_sim)):
            ix = rng.integers(0, self.n, size=n)
            xb = self._x_raw[ix]
            tb, db = self._sim_response(xb, censored, rng, truth)
            try:
                lr, lr_model = self._lrt_terms(xb, tb, db)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
            ok += 1
            hits += (chi2.sf(lr, 1) < alpha).astype(float)
            model_hits += float(chi2.sf(lr_model, p) < alpha)
        if ok < max(10, int(n_sim) // 5):
            raise ValueError(
                f'power_analysis: only {ok} of {n_sim} simulations converged. '
                f'Increase n or simplify the model.')
        power = pd.Series(hits / ok, index=list(self.feature_names),
                          name='power')
        power['(model)'] = model_hits / ok
        if print_report:
            kind = ('Weibull AFT' if self.model == 'weibull_aft'
                    else f'Cox PH ({self.ties})')
            print(f'Monte-Carlo power  —  {kind}')
            print(f'n = {n}   simulations = {ok}/{int(n_sim)}   '
                  f'alpha = {alpha:g}   test = drop-one likelihood ratio')
            print(power.round(4).to_string())
        return power

    def sample_size(self, target_power=0.8, term=None, n_grid=None,
                    n_sim=300, alpha=None, seed=0, print_report=True):
        """Smallest sample size on a search grid at which the Monte-Carlo
        power reaches `target_power`. With `term=None` the criterion is the
        weakest covariate (min power over all terms); pass a covariate
        name to target just that one. Returns the winning n, or None if
        the grid tops out below the target. `n_grid` overrides the default
        geometric grid."""
        self._check_fitted()
        p = len(self.feature_names)
        if term is not None and term not in self.feature_names:
            raise ValueError(f'unknown term {term!r}.')
        if n_grid is None:
            lo = max(p + 2, 8)
            hi = max(int(4 * self.n), lo * 4)
            n_grid = np.unique(np.round(np.geomspace(lo, hi, 8)).astype(int))
        else:
            n_grid = np.unique(np.asarray(n_grid, dtype=int))
        curve = {}
        hit = None
        for nn in n_grid:
            pw = self.power_analysis(n=int(nn), n_sim=n_sim, alpha=alpha,
                                     seed=seed, print_report=False)
            val = (float(pw[term]) if term is not None
                   else float(pw.drop('(model)').min()))
            curve[int(nn)] = val
            if hit is None and val >= target_power:
                hit = int(nn)
        if print_report:
            lbl = term if term is not None else 'weakest covariate'
            print(f'Sample size for power >= {target_power:g}  ({lbl}, '
                  f'{self.model}, n_sim = {int(n_sim)})')
            for nn, v in curve.items():
                mark = '   <-- target reached' if hit == nn else ''
                print(f'  n = {nn:>6}   power = {v:.3f}{mark}')
            if hit is None:
                print(f'  target not reached up to n = {int(n_grid[-1])}')
        return hit

    # ==================================================================
    # goodness of fit
    # ==================================================================
    @staticmethod
    def _nelson_aalen(z, d):
        """Nelson-Aalen cumulative hazard of the sample (z, d), evaluated
        at every input point (ties share the summed 1/at-risk increment)."""
        z = np.asarray(z, dtype=float)
        d = np.asarray(d, dtype=float)
        order = np.argsort(z, kind='mergesort')
        ds = d[order]
        n = len(z)
        inc = np.where(ds > 0, 1.0 / (n - np.arange(n)), 0.0)
        H = np.empty(n)
        H[order] = np.cumsum(inc)
        return H

    @staticmethod
    def _kaplan_meier(z, d):
        """Pooled Kaplan-Meier estimate: (distinct times, S at those
        times), a right-continuous step function."""
        z = np.asarray(z, dtype=float)
        d = np.asarray(d, dtype=float)
        order = np.argsort(z, kind='mergesort')
        zs, ds = z[order], d[order]
        uniq = np.unique(zs)
        S = np.empty(len(uniq))
        s = 1.0
        for i, tau in enumerate(uniq):
            lo = np.searchsorted(zs, tau, side='left')
            hi = np.searchsorted(zs, tau, side='right')
            n_risk = len(zs) - lo
            n_ev = int(ds[lo:hi].sum())
            if n_ev and n_risk:
                s *= 1.0 - n_ev / n_risk
            S[i] = s
        return uniq, S

    def _cumhaz_at_own_time(self):
        """Fitted cumulative hazard H(t_i | x_i) for every unit - the raw
        material for the Cox-Snell / martingale / deviance residuals."""
        if self.model == 'weibull_aft':
            lp = ((self.intercept if self.fit_intercept else 0.0)
                  + self._x_raw @ self.coef)
            return np.exp(np.clip((np.log(self._t) - lp) / self.sigma,
                                  -self._R_CAP, self._R_CAP))
        ev, H0 = self.baseline_cumhaz
        pos = np.searchsorted(ev, self._t, side='right') - 1
        H0i = np.where(pos >= 0, H0[np.clip(pos, 0, len(H0) - 1)], 0.0)
        return H0i * np.exp((self._x_raw - self._x_mean) @ self.coef)

    def residuals(self, kind=None):
        """Model residuals, one value per unit (same order as the input).

        kind=None returns a DataFrame with all three columns; otherwise a
        single array for 'cox_snell', 'martingale' or 'deviance'.

        - 'cox_snell'  r_i = H(t_i | x_i). If the model is right, (r_i, d_i)
          is a right-censored sample from Exp(1): the Nelson-Aalen
          cumulative hazard of r plotted against r follows the 45deg line
          (see plot_gof / goodness_of_fit).
        - 'martingale' d_i - r_i, in (-inf, 1]. Plotted against a covariate
          its smoother shows whether that covariate enters with the right
          functional form.
        - 'deviance'   a symmetrised martingale residual, roughly N(0, 1)
          under the model; large |values| flag poorly-fit units.
        """
        self._check_fitted()
        d = self._event
        cs = self._cumhaz_at_own_time()
        mart = d - cs
        with np.errstate(divide='ignore', invalid='ignore'):
            dev = np.sign(mart) * np.sqrt(np.clip(
                -2.0 * (mart + d * np.log(np.where(d > 0, cs, 1.0))),
                0.0, np.inf))
        out = {'cox_snell': cs, 'martingale': mart, 'deviance': dev}
        if kind is None:
            return pd.DataFrame(out)
        if kind not in out:
            raise ValueError("kind must be None, 'cox_snell', 'martingale' "
                             "or 'deviance'.")
        return out[kind]

    # Thresholds that turn the raw diagnostics into a good / marginal /
    # poor verdict. The absolute-fit call is the primary one: it is the
    # departure of the Cox-Snell cumulative-hazard curve from the 45deg
    # line (slope away from 1, and the largest vertical gap) over the
    # lower 90% of the residuals - i.e. whether the assumed time model
    # (Weibull error for AFT, estimated baseline for Cox) reproduces the
    # observed event pattern. The concordance bands are the usual
    # rank-discrimination reading (0.5 = coin toss, 1 = perfect order).
    _GOF_SLOPE_GOOD, _GOF_SLOPE_POOR = 0.10, 0.25      # |slope - 1|
    _GOF_DEV_GOOD, _GOF_DEV_POOR = 0.15, 0.35          # max |NA(r) - r|
    _GOF_PH_POOR, _GOF_PH_MARGINAL = 0.01, 0.05        # global PH p-value

    def _gof_verdict(self, slope, max_dev, ph_p):
        """(rating dict, one-line overall verdict string) from the raw
        goodness-of-fit numbers. rating values are 'good' / 'marginal' /
        'poor' (absolute fit, ph) or a word band (discrimination)."""
        sd = abs(slope - 1.0) if np.isfinite(slope) else np.inf
        md = max_dev if np.isfinite(max_dev) else np.inf
        if sd <= self._GOF_SLOPE_GOOD and md <= self._GOF_DEV_GOOD:
            abs_fit = 'good'
        elif sd > self._GOF_SLOPE_POOR or md > self._GOF_DEV_POOR:
            abs_fit = 'poor'
        else:
            abs_fit = 'marginal'

        c = float(self.concordance)
        disc = ('weak' if c < 0.60 else 'modest' if c < 0.70
                else 'good' if c < 0.80 else 'strong')
        informative = self.lr_pvalue < 0.05

        ph = None
        if ph_p is not None and np.isfinite(ph_p):
            ph = ('poor' if ph_p < self._GOF_PH_POOR
                  else 'marginal' if ph_p < self._GOF_PH_MARGINAL else 'good')

        model_name = ('Weibull time model' if self.model == 'weibull_aft'
                      else 'Cox baseline hazard')
        overall = abs_fit
        if ph == 'poor' and overall == 'good':
            overall = 'marginal'
        if overall == 'poor':
            head = (f'POOR FIT - the assumed {model_name} is not consistent '
                    f'with the data')
        elif overall == 'marginal':
            head = (f'MARGINAL FIT - some departure from the assumed '
                    f'{model_name}; inspect plot_gof()')
        elif informative:
            head = (f'GOOD FIT - {model_name} consistent with the data, '
                    f'covariates informative')
        else:
            head = (f'GOOD FIT of the {model_name}, but the covariates add '
                    f'little (LR p = {self.lr_pvalue:.2g})')
        if ph == 'poor':
            head += f'; proportional-hazards assumption questionable (PH p = {ph_p:.2g})'
        elif ph == 'marginal':
            head += f'; watch the proportional-hazards assumption (PH p = {ph_p:.2g})'
        rating = {'absolute_fit': abs_fit, 'discrimination': disc}
        if ph is not None:
            rating['proportional_hazards'] = ph
        return rating, head

    def goodness_of_fit(self, decimals=4, print_report=True):
        """Overall fit summary as a pandas Series (and, by default, a
        printed report grouped into discrimination / relative fit /
        absolute fit, each metric tagged good / marginal / poor and closed
        by a one-line overall verdict).

        - concordance (Harrell's C)         - ranking / discrimination
          (0.5 chance ... 1.0 perfect; <0.6 weak, 0.6-0.7 modest,
          0.7-0.8 good, >0.8 strong).
        - log-likelihood, AIC               - relative fit / model choice.
        - LR test vs. the covariate-free null (p < 0.05 -> the covariates
          carry information).
        - cox_snell_slope                   - slope through the origin of
          the Nelson-Aalen cumulative hazard of the Cox-Snell residuals on
          the residuals themselves, over their lower 90%; ~1 for a good
          absolute fit (|slope - 1| <= 0.10 good, > 0.25 poor).
        - cox_snell_max_dev                 - largest |NA(r) - r| over the
          same range (<= 0.15 good, > 0.35 poor).
        - ph_pvalue (Cox only)              - global proportional-hazards
          test p-value; low means at least one effect drifts over time.

        The returned Series additionally carries the string ratings
        'absolute_fit', 'discrimination' (and 'proportional_hazards' for
        Cox) and the overall 'verdict'.
        """
        self._check_fitted()
        cs = self.residuals('cox_snell')
        Hna = self._nelson_aalen(cs, self._event)
        cut = float(np.quantile(cs, 0.90))
        m = (cs <= cut) & (cs > 0)
        if m.any():
            slope = float(np.sum(cs[m] * Hna[m]) / np.sum(cs[m] ** 2))
            max_dev = float(np.max(np.abs(Hna[m] - cs[m])))
        else:
            slope = max_dev = np.nan

        ph_p = None
        if self.model == 'cox_ph':
            try:
                ph_p = float(self.check_ph(print_report=False)
                             .loc['GLOBAL', 'p'])
            except Exception:
                ph_p = np.nan

        rating, verdict = self._gof_verdict(slope, max_dev, ph_p)

        stats = pd.Series({
            'n': float(self.n),
            'events': float(self.n_events),
            'concordance': float(self.concordance),
            'loglik': float(self.loglik),
            'aic': float(self.aic),
            'lr_stat': float(self.lr_stat),
            'lr_pvalue': float(self.lr_pvalue),
            'cox_snell_slope': slope,
            'cox_snell_max_dev': max_dev,
        }, dtype=object)
        if ph_p is not None:
            stats['ph_pvalue'] = ph_p
        stats['discrimination'] = rating['discrimination']
        stats['absolute_fit'] = rating['absolute_fit']
        if 'proportional_hazards' in rating:
            stats['proportional_hazards'] = rating['proportional_hazards']
        stats['verdict'] = verdict

        if print_report:
            def tag(r):
                return {'good': '[ good ]', 'marginal': '[ marg ]',
                        'poor': '[ POOR ]'}.get(r, f'[{r:^6s}]')
            name = ('Weibull AFT' if self.model == 'weibull_aft'
                    else f'Cox PH ({self.ties})')
            lp_flag = ('covariates matter' if self.lr_pvalue < 0.05
                       else 'covariates add little')
            print(f'Goodness of fit - {name}')
            print(f'n = {self.n}   events = {self.n_events}')
            print('discrimination')
            print(f'  concordance (Harrell C)     = {self.concordance:.{decimals}f}'
                  f'   {tag(rating["discrimination"])}  '
                  f'(0.5 chance .. 1.0 perfect)')
            print('relative fit')
            print(f'  AIC                         = {self.aic:.{decimals}f}')
            print(f'  LR test vs. null            = chi2 {self.lr_stat:.{decimals}f}'
                  f'   p = {self.lr_pvalue:.3g}   -> {lp_flag}')
            print(f'absolute fit (time model / link)   '
                  f'{tag(rating["absolute_fit"])}')
            print(f'  Cox-Snell NA slope          = {slope:.{decimals}f}'
                  f'   (want ~1.0)')
            print(f'  Cox-Snell max |NA(r) - r|   = {max_dev:.{decimals}f}'
                  f'   (want ~0.0)')
            if ph_p is not None and np.isfinite(ph_p):
                print(f'proportional hazards   {tag(rating["proportional_hazards"])}')
                print(f'  global Schoenfeld test      = p {ph_p:.3g}')
            print(f'overall: {verdict}')
        return stats

    gof = goodness_of_fit

    def _schoenfeld(self):
        """Raw Schoenfeld residuals for the Cox fit: (event_times, matrix
        of shape (#events, p)). Row k is x_i - xbar(tau_k) for the unit(s)
        failing at tau_k, xbar the risk-set average of exp(x'beta)-weighted
        covariates (Breslow handling of ties)."""
        X = self._x_raw
        t, dmask = self._t, self._event.astype(bool)
        eta = X @ self.coef
        w = np.exp(eta - eta.max())
        order = np.argsort(t, kind='mergesort')
        ts, Xs, ws, dm = t[order], X[order], w[order], dmask[order]
        S0 = 0.0
        S1 = np.zeros(X.shape[1])
        et, rows = [], []
        for tau in np.unique(ts)[::-1]:
            lo = np.searchsorted(ts, tau, side='left')
            hi = np.searchsorted(ts, tau, side='right')
            S0 += ws[lo:hi].sum()
            S1 += ws[lo:hi] @ Xs[lo:hi]
            em = dm[lo:hi]
            if em.any():
                xbar = S1 / S0
                for xi in Xs[lo:hi][em]:
                    et.append(tau)
                    rows.append(xi - xbar)
        if not rows:
            return np.zeros(0), np.zeros((0, X.shape[1]))
        return np.array(et[::-1]), np.array(rows[::-1])

    def check_ph(self, transform='km', decimals=4, print_report=True):
        """Cox only. Test the proportional-hazards assumption per covariate
        by looking for a trend of the scaled Schoenfeld residuals against a
        function of time (Grambsch & Therneau, 1994).

        transform : 'km' (1 - pooled Kaplan-Meier, the default), 'rank',
                    'log' or 'identity'.

        Returns a DataFrame indexed by covariate with columns 'test_stat'
        (chi-square, 1 df), 'p', 'rho_time' (correlation of the scaled
        residual with the time transform), plus a 'GLOBAL' row (chi-square
        with p df). A small p means the effect of that covariate changes
        over time - consider stratifying on it, adding a covariate-by-time
        interaction, or binning it.
        """
        self._check_fitted()
        if self.model != 'cox_ph':
            raise ValueError("check_ph is for model='cox_ph'.")
        et, S = self._schoenfeld()
        D, p = S.shape
        if D < max(4, p + 2):
            raise ValueError('too few events for a proportional-hazards test.')
        if transform == 'km':
            tk, Sk = self._kaplan_meier(self._t, self._event)
            pos = np.searchsorted(tk, et, side='right') - 1
            g = 1.0 - np.where(pos >= 0, Sk[np.clip(pos, 0, len(Sk) - 1)], 1.0)
        elif transform == 'rank':
            g = np.argsort(np.argsort(et)).astype(float)
        elif transform == 'log':
            g = np.log(et)
        elif transform == 'identity':
            g = et.astype(float)
        else:
            raise ValueError("transform must be 'km', 'rank', 'log' or "
                             "'identity'.")
        G = g - g.mean()
        SGG = float(G @ G)
        V = np.asarray(self.cov, dtype=float)
        u = S.T @ G
        Vu = V @ u
        vd = np.clip(np.diag(V), 1e-300, None)
        with np.errstate(divide='ignore', invalid='ignore'):
            T_p = D * Vu ** 2 / (SGG * vd)
        p_p = chi2.sf(T_p, 1)
        Sstar = self.coef[None, :] + D * (S @ V)
        rho = np.array([
            float(np.corrcoef(G, Sstar[:, j])[0, 1])
            if Sstar[:, j].std() > 0 else 0.0 for j in range(p)])
        Tg = float(D / SGG * (u @ V @ u)) if SGG > 0 else np.nan
        pg = float(chi2.sf(Tg, p)) if np.isfinite(Tg) else np.nan
        tab = pd.DataFrame({'test_stat': T_p, 'p': p_p, 'rho_time': rho},
                           index=list(self.feature_names))
        tab.loc['GLOBAL'] = [Tg, pg, np.nan]
        if print_report:
            print(f'Proportional-hazards test  (scaled Schoenfeld, g = '
                  f'{transform})')
            print(f'events = {D}   H0: each coefficient is constant over '
                  f'time')
            print(tab.round(decimals).to_string())
            per = tab.drop('GLOBAL')
            bad = list(per.index[per['p'] < 0.05])
            if bad:
                print(f'  -> PH questionable for: {", ".join(bad)}  '
                      f'(stratify, add a time interaction, or bin it).')
        return tab

    def _render_plot(self, fig, show=None, save=False):
        """Shared tail for the plotting methods. show=None falls back to
        the constructor's show=; the public plot* methods pass show=True
        by default, so a bare `r.plot_*()` in a notebook draws the figure
        and returns None - appearing exactly once, with no trailing ';'
        or plt.show() needed. Pass show=False to suppress the draw and get
        the Figure back for programmatic composition."""
        shown = self.show if show is None else bool(show)
        _finish_plot(fig, shown, bool(save) or self.save)
        return None if shown else fig

    def plot_gof(self, show=True):
        """Two diagnostic panels: (left) the Cox-Snell residuals - their
        Nelson-Aalen cumulative hazard against the residuals, which should
        track the dashed 45deg line; (right) martingale residuals against
        the linear predictor with a running-mean smoother, which should sit
        flat on zero if every covariate enters with the right shape."""
        self._check_fitted()
        _apply_plot_style(self.plot_style)
        cs = self.residuals('cox_snell')
        mart = self.residuals('martingale')
        ev = self._event.astype(bool)
        Hna = self._nelson_aalen(cs, self._event)

        # headline verdict, reusing goodness_of_fit's logic
        gof = self.goodness_of_fit(print_report=False)
        verdict = str(gof['verdict'])

        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))
        ax1 = fig.add_subplot(1, 2, 1)
        o = np.argsort(cs)
        hi = float(np.quantile(cs, 0.98)) or 1.0
        ax1.step(cs[o], Hna[o], where='post', color=PREDICTR_FIT_COLOR,
                 linewidth=1.5, zorder=2)
        ax1.plot(cs[ev], Hna[ev], 'o', ms=3, color=PREDICTR_FIT_COLOR, alpha=0.5,
                 zorder=3)
        ax1.plot([0.0, hi], [0.0, hi], linestyle='--', color='#8b8b8b',
                 linewidth=1, zorder=1)
        ax1.set_xlim(0.0, hi)
        ax1.set_ylim(0.0, hi)
        ax1.set_xlabel('Cox-Snell residual  r', color='black',
                       fontsize=self.xy_fontsize)
        ax1.set_ylabel('Nelson-Aalen cum. hazard of r', color='black',
                       fontsize=self.xy_fontsize)
        ax1.set_title('Cox-Snell (follow the 45° line)',
                      fontsize=self.plot_title_fontsize)
        ax1.tick_params(labelsize=self.tick_fontsize)
        ax1.grid(True, which='major')

        ax2 = fig.add_subplot(1, 2, 2)
        lp = (self._x_raw - self._x_mean) @ self.coef
        oo = np.argsort(lp)
        ax2.axhline(0.0, color='#8b8b8b', linestyle='--', linewidth=1,
                    zorder=1)
        ax2.plot(lp[oo], mart[oo], 'o', ms=3, color=PREDICTR_FIT_COLOR, alpha=0.5,
                 zorder=2)
        k = max(5, len(lp) // 8)
        if len(lp) > k:
            csum = np.cumsum(np.insert(mart[oo], 0, 0.0))
            sm = (csum[k:] - csum[:-k]) / k
            ax2.plot(lp[oo][k // 2: k // 2 + len(sm)], sm, color='crimson',
                     linewidth=1.6, zorder=3)
        ax2.set_xlabel(r'linear predictor  $(x-\bar{x})^\top\beta$',
                       color='black', fontsize=self.xy_fontsize)
        ax2.set_ylabel('martingale residual', color='black',
                       fontsize=self.xy_fontsize)
        ax2.set_title('Functional form (smoother ~ flat at 0)',
                      fontsize=self.plot_title_fontsize)
        ax2.tick_params(labelsize=self.tick_fontsize)
        ax2.grid(True, which='major')

        colour = ('#1a7f37' if verdict.startswith('GOOD') else
                  '#b54708' if verdict.startswith('MARGINAL') else '#b42318')
        fig.suptitle(verdict, fontsize=max(9, self.plot_title_fontsize - 2),
                     color=colour)
        plt.tight_layout(rect=(0, 0, 1, 0.94))
        if self.save:
            try:
                plt.savefig(self.save_path)
            except Exception:
                raise ValueError('Path is faulty.')
        return self._render_plot(fig, show)

    # ==================================================================
    # descriptive (non-parametric, covariate-free) estimators
    #
    # Kaplan-Meier survival S(t) and Nelson-Aalen cumulative hazard H(t)
    # for the data as supplied - the model-free description the fitted
    # AFT / Cox curves are compared against. No fit() required.
    #
    # At each distinct observed time tau with d events out of n at risk:
    #   S(tau) = prod (1 - d/n)         (Kaplan-Meier product limit)
    #   H(tau) = sum  d/n               (Nelson-Aalen)
    # Pointwise standard errors:
    #   Var(S)/S^2 = sum d / (n (n - d))            (Greenwood)
    #   Var(H)     = sum d / n^2
    # The KM interval is formed on the scale ln(-ln S) (so it stays in
    # (0, 1)); the NA interval on ln H (so it stays positive).
    # ==================================================================
    def _group_indices(self, by):
        """[(label, row indices)] for by=None (one group, label None), a
        raw covariate column name, or a length-n array of group labels."""
        if by is None:
            return [(None, np.arange(self.n))]
        if isinstance(by, str):
            if by not in self._x_frame.columns:
                raise ValueError(
                    f"by='{by}' is not a covariate column; choose from "
                    f"{list(self._x_frame.columns)} or pass a length-n "
                    f"array of group labels.")
            lab = self._x_frame[by].to_numpy()
            pre = by
        else:
            lab = np.asarray(by)
            if lab.shape[0] != self.n:
                raise ValueError(f'by must have length n = {self.n}.')
            pre = 'group'
        return [(f'{pre}={v}', np.where(lab == v)[0]) for v in pd.unique(lab)]

    def _life_table(self, idx, cl):
        """KM + NA life table for the units in idx (shared engine with
        Analysis): one row per distinct observed time, pointwise limits
        at level cl."""
        idx = np.asarray(idx)
        return _km_na_life_table(self._t[idx], self._event[idx], cl)

    def _descriptive(self, by, cl, cols):
        cl = float(self.cl if cl is None else cl)
        if not 0.0 < cl < 1.0:
            raise ValueError('cl must be between 0 and 1.')
        frames = []
        for label, idx in self._group_indices(by):
            tab = self._life_table(idx, cl)[cols].copy()
            if label is not None:
                tab.insert(0, 'group', label)
            frames.append(tab)
        return pd.concat(frames, ignore_index=True)

    def kaplan_meier(self, by=None, cl=None):
        """Kaplan-Meier survival estimate S(t) for the data as supplied -
        non-parametric, covariate-free (no fit() needed). This is the
        model-free reference the fitted survival curves are compared
        against.

        by : str or array-like, optional
            Split into groups: a raw covariate column name (e.g.
            'material'), or a length-n array of group labels. One curve
            per level.
        cl : float, optional
            Confidence level for the pointwise band (default: self.cl).
            Greenwood standard error; interval on the ln(-ln S) scale.

        Returns a DataFrame: [group,] time, n_risk, n_event, n_censor,
        surv, surv_se, surv_lower, surv_upper. Right-continuous step
        function - surv is the value from `time` up to the next row.
        """
        return self._descriptive(
            by, cl, ['time', 'n_risk', 'n_event', 'n_censor',
                     'surv', 'surv_se', 'surv_lower', 'surv_upper'])

    def nelson_aalen(self, by=None, cl=None):
        """Nelson-Aalen cumulative-hazard estimate H(t) for the data as
        supplied - non-parametric, covariate-free (no fit() needed).

        by, cl : as for kaplan_meier(); the NA interval is on the ln H
        scale with Var(H) = sum d / n^2.

        Returns a DataFrame: [group,] time, n_risk, n_event, n_censor,
        cumhaz, cumhaz_se, cumhaz_lower, cumhaz_upper.
        """
        return self._descriptive(
            by, cl, ['time', 'n_risk', 'n_event', 'n_censor',
                     'cumhaz', 'cumhaz_se', 'cumhaz_lower', 'cumhaz_upper'])

    def _plot_descriptive(self, kind, by, ci, cl, show):
        cl = float(self.cl if cl is None else cl)
        if not 0.0 < cl < 1.0:
            raise ValueError('cl must be between 0 and 1.')
        tables = [(label, self._life_table(idx, cl))
                  for label, idx in self._group_indices(by)]
        fig = _plot_km_na(self, tables, kind, ci)
        return self._render_plot(fig, show)

    def plot_km(self, by=None, ci=True, cl=None, show=True):
        """Step plot of the Kaplan-Meier survival estimate (one line per
        group with by=), pointwise band when ci=True, censoring marks as
        vertical ticks. show=True (default) draws it and returns None;
        show=False returns the Figure instead."""
        return self._plot_descriptive('km', by, ci, cl, show)

    def plot_na(self, by=None, ci=True, cl=None, show=True):
        """Step plot of the Nelson-Aalen cumulative-hazard estimate (one
        line per group with by=), pointwise band when ci=True. show as for
        plot_km()."""
        return self._plot_descriptive('na', by, ci, cl, show)

    # ==================================================================
    # named life-stress (aging) laws: results, prediction, diagnostics
    # ==================================================================
    def _stress_term_index(self):
        """Column indices of every life-stress design term, in
        self.feature_names / self.coef order."""
        return [self.feature_names.index(c)
                for raw in self._stress_spec
                for c in self._stress_terms[raw]]

    def _build_stress_params(self):
        """Fill self.stress_params: one row per physical parameter derived
        from the fitted life-stress coefficients (honours bounds= via
        ci_lower/ci_upper). Weibull AFT only - on the Cox hazard scale the
        coefficients are not life-stress-law parameters."""
        if not self._stress_spec or self.model != 'weibull_aft':
            self.stress_params = None
            return
        fn = self.feature_names
        lo = self.ci_lower if self.ci_lower is not None else np.full(
            len(self.coef), np.nan)
        hi = self.ci_upper if self.ci_upper is not None else np.full(
            len(self.coef), np.nan)
        rows = []
        for raw, law in self._stress_spec.items():
            j = [fn.index(c) for c in self._stress_terms[raw]]
            c, s, l, u = (self.coef[j], self.se_coef[j], lo[j], hi[j])
            if law == 'arrhenius':
                rows.append((raw, law, 'Ea_eV', c[0], s[0], l[0], u[0]))
            elif law == 'eyring':
                rows.append((raw, law, 'Ea_eV', c[0], s[0], l[0], u[0]))
                rows.append((raw, law, 'lnT_coef', c[1], s[1], l[1], u[1]))
            elif law == 'inverse_power':
                rows.append((raw, law, 'n', -c[0], s[0], -u[0], -l[0]))
            elif law == 'coffin_manson':
                rows.append((raw, law, 'fatigue_exponent',
                             -c[0], s[0], -u[0], -l[0]))
            else:  # exponential
                rows.append((raw, law, 'rate_coef', c[0], s[0], l[0], u[0]))
        self.stress_params = pd.DataFrame(
            rows, columns=['stress', 'law', 'parameter', 'value', 'se',
                           'ci_lower', 'ci_upper'])

    def acceleration_factor(self, from_, to_, cl=None):
        """Weibull AFT + stress_model: the lifetime ratio
        life(to_) / life(from_) for two operating points given in raw
        stress units (dicts, e.g. {'temp': 55, 'volt': 3.3}). > 1 means
        `to_` lasts longer. Returns (factor, lower, upper) at level cl
        (delta method on the log ratio)."""
        self._check_fitted()
        if self.model != 'weibull_aft':
            raise ValueError("acceleration_factor is for model='weibull_aft'.")
        if not self._stress_spec:
            raise ValueError('acceleration_factor needs a stress_model.')
        a = self._prep_X(from_)[0]
        b = self._prep_X(to_)[0]
        dvec = b - a
        d = float(dvec @ self.coef)
        var = float(dvec @ self._coef_cov @ dvec)
        cl = self.cl if cl is None else float(cl)
        z = norm.ppf(1.0 - (1.0 - cl) / 2.0)
        se = np.sqrt(max(var, 0.0))
        return (float(np.exp(d)), float(np.exp(d - z * se)),
                float(np.exp(d + z * se)))

    def _marginal_weibull(self, idx):
        """(location, sigma) of an intercept-only Weibull AFT fitted to the
        units in `idx`. None if fewer than one failure there."""
        idx = np.asarray(idx)
        d = self._event[idx]
        if d.sum() < 1:
            return None
        pr, _, _ = self._aft_fit(np.ones((len(idx), 1)),
                                 np.log(self._t[idx]), d)
        return float(pr[0]), float(np.exp(pr[1]))

    def _stress_groups(self):
        """{rounded stress-term tuple: row indices} over the training data."""
        sidx = self._stress_term_index()
        key = np.round(self._x_raw[:, sidx], 9)
        groups = {}
        for i in range(self.n):
            groups.setdefault(tuple(key[i]), []).append(i)
        return {k: np.asarray(v) for k, v in groups.items()}

    def check_shape(self, decimals=4, print_report=True):
        """Weibull AFT + stress_model: fit a separate Weibull to each
        tested stress level and compare the shape estimates - the scale-
        accelerated (SAFT) model assumes one common shape. Returns a
        DataFrame (n, events, beta per level) with a 'pooled' row."""
        self._check_fitted()
        if self.model != 'weibull_aft' or not self._stress_spec:
            raise ValueError("check_shape needs model='weibull_aft' with a "
                             "stress_model.")
        rows = []
        for _, ii in self._stress_groups().items():
            d = self._event[ii]
            beta = np.nan
            if d.sum() >= 2:
                mw = self._marginal_weibull(ii)
                beta = 1.0 / mw[1] if mw else np.nan
            rows.append((len(ii), int(d.sum()), beta))
        tab = pd.DataFrame(rows, columns=['n', 'events', 'beta'])
        tab.loc['pooled'] = [self.n, self.n_events, self.beta]
        fin = tab['beta'].iloc[:-1].dropna()
        ratio = float(fin.max() / fin.min()) if len(fin) >= 2 else np.nan
        consistent = not (ratio > 2.0)
        if print_report:
            print('Weibull shape (beta) per tested stress level  '
                  '(SAFT assumes one common shape)')
            print(tab.round(decimals).to_string())
            if len(fin) >= 2:
                verdict = ('consistent' if consistent else
                           'SHAPES DIFFER - the constant-shape assumption '
                           'looks violated')
                print(f'max/min = {ratio:.2f}  ->  {verdict}')
        return tab

    def plot_stress_life(self, q=0.5, show=True):
        """Weibull AFT + stress_model: the life-stress diagnostic. For each
        tested stress level the marginal B(100q) life is plotted (log y)
        against that level's stress linear predictor, with the fitted AFT
        line. Roughly collinear points support the aging law."""
        self._check_fitted()
        if self.model != 'weibull_aft' or not self._stress_spec:
            raise ValueError("plot_stress_life needs model='weibull_aft' "
                             "with a stress_model.")
        sidx = self._stress_term_index()
        pidx = [j for j in range(len(self.feature_names)) if j not in sidx]
        wq = float(np.log(-np.log(1.0 - q)))
        xs, ys = [], []
        for _, ii in self._stress_groups().items():
            mw = self._marginal_weibull(ii)
            if mw is None:
                continue
            loc, sig = mw
            xs.append(float(self._x_raw[ii[0], sidx] @ self.coef[sidx]))
            ys.append(np.exp(loc + sig * wq))
        xs, ys = np.asarray(xs), np.asarray(ys)
        order = np.argsort(xs)

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))
        ax = plt.gca()
        ax.plot(xs[order], ys[order], 'o', color=PREDICTR_FIT_COLOR, markersize=6,
                label=f'marginal B{q * 100:g} per stress level')
        plain_mean = (float(self._x_raw[:, pidx].mean(axis=0)
                            @ self.coef[pidx]) if pidx else 0.0)
        gx = np.linspace(xs.min(), xs.max(), 100)
        gy = np.exp(self.intercept + gx + plain_mean + self.sigma * wq)
        ax.plot(gx, gy, '-', color='crimson', linewidth=1.6,
                label='AFT life-stress line')
        ax.set_yscale('log')
        ax.set_xlabel(r'stress linear predictor  $\sum_j g_j(S)\,\theta_j$',
                      color='black', fontsize=self.xy_fontsize)
        ax.set_ylabel(f'B{q * 100:g} life', color='black',
                      fontsize=self.xy_fontsize)
        ax.set_title('Life-stress relationship', color='black',
                     fontsize=self.plot_title_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(True, which='both')
        if self.show_legend:
            ax.legend(fontsize=self.legend_fontsize)
        plt.tight_layout()
        if self.save:
            try:
                plt.savefig(self.save_path)
            except Exception:
                raise ValueError('Path is faulty.')
        return self._render_plot(fig, show)

    # ------------------------------------------------------------------
    # summary
    # ------------------------------------------------------------------
    def _build_summary(self):
        cl_lbl = f'{self.cl:.2f}'
        # exp(coef) can legitimately overflow to inf when a covariate is on
        # a tiny scale (e.g. 1/T in Kelvin): a one-unit change is then
        # enormous. inf is the intended cell value, not an error.
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            expc = np.exp(self._coef_row)
            exp_lo = np.exp(self._lo_row)
            exp_hi = np.exp(self._hi_row)
            z = self._coef_row / self._se_row
        pval = 2.0 * norm.sf(np.abs(z))
        exp_lbl = ('exp(coef) [time ratio]' if self.model == 'weibull_aft'
                   else 'exp(coef) [hazard ratio]')
        data = {
            'coef': self._coef_row,
            exp_lbl: expc,
            'se(coef)': self._se_row,
            'z': z,
            'p': pval,
            f'coef lower {cl_lbl}': self._lo_row,
            f'coef upper {cl_lbl}': self._hi_row,
            f'exp(coef) lower {cl_lbl}': exp_lo,
            f'exp(coef) upper {cl_lbl}': exp_hi,
        }
        self.summary_ = pd.DataFrame(data, index=self._row_names)

    def summary(self, decimals=4, print_report=True):
        """Return the coefficient table as a DataFrame. With
        print_report=True (default) also print the full report - header
        block (model, sample sizes, log-likelihood, AIC, likelihood-ratio
        test, plus sigma/shape for AFT or concordance for Cox) and the
        table."""
        if self.summary_ is None:
            raise ValueError('Call fit() before summary().')
        table = self.summary_.round(decimals)
        if print_report:
            name = ('Weibull AFT' if self.model == 'weibull_aft'
                    else f'Cox PH (ties: {self.ties})')
            lines = [
                f'Lifetime regression — {name}',
                f'n = {self.n}   events = {self.n_events}   '
                f'censored = {self.n - self.n_events}',
                f'covariates: {", ".join(self.feature_names)}',
            ]
            if self.bounds is not None:
                bname = {'fb': 'Wald', 'lrb': 'profile-likelihood',
                         'npbb': 'non-parametric bootstrap',
                         'pbb': 'parametric bootstrap'}[self.bounds]
                extra = (f', n_boot = {self.n_boot}'
                         if self.bounds in ('npbb', 'pbb') else '')
                lines.append(f'bounds: {bname} ({self.bounds}), '
                             f'{self.bounds_type} @ {self.cl * 100:g}%{extra}')
            lines.append(f'log-likelihood = {self.loglik:.{decimals}f}   '
                         f'AIC = {self.aic:.{decimals}f}')
            lines.append(f'LR test vs. null:  chi2 = {self.lr_stat:.{decimals}f}   '
                         f'df = {len(self.feature_names)}   '
                         f'p = {self.lr_pvalue:.3g}')
            if self.model == 'weibull_aft':
                lines.append(f'sigma = {self.sigma:.{decimals}f} '
                             f'(se {self.se_sigma:.{decimals}f})   ->   '
                             f'Weibull shape = {self.beta:.{decimals}f}')
            lines.append(f'concordance = {self.concordance:.{decimals}f}')
            if self.stress_params is not None:
                has_ci = self.bounds is not None
                head = (f'{self.cl * 100:g}% CI' if has_ci
                        else 'no bounds requested')
                lines.append(f'life-stress model (physical parameters, '
                             f'{head}):')
                for _, r in self.stress_params.iterrows():
                    tail = (f'  ({r["ci_lower"]:.{decimals}f}, '
                            f'{r["ci_upper"]:.{decimals}f})' if has_ci else '')
                    lines.append(
                        f'  {r["stress"]:<12s} {r["law"]:<14s} '
                        f'{r["parameter"]} = {r["value"]:.{decimals}f}{tail}')
            print('\n'.join(lines))
            print()
            print(table.to_string())
        return table

    # ------------------------------------------------------------------
    # prediction
    # ------------------------------------------------------------------
    def _check_fitted(self):
        if self.coef is None:
            raise ValueError('Call fit() before predicting.')

    def _lin_pred(self, X, centered):
        """Linear predictor for new covariate rows. centered=True measures
        each row against the training-set mean covariate vector (used for
        the ratio predictions and for Cox survival, whose baseline is at
        the mean); centered=False adds the AFT intercept and is the
        absolute location on the log-time scale."""
        Xn = self._prep_X(X)
        if centered:
            return (Xn - self._x_mean) @ self.coef
        eta = Xn @ self.coef
        if self.model == 'weibull_aft':
            eta = eta + self.intercept
        return eta

    def predict_time_ratio(self, X):
        """AFT: exp((x - xbar)'coef) - the multiplicative effect on the
        predicted lifetime relative to a unit at the mean covariate
        values."""
        self._check_fitted()
        if self.model != 'weibull_aft':
            raise ValueError("predict_time_ratio is for model='weibull_aft'.")
        return np.exp(self._lin_pred(X, centered=True))

    def predict_hazard_ratio(self, X):
        """Cox: exp((x - xbar)'beta) - the hazard relative to a unit at the
        mean covariate values."""
        self._check_fitted()
        if self.model != 'cox_ph':
            raise ValueError("predict_hazard_ratio is for model='cox_ph'.")
        return np.exp(self._lin_pred(X, centered=True))

    def predict_quantile(self, X, q=0.5, ci=False, cl=None, bounds=None):
        """AFT: the time t with F(t | x) = q, per covariate row of X. With
        ci=True also return the (lower, upper) confidence limits - delta
        method for bounds='fb', profile likelihood for bounds='lrb'
        (defaults to the model's own bounds= setting)."""
        self._check_fitted()
        if self.model != 'weibull_aft':
            raise ValueError("predict_quantile is for model='weibull_aft'.")
        Xn = self._prep_X(X)
        wq = float(np.log(-np.log(1.0 - q)))
        lp = (self.intercept if self.fit_intercept else 0.0) + Xn @ self.coef
        tq = np.exp(lp + self.sigma * wq)
        if not ci:
            return tq
        cl, method = self._resolve_band_opts(cl, bounds)
        z = norm.ppf(1.0 - (1.0 - cl) / 2.0)
        lo = np.empty(len(Xn))
        hi = np.empty(len(Xn))
        for i, xr in enumerate(Xn):
            if method in ('npbb', 'pbb'):
                lo[i], _, hi[i] = self._boot_bq(xr, q, cl, method)
            elif method == 'lrb' and self.fit_intercept:
                lo[i], hi[i] = self._aft_quantile_profile_ci(xr, wq, cl)
            else:
                g = np.concatenate(([1.0] if self.fit_intercept else [],
                                    xr, [self.sigma * wq]))
                se = float(np.sqrt(max(g @ self.cov @ g, 0.0)))
                lo[i] = np.exp(np.log(tq[i]) - z * se)
                hi[i] = np.exp(np.log(tq[i]) + z * se)
        return tq, lo, hi

    def _aft_quantile_profile_ci(self, xrow, wq, cl):
        xrow = np.asarray(xrow, dtype=float).ravel()
        y, d, Z = self._y, self._event, self._Z
        p = len(self.coef)
        thr = chi2.ppf(cl, 1)
        ll_max = -self._aft_negll(self._aft_params, Z, y, d)
        free0 = np.concatenate([self._aft_params[1:1 + p],
                                [self._aft_params[-1]]])
        lnt_hat = self.intercept + xrow @ self.coef + self.sigma * wq

        def restr(free, lnt0):
            coef_f, logsig_f = free[:p], free[p]
            sig_f = np.exp(logsig_f)
            b0 = lnt0 - xrow @ coef_f - sig_f * wq
            full = np.concatenate([[b0], coef_f, [logsig_f]])
            gf = self._aft_neggrad(full, Z, y, d)
            return (self._aft_negll(full, Z, y, d),
                    np.concatenate([gf[1:1 + p] - gf[0] * xrow,
                                    [gf[p + 1] - gf[0] * sig_f * wq]]))

        def prof(lnt0):
            r = optimize.minimize(restr, free0, args=(lnt0,), jac=True,
                                  method='L-BFGS-B',
                                  options={'ftol': 1e-13, 'gtol': 1e-9})
            return -r.fun

        def g(lnt0):
            return 2.0 * (ll_max - prof(lnt0)) - thr

        span = max(abs(self.sigma), 0.5)
        out = []
        for sgn in (-1.0, 1.0):
            step, v, it = span, lnt_hat + sgn * span, 0
            while g(v) < 0 and it < 40:
                step *= 1.6
                v += sgn * step
                it += 1
            try:
                out.append(float(np.exp(optimize.brentq(
                    g, min(lnt_hat, v), max(lnt_hat, v), xtol=1e-6,
                    maxiter=100))))
            except (ValueError, RuntimeError):
                out.append(float(np.exp(lnt_hat + sgn * 6.0 * span)))
        return out[0], out[1]

    def _cox_quantile_ci(self, xrow, q, cl):
        """B_q life for a Cox profile: (lower, point, upper) - the earliest
        event time at which S(t | x) (resp. the band edges) is at or below
        1 - q. Step convention, so the point value lands exactly where the
        plotted step survival curve drops through the line. inf if a curve
        never reaches 1 - q within the data."""
        pc = self._cox_band_pieces()
        ev = pc['ev']
        S, lo, hi, _, _ = self._cox_delta_band(xrow, ev, cl)
        target = 1.0 - q
        # lo <= S <= hi, so the lower edge reaches the line first: the
        # triplet comes out ordered (u <= m <= o).
        return (self._cross_time(ev, lo, target),
                self._cross_time(ev, S, target),
                self._cross_time(ev, hi, target))

    def _bq_life(self, xrow, bq, cl, method):
        if method in ('npbb', 'pbb'):
            return self._boot_bq(np.asarray(xrow), bq, cl, method)
        if self.model == 'weibull_aft':
            tq, lo, hi = self.predict_quantile(np.asarray(xrow)[None, :], q=bq,
                                               ci=True, cl=cl, bounds=method)
            return float(lo[0]), float(tq[0]), float(hi[0])
        return self._cox_quantile_ci(np.asarray(xrow), bq, cl)

    def predict_median(self, X):
        """AFT: the median predicted lifetime for each covariate row."""
        return self.predict_quantile(X, 0.5)

    def predict_survival(self, X, times=None, ci=False, cl=None, bounds=None,
                         simultaneous=False):
        """Predicted survival S(t | x).

        ci=False (default): a DataFrame indexed by time, one column per
        covariate row of X. With times=None a dense grid from just above 0
        to slightly past the largest observed time is used.

        ci=True: a dict with keys 'surv', 'lower', 'upper' (pointwise
        confidence band - delta method for bounds='fb' / profile
        likelihood for AFT bounds='lrb', both on the complementary-log-log
        scale, or pointwise resample percentiles for bounds='npbb'/'pbb'),
        'lower_sim'/'upper_sim' (simultaneous band over the observed time
        range, None unless simultaneous=True; not available for the
        bootstrap methods), 'in_data_range' (bool array over the time grid),
        'method', 'cl'. cl overrides the confidence level."""
        self._check_fitted()
        if times is None:
            hi_t = float(self._t.max()) * 1.02
            times = np.linspace(hi_t * 1e-4, hi_t, 200)
        times = np.atleast_1d(np.asarray(times, dtype=float))
        Xn = self._prep_X(X)
        if self.model == 'weibull_aft':
            lp = (self.intercept if self.fit_intercept else 0.0) + Xn @ self.coef
            S = np.exp(-np.exp((np.log(times)[:, None] - lp[None, :]) / self.sigma))
        else:
            lpc = (Xn - self._x_mean) @ self.coef
            bt, bs = self.baseline_surv
            pos = np.searchsorted(bt, times, side='right') - 1
            S0 = np.where(pos >= 0, bs[np.clip(pos, 0, len(bs) - 1)], 1.0)
            S = S0[:, None] ** np.exp(lpc)[None, :]
        idx = pd.Index(times, name='time')
        cols = (list(X.index) if isinstance(X, pd.DataFrame)
                else list(range(len(Xn))))
        if not ci:
            return pd.DataFrame(S, index=idx, columns=cols)

        cl, method = self._resolve_band_opts(cl, bounds)
        low = np.empty_like(S)
        up = np.empty_like(S)
        slo = np.full_like(S, np.nan)
        shi = np.full_like(S, np.nan)
        any_sim = False
        for j, xr in enumerate(Xn):
            s_j, lo_j, hi_j, sl_j, sh_j = self._survival_band(
                xr, times, cl, method, simultaneous)
            S[:, j], low[:, j], up[:, j] = s_j, lo_j, hi_j
            if sl_j is not None:
                any_sim = True
                slo[:, j], shi[:, j] = sl_j, sh_j
        return {
            'time': times,
            'surv': pd.DataFrame(S, index=idx, columns=cols),
            'lower': pd.DataFrame(low, index=idx, columns=cols),
            'upper': pd.DataFrame(up, index=idx, columns=cols),
            'lower_sim': pd.DataFrame(slo, index=idx, columns=cols) if any_sim else None,
            'upper_sim': pd.DataFrame(shi, index=idx, columns=cols) if any_sim else None,
            'in_data_range': (times >= self._t.min()) & (times <= self._t.max()),
            'method': method, 'cl': cl, 'simultaneous': bool(simultaneous),
        }

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    def plot(self, show=True):
        """Forest plot of the coefficients with their confidence bounds
        (log-time scale for AFT, log-hazard scale for Cox). A dashed line
        marks no effect (coef = 0)."""
        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))
        ax = plt.gca()

        names = list(self.feature_names)
        vals = np.asarray(self.coef, dtype=float)
        lo = np.asarray(self.ci_lower, dtype=float)
        hi = np.asarray(self.ci_upper, dtype=float)
        ypos = np.arange(len(names))[::-1]

        left = np.where(np.isfinite(lo), vals - lo, 0.0)
        right = np.where(np.isfinite(hi), hi - vals, 0.0)
        ax.errorbar(vals, ypos, xerr=np.vstack([left, right]), fmt='o',
                    color=PREDICTR_FIT_COLOR, ecolor='royalblue', elinewidth=1.5,
                    capsize=3, markersize=5, zorder=3)
        ax.axvline(0.0, color='#8b8b8b', linestyle='--', linewidth=1, zorder=1)
        ax.set_yticks(ypos)
        ax.set_yticklabels(names, fontsize=self.tick_fontsize)
        ax.set_xlabel('coefficient (log-time scale)' if self.model == 'weibull_aft'
                      else 'coefficient (log-hazard scale)',
                      color='black', fontsize=self.xy_fontsize)
        ax.set_title(self.plot_title, color='black',
                     fontsize=self.plot_title_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(True, which='major', axis='x')
        plt.tight_layout()
        if self.save:
            try:
                plt.savefig(self.save_path)
            except Exception:
                raise ValueError('Path is faulty.')
        return self._render_plot(fig, show)

    def _draw_extrapolation_marker(self, ax, t_hi, t_max, **over):
        """Dotted vertical at the last observed time t_hi, tagged with a
        small 't_max' label reading bottom-to-top in the gap near the top
        of the axes. Everything to the right of the line is extrapolation."""
        st = dict(self._EXTRAP_STYLE, **over)
        x0, x1 = ax.get_xlim()
        if not (x0 <= t_hi <= x1):
            return
        ax.axvline(t_hi, color=st['line_color'], linestyle=st['line_style'],
                   linewidth=st['line_width'], zorder=1)
        if not st.get('label'):
            return
        fs = st['label_fs'] or max(6.0, self.tick_fontsize - 3)
        ax.text(t_hi, float(st['label_y']), st['label'],
                transform=ax.get_xaxis_transform(), rotation=90,
                rotation_mode='anchor', ha='center', va='center',
                fontsize=fs, color=st['label_color'], zorder=5, clip_on=False,
                bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                          edgecolor='none'))

    def plot_survival(self, profiles, times=None, labels=None, ci=False,
                      cl=None, bounds=None, simultaneous=False, target_bq=None,
                      km_overlay=False, extrap_style=None, show=True,
                      save=False, **kwargs):
        """Plot predicted survival curves S(t | x) for one or more covariate
        profiles (DataFrame / list of dicts / 2D array-like, one row each).

        ci=True adds the pointwise confidence band (bounds='fb' delta
        method, bounds='lrb' profile likelihood for AFT; defaults to the
        model's bounds= setting, cl overrides the level). Up to the last
        observed time the band is solid + filled; past it (extrapolation)
        it is dashed + hatched, flagged by a dotted vertical line at the
        last observed time carrying a small vertical "t_max" label near the
        top. simultaneous=True overlays the
        wider simultaneous band (valid over the whole observed time range)
        as dotted lines. The legend lists, top first, the bounds method
        and level ("Fisher bounds @ 90%" / "Likelihood-ratio bounds @
        90%"), then - with a dotted sample - "Simultaneous band" when
        simultaneous=True, then the profile curves.

        target_bq=p (e.g. 0.1 for B10) draws the horizontal line
        S = 1 - p and adds the B(100p) life to the legend as
        "B10: lower / median / upper" (just the median when ci=False).

        km_overlay=True adds the pooled Kaplan-Meier estimate of the whole
        sample as a grey step line - a model-free reference for the
        overall level of S(t)."""
        self._check_fitted()
        if times is None:
            hi_t = float(self._t.max()) * 1.02
            times = np.linspace(hi_t * 1e-4, hi_t, 200)
        times = np.atleast_1d(np.asarray(times, dtype=float))
        Xn = self._prep_X(profiles)
        n_prof = len(Xn)
        if labels is None:
            labels = ([str(i) for i in profiles.index]
                      if isinstance(profiles, pd.DataFrame)
                      else [f'profile {i + 1}' for i in range(n_prof)])

        band = None
        cl_ = method_ = None
        if ci or target_bq is not None:
            cl_, method_ = self._resolve_band_opts(cl, bounds)
        if ci:
            band = self.predict_survival(profiles, times=times, ci=True, cl=cl,
                                         bounds=bounds, simultaneous=simultaneous)
            S_all = band['surv']
        else:
            S_all = self.predict_survival(profiles, times=times)

        _apply_plot_style(self.plot_style)
        fig = plt.figure(figsize=self.fig_size, num=next(_figure_counter))
        ax = plt.gca()
        colors, linestyles = _categorical_style(n_prof)
        drawstyle = 'steps-post' if self.model == 'cox_ph' else 'default'
        t_hi = float(self._t.max())
        ext = times > t_hi                       # extrapolation (past last data)
        step = 'post' if drawstyle != 'default' else None

        def _fmt(v):
            return ('>%.3g' % t_hi) if not np.isfinite(v) else ('%.3g' % v)

        bq_marks = []
        prof_handles = []
        for j in range(n_prof):
            col = colors[j]
            S = S_all.iloc[:, j].to_numpy()
            lab = labels[j]
            if target_bq is not None:
                if ci:
                    u, m, o = self._bq_life(Xn[j], target_bq, cl_, method_)
                    lab = (f'{lab} — B{target_bq * 100:g}: '
                           f'{_fmt(u)} / {_fmt(m)} / {_fmt(o)}')
                else:
                    m = float(np.atleast_1d(
                        self.predict_quantile(Xn[j][None, :], q=target_bq)
                        if self.model == 'weibull_aft'
                        else self._cox_quantile_ci(Xn[j], target_bq,
                                                   self.cl)[1])[0])
                    lab = f'{lab} — B{target_bq * 100:g}: {_fmt(m)}'
                bq_marks.append((m, col))
            prof_handles.append(ax.plot(
                times, S, color=col, linestyle=linestyles[j],
                drawstyle=drawstyle, linewidth=1.8, label=lab)[0])
            if ci:
                lo = band['lower'].iloc[:, j].to_numpy()
                hi = band['upper'].iloc[:, j].to_numpy()
                # one continuous fill, then a hatch overlay on the
                # extrapolated part (interpolate so there is no gap)
                ax.fill_between(times, lo, hi, color=col, alpha=0.15,
                                linewidth=0.0, step=step)
                if ext.any():
                    ax.fill_between(times, lo, hi, where=ext, facecolor='none',
                                    edgecolor=col, hatch='////', alpha=0.45,
                                    linewidth=0.0, step=step, interpolate=True)
                em = ext.copy()
                if ext.any():
                    em[max(int(np.argmax(ext)) - 1, 0)] = True
                for arr in (lo, hi):
                    ax.plot(times, arr, color=col, linewidth=0.9,
                            linestyle='-', drawstyle=drawstyle)
                    if ext.any():
                        ax.plot(np.where(em, times, np.nan),
                                np.where(em, arr, np.nan), color=col,
                                linewidth=0.9, linestyle='--',
                                drawstyle=drawstyle)
                if band['lower_sim'] is not None:
                    ax.plot(times, band['lower_sim'].iloc[:, j].to_numpy(),
                            color=col, linewidth=1.0, linestyle=':',
                            drawstyle=drawstyle)
                    ax.plot(times, band['upper_sim'].iloc[:, j].to_numpy(),
                            color=col, linewidth=1.0, linestyle=':',
                            drawstyle=drawstyle)

        if km_overlay:
            tk, Sk = self._kaplan_meier(self._t, self._event)
            prof_handles.append(ax.plot(
                tk, Sk, drawstyle='steps-post', color='#8b8b8b',
                linewidth=1.3, label='Kaplan-Meier (pooled)')[0])

        # legend, ordered top -> bottom: bounds method, then (dotted
        # sample) the simultaneous band, then the profile curves
        head_handles = []
        if ci:
            bname = {'fb': 'Fisher bounds',
                     'lrb': 'Likelihood-ratio bounds',
                     'npbb': 'Non-parametric bootstrap bounds',
                     'pbb': 'Parametric bootstrap bounds'}[method_]
            head_handles.append(mpl.lines.Line2D(
                [], [], linestyle='none', marker='',
                label=f'{bname} @ {cl_ * 100:g}%'))
            if simultaneous:
                head_handles.append(mpl.lines.Line2D(
                    [], [], color='#8b8b8b', linewidth=1.0, linestyle=':',
                    label='Simultaneous band'))
        if target_bq is not None:
            ax.axhline(1.0 - target_bq, color='#8b8b8b', linestyle='--',
                       linewidth=1, zorder=1)
            for m, col in bq_marks:
                if np.isfinite(m):
                    ax.plot([m], [1.0 - target_bq], marker='o', color=col,
                            markersize=5, zorder=4)

        ax.set_xlabel(f'{self.x_label}'
                      f'{" in " + self.unit if self.unit != "-" else ""}',
                      color='black', fontsize=self.xy_fontsize)
        ax.set_ylabel('Survival probability', color='black',
                      fontsize=self.xy_fontsize)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        ax.set_title('Predicted Survival', color='black',
                     fontsize=self.plot_title_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(True, which='major')
        if self.show_legend:
            ax.legend(handles=head_handles + prof_handles,
                      fontsize=self.legend_fontsize)
        plt.tight_layout()
        # drawn last, after tight_layout, so the axes box (and thus the
        # room left of the right edge) is final when the label is sized
        if ci and ext.any():
            self._draw_extrapolation_marker(ax, t_hi, float(times.max()),
                                            **(extrap_style or {}))
        if save:
            try:
                plt.savefig(kwargs['path'])
            except Exception:
                raise ValueError('Path is faulty.')
        return self._render_plot(fig, show, save)


if __name__ == '__main__':
    data1 = np.random.normal(loc=10.0, scale=2.0, size=15)
    data2 = np.random.normal(loc=6.0, scale=1.0, size=12)
    n1 = Analysis(df=data1, dist='lognormal', bounds='lrb'); n1.mle()
    n2 = Analysis(df=data2, dist='lognormal', bounds='lrb'); n2.mle()
    PlotAll(objects={'Los 1': n1, 'Los 2': n2}).mult_lognormal(plot_ranks=False)
