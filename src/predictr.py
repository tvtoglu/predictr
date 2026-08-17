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
from scipy.special import gamma
from scipy.stats import norm, chi2, beta, linregress, trim_mean
from scipy.stats.distributions import weibull_min
from scipy.spatial import ConvexHull


# predictr's categorical palette (up to 6 clearly distinct series) and the
# linestyles used to extend it beyond 6 - see _categorical_style() below.
PREDICTR_PALETTE = ['#008b8b', '#008fd5', '#810f7c', '#8b8b8b', '#fc4f30', '#e5c494']
PREDICTR_LINESTYLES = ['-', '--', ':', '-.']


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


# predictr's default matplotlib style. Kept in sync with src/predictr.mplstyle
# (the human-editable reference copy) but embedded here as a plain dict so
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


class Analysis:
    """
    Analysis provides parameter estimations, confidence bounds
    computations, bias corrections, and plotting of the data.
    """

    # Distributions Analysis can be configured for via dist=. Only 'weibull'
    # is implemented so far; this set is the single place future
    # distributions get registered as they're added.
    SUPPORTED_DISTRIBUTIONS = {'weibull', 'normal'}

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
        """
        return norm.pdf(z) / (1 - norm.cdf(z))

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
                 plot_title=None, plot_ranks=True,
                 fig_size=(6, 7), show_legend=True, legend_fontsize=9, save=False, **kwargs):
        """
        Parameters
        ----------
        df : list
            Contains failures. The default is None.
        ds : list
            Contains suspensions. The default is None.
        dist : string, optional
            Sets the distribution to fit: 'weibull' or 'normal'. More
            distributions will be added over time. The default is
            'weibull'. For 'normal', only bounds='fb' (Fisher bounds) and
            bcm=None are currently supported.
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

        # Plot related attributes
        self.show = show
        self.plot_style = plot_style
        if plot_title is None:
            plot_title = ('Weibull Probability Plot' if dist == 'weibull'
                           else f'{dist.capitalize()} Probability Plot')
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
        self.beta_hrbu, self.eta_hrbu = None, None
        self.f, self.f_inv= None, None
        self.k_a_bound, self.se_beta,self.se_eta = None, None, None
        self.se_mu, self.se_sigma = None, None
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
            and self.bounds != 'npbb'):
            raise ValueError(f'"{self.bounds}" is not supported by mle')
        if self.dist == 'normal' and self.bounds not in (None, 'fb'):
            raise ValueError(f'"{self.bounds}" is not yet supported for '
                              f'dist="normal". Only bounds="fb" is '
                              f'currently implemented.')
        if self.dist == 'normal' and self.bcm is not None:
            raise ValueError('Bias-correction methods (bcm) are not yet '
                              'supported for dist="normal".')
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

        elif self.dist == 'normal':
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

            if self.ds is None:
                # Closed-form solution of the score equations above with
                # no suspension terms.
                data = np.asarray(self.df, dtype=float)
                self.mu = data.mean()
                self.sigma = np.sqrt(np.mean((data - self.mu) ** 2))
            else:
                self.censoring = 'right-censored'
                mu_init = float(np.mean(self.df))
                sigma_init = float(np.std(self.df)) if len(self.df) > 1 else 1.0
                self.mu, self.sigma = optimize.fsolve(
                    normal_score, x0=[mu_init, sigma_init],
                    fprime=normal_score_jac, args=(self.df, self.ds))
                if self.sigma <= 0:
                    raise ValueError('Sigma estimation is not valid. '
                                     'Check the input data.')

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
        if self.bounds == 'fb':
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
        plt.figure(figsize=self.fig_size)

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
                            0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                            0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.ylim(bottom=0.001, top=0.999)
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
                     color='mediumblue', linestyle='-',
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
                                 yerr=[yerr_lower, yerr_upper], fmt='none', ecolor= 'royalblue',
                                 capsize=5, markeredgewidth=2, alpha=0.5)
                    plt.fill_between(x = self.df,
                                     y1=weibull_prob_paper(self.bounds_lower),
                                     y2=weibull_prob_paper(self.bounds_upper),
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 ecolor= 'royalblue', capsize=5, markeredgewidth=2, alpha=0.5)
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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                                 ecolor= 'royalblue', capsize=5, markeredgewidth=2, alpha=0.5)
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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                                 color='royalblue', linestyle='-', linewidth=1)

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
                             markerfacecolor='mediumblue', markeredgecolor='mediumblue',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)
            else:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank_cens()), marker='o',
                             markerfacecolor='mediumblue', markeredgecolor='mediumblue',
                             markersize=4, alpha=.5, linestyle='None', zorder= 3)

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
            fig = plt.gcf()
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

        if self.show:
            plt.show()

    def fisher_bounds(self):
        """
        Computes Fisher bounds for the fitted distribution (self.dist).
        """
        if self.dist == 'normal':
            self._fisher_bounds_normal()
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

    def lrb(self):
        """
        # Goal: Find all solution pairs (beta, eta) for
        # L(beta, eta) = exp(chi ** 2) / -2) * L(beta_mle, eta_mle)
        """

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
            self._plot_normal()
            return

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
        plt.figure(figsize=self.fig_size)

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
                            0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                            0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.ylim(bottom=0.001, top=0.999)
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
                         color='mediumblue', linestyle='-',
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                         color='mediumblue', linestyle='-',
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
                                 color='royalblue',
                                 linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                         color='mediumblue', linestyle='-',
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                         color='mediumblue', linestyle='-',
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                         color='mediumblue', linestyle='-',
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
                                 color='royalblue', linestyle='-', linewidth=1)
                    plt.semilogx(self.bounds_upper, weibull_prob_paper(self.unrel),
                                 color='royalblue', linestyle='-', linewidth=1, label='_nolegend_')
                    plt.fill_betweenx(y=weibull_prob_paper(self.unrel),
                                     x1=self.bounds_lower,
                                     x2=self.bounds_upper,
                                     alpha=0.1, color = 'royalblue', label='_nolegend_')
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                                 color='royalblue', linestyle='-', linewidth=1)
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
                             markerfacecolor='mediumblue', markeredgecolor='mediumblue',
                             markersize=4, alpha=.5, linestyle='None', zorder=3)
            else:
                plt.semilogx(self.df, weibull_prob_paper(self.median_rank_cens()), marker='o',
                             markerfacecolor='mediumblue', markeredgecolor='mediumblue',
                             markersize=4, alpha=.5, linestyle='None', zorder= 3)

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
            fig = plt.gcf()
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

        if self.show:
            plt.show()

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
        plt.figure(figsize=self.fig_size)

        # Y-axis: probit scale, same tick percentages as the Weibull plot
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(normal_ticks))
        y_ticks = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
                            0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                            0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        z_ticks = norm.ppf(y_ticks)
        plt.ylim(bottom=z_ticks[0], top=z_ticks[-1])
        plt.yticks(z_ticks, color='black')
        ax.set_yticks([0.0], minor=True)
        plt.grid(True, which='minor', axis='y', linestyle='--')

        # PREDICTR_PALETTE (the categorical 6-color palette) is for telling
        # multiple *different* Analysis objects apart, e.g. in PlotAll -
        # not for the elements of a single result. A single fit instead
        # uses predictr's own single-result convention, exactly like
        # Weibull's plot() above: mediumblue for the MLE line/markers,
        # royalblue (thinner, plus a light fill) for the bounds.
        x_data = np.array(self.df, dtype=float)
        susp_num = len(self.ds) if self.ds is not None else 0

        # Fitted MLE line - directly linear in x, no inverse-CDF sampling
        # needed the way Weibull's inverse_weibull() draws its line.
        pad = 0.2 * (x_data.max() - x_data.min()) if x_data.max() > x_data.min() else max(abs(x_data.min()), 1.0)
        x_line = np.linspace(x_data.min() - pad, x_data.max() + pad, 200)
        y_line = (x_line - self.mu) / self.sigma
        plt.plot(x_line, y_line, color='mediumblue', linestyle='-',
                 linewidth=1.5, zorder=2)

        leg_title = 'MLE'
        leg_text = ('n = {} (f: {} | s: {})\n'.format(len(self.df) + susp_num,
                                                        len(self.df), susp_num)
                    + r'$\widehat\mu={:.3f}$ | '.format(self.mu)
                    + r'$\widehat\sigma={:.3f}$'.format(self.sigma))
        legend_labels = (leg_text,)

        if self.bounds == 'fb' and (self.bounds_lower is not None
                                     or self.bounds_upper is not None):
            z_p = norm.ppf(self.unrel)
            if self.bounds_type == '2s':
                plt.plot(self.bounds_lower, z_p, color='royalblue',
                         linestyle='-', linewidth=1)
                plt.plot(self.bounds_upper, z_p, color='royalblue',
                         linestyle='-', linewidth=1, label='_nolegend_')
                plt.fill_betweenx(y=z_p, x1=self.bounds_lower, x2=self.bounds_upper,
                                   alpha=0.1, color='royalblue', label='_nolegend_')
                bt_legend = '2s'
            elif self.bounds_type == '1su':
                plt.plot(self.bounds_upper, z_p, color='royalblue',
                         linestyle='-', linewidth=1)
                bt_legend = '1su'
            elif self.bounds_type == '1sl':
                plt.plot(self.bounds_lower, z_p, color='royalblue',
                         linestyle='-', linewidth=1)
                bt_legend = '1sl'
            legend_labels = (leg_text, '\nFisher bounds:\n{} @{}%'.format(
                bt_legend, self.cl * 100))

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
            plt.plot(x_data, y_data, marker='o', markerfacecolor='mediumblue',
                     markeredgecolor='mediumblue', markersize=4, alpha=.5,
                     linestyle='None', zorder=3)

        plt.tight_layout()

        # Same legend-overflow-safety as plot()/plot_mrr() above: widen the
        # figure if the legend box would otherwise get clipped at the right
        # edge. See the comment on the first occurrence of this pattern
        # (plot()) for why it must run after tight_layout().
        legend = plt.gca().get_legend()
        if legend is not None:
            fig = plt.gcf()
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
        if self.show:
            plt.show()

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
                     **kwargs):
        """
        Plots multiple Analysis class objects in one figure

        """

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
        plt.figure(figsize=fig_size)

        # Y-Axis
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(weibull_ticks))
        y_ticks = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
                            0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                            0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        lny_ticks = np.log(-np.log(1 - y_ticks))
        plt.ylim(bottom=0.001, top=0.999)
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
        plt.tight_layout()
        plt.grid(True, which='both')
        plt.legend(fontsize= legend_fontsize)

        # Save plot
        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

        plt.show()

    def contour_plot(self, show=True, style='hull', show_weibull=False, show_legend=True, color=None, x_label=r'$\widehat\beta$',
                     y_label=None, plot_title='Contour Plot', xy_fontsize=12,
                     plot_title_fontsize=14, legend_fontsize=9, fig_size=(6.4, 4.8), save=False,
                     scale_mode='auto', log_ratio_threshold=10, cl_set=None,
                     curve_fill=True, fill_alpha=0.25, **kwargs):
        """
        Plots the contour plot when likelihood ratio bounds are being used.
        Multiple objects can be used as well.

        Parameters
        ----------
        scale_mode : {'auto', 'linear', 'log'}, optional
            Scaling of the eta (y) axis. 'auto' (default) inspects the eta_lrb
            range across all objects and switches to a logarithmic scale when
            it spans more than `log_ratio_threshold`x in magnitude. This
            avoids small-scale ellipses being crushed into a flat line next
            to large-scale ones when several objects with very different eta
            magnitudes are plotted together. 'linear' and 'log' force the
            respective scale regardless of the data.
        log_ratio_threshold : float, optional
            Ratio of max(eta) / min(eta) across all objects above which
            'auto' switches to a logarithmic scale. Default = 10.
        cl_set : list of float, optional
            Confidence levels to draw per dataset, e.g. [0.95, 0.9, 0.8].
            If None or empty (default), each object's own `cl` attribute is
            used, i.e. exactly one curve per object, unchanged from before.
            When given, `cl_set` is sorted from largest to smallest and, for
            every object, its contour is recomputed and drawn once per
            confidence level - via a temporary deep copy of the object, so
            the original object's own `cl`/`beta_lrb`/`eta_lrb` are left
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

        # Configure plot
        _apply_plot_style(self.plot_style)
        fig, ax = plt.subplots(figsize=fig_size)
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

        # Collect the (cl, beta_lrb, eta_lrb) curve(s) to draw per object. If
        # cl_set is given, the LRB region is recomputed for every requested
        # confidence level on a deep copy of the object, leaving the
        # original object untouched; otherwise the object's own cl/beta_lrb/
        # eta_lrb (as computed by .mle()/.lrb()) is used, as before.
        curves = {}
        for key, val in self.objects.items():
            if cl_set:
                object_curves = []
                for cl_value in cl_set:
                    val_cl = copy.deepcopy(val)
                    val_cl.cl = cl_value
                    val_cl.lrb()
                    object_curves.append((cl_value, val_cl.beta_lrb, val_cl.eta_lrb))
                curves[key] = object_curves
            else:
                curves[key] = [(getattr(val, 'cl'), getattr(val, 'beta_lrb'), getattr(val, 'eta_lrb'))]

        # Analyze the eta_lrb data across all objects/confidence levels to decide on the y-axis scale
        all_eta = np.concatenate([np.asarray(eta) for object_curves in curves.values()
                                  for _, _, eta in object_curves])
        if scale_mode == 'log':
            use_log = True
        elif scale_mode == 'linear':
            use_log = False
        else:
            use_log = (all_eta.max() / all_eta.min()) > log_ratio_threshold
        if use_log and np.any(all_eta <= 0):
            raise ValueError('Logarithmic eta scale requires strictly positive eta_lrb values.')

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

            # Point estimate marker (bias-corrected if a bcm was used), always
            # one tick darker (t=1.0) than the dataset's darkest curve (which
            # tops out at CURVE_T_MAX < 1.0)
            point_beta, point_eta = Analysis._point_estimate(val)
            ax.scatter(point_beta, point_eta, s=40,
                        c=[shade_of(base_color, 1.0)], marker='o', zorder=5)

        ax.set_xlabel(x_label, fontsize=xy_fontsize)
        if y_label is None:
            y_label = r'$\log(\widehat\eta)$' if use_log else r'$\widehat\eta$'
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
            return fig
        if show:
            plt.show()


    def weibull_pdf(self, beta=None, eta=None, linestyle=['-', '--', ':', '-.'], labels = None,
                    x_label = None, y_label=None, xy_fontsize=12, tick_fontsize=10,
                    legend_fontsize=9,
                    plot_title='Weibull PDF', plot_title_fontsize=14, x_bounds=None,
                    fig_size=None, color=None, save=False, plot_style='predictr', **kwargs):
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
        plt.figure(figsize=fig_size)

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

        plt.show()

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
        x.plot()

        # Save plot
        if save:
            try:
                plt.savefig(kwargs['path'])
            except:
                raise ValueError('Path is faulty.')

if __name__ == '__main__':
    failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]

    #
    #PlotAll(objects).contour_plot()
    
    contour_decision = False
    objects = {}
    if contour_decision == True:
        for cl in [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
            key = f"variable_{str(cl).replace('.', '_')}"
            objects[key] = Analysis(df=failures_a, bounds='lrb', bounds_type='2s', show=False, cl=cl)
            objects[key].mle()
        PlotAll(objects).contour_plot()


    df = [0.670659, 0.976145, 1.41494, 0.859942, 0.468364, 1.17272, 0.648734, 0.972926, 0.851652, 1.08389]
    df_caf = [1400699, 45477, 49358, 53379, 70695, 74721, 116451]
    ds_caf = [3_000_000] * 3


    x = Analysis(df=df, bounds='fb', dist='normal', show=True)
    x.mle()
"""
    y = Analysis(df=failures_a, bounds='lrb')
    y.mle()

    obj = {'x': x, 'y':y }
    PlotAll(obj).contour_plot(curve_fill= True, scale_mode='auto', cl_set=[.7, 0.8, 0.9])

    failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
    prototype_a = Analysis(df=failures, bounds='bbb', show=True)
    prototype_a.mrr()
"""