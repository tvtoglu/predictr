# predictr. 
predict + reliability, in other words: A tool to predict the reliability.  

**predictr** is a Python package for Weibull-based life data analysis (reliability engineering). It covers parameter estimation, bias-correction, confidence bounds, and publication-ready Weibull plots in a single, consistent API.

![](https://img.shields.io/pypi/v/predictr?color=blue&style=flat&label=pypi)
[![Downloads](https://pepy.tech/badge/predictr)](https://pepy.tech/project/predictr)
![](https://img.shields.io/pypi/pyversions/predictr)
![](https://img.shields.io/pypi/l/predictr)
![](https://img.shields.io/github/stars/tvtoglu/predictr?style=flat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21901236.svg)](https://doi.org/10.5281/zenodo.21901236)


## Installation

```bash
pip install predictr
```

Requires Python >= 3.6.

## Quick start

```python
from predictr import Analysis

failures = [0.4508831, 0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]

weibull = Analysis(df=failures, bounds='fb', show=True)
weibull.mle()

print(weibull.beta, weibull.eta)  # shape and scale estimates
```

This fits a two-parameter Weibull distribution via Maximum Likelihood Estimation, adds Fisher confidence bounds, and renders the probability plot below.

<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored.png" alt="Weibull probability plot with Fisher confidence bounds" width="380">

## Main features

**Parameter estimation**
- Uncensored and type I / type II right-censored two-parameter Weibull distribution
- Maximum Likelihood Estimation (MLE) and Median Rank Regression (MRR)
- Bx-life calculator

**Bias-correction**
- C4 method (reduced bias adjustment)
- Hirose and Ross method
- Parametric and non-parametric bootstrap correction (mean, median, trimmed mean)

**Confidence bounds**
- Fisher bounds
- Likelihood Ratio bounds
- Beta-Binomial bounds
- Monte Carlo Pivotal bounds
- Parametric and non-parametric bootstrap bounds

**Plots**
- Weibull probability plots with all relevant statistics in the legend
- Multiple Weibull fits overlaid in one figure, for design comparisons
- Contour plots for the joint confidence region of shape and scale

| Comparing multiple fits | Confidence region contour |
|:---:|:---:|
| <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_MLE_2s.png" alt="Multiple Weibull fits in one plot" width="260"> | <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Contour_plot_LRB.png" alt="Contour plot of shape and scale confidence region" width="260"> |

See the [class documentation](https://tvtoglu.github.io/predictr/classes/) for the full method and parameter reference, including censored-data and bias-correction examples.

## Documentation and links

- [Documentation](https://tvtoglu.github.io/predictr/)
- [PyPI](https://pypi.org/project/predictr/)
- [Changelog](https://tvtoglu.github.io/predictr/CHANGELOG/)
- [Citation / Zenodo](https://doi.org/10.5281/zenodo.4433164)
- [Discussions](https://github.com/tvtoglu/predictr/discussions)

## Citing predictr

If you use predictr in academic work, please cite it via its [Zenodo DOI](https://doi.org/10.5281/zenodo.4433164). See [docs/citation.md](https://tvtoglu.github.io/predictr/citation/) for details.

## License

MIT — see [LICENSE.txt](LICENSE.txt).

## Contacte me

If you have any questions and / or suggestions, don't hesitate to contact me.
