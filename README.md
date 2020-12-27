# predictr
Weibull Analysis Utilities

![](https://img.shields.io/pypi/v/predictr?color=blue&label=pypi)
![](https://img.shields.io/github/stars/tvtoglu/predictr?style=flat-square)

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install predictr.

```bash
pip install predictr
```

## Usage
### Default Parameter values
df: list = None -> failures, e.g. df = [100, 120, 80, 300]<br>

ds: list = None -> suspensions (right-censored), ds = [300, 400, 400]<br>

show: bool = False -> If True, the Weibull probability plot will be plotted.<br>
plot_style = 'ggplot' -> Choose a style according to your needs. See https://matplotlib.org/3.1.0/gallery/style_sheets/style_sheets_reference.html for styles.<br>
bounds=None -> Use following table to configure everything related to confidence bounds, e.g. if you want to use Monte-Carlo pivotal bounds for the Median Rank Regression: bounds = 'mcpb'.<br>

| confidence bounds               | mle() | mrr() | uncensored data | censored data |        type        | argument value |
|---------------------------------|:-----:|:-----:|:---------------:|:-------------:|:------------------:|:--------------:|
| Beta-Binomial Bounds            |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |      'bbb'     |
| Monte-Carlo Pivotal Bounds      |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |     'mcpb'     |
| Non-Parametric Bootstrap Bounds |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |     'npbb'     |
| Parametric Bootstrap Bounds     |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |      'pbb'     |
| Fisher Bounds                   |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |    'fisher'    |
| Likelihood Ratio Bounds         |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'lrb'     |

bounds_type = '2s' -> '2s': two-sided confidence bounds, '1su': upper confidence bounds, '1sl': lower confidence bounds. E.g. bounds_type = '1sl'.<br>
cl=0.9 -> configure the confidence level in the intervall (0, 1.0)
bcm=None, 
bs_size=5000, 
est_type='median',
                 unit='-'
### To use the Maximum Likelihood Estimation (MLE)
#### Uncensored sample
example:
sample = 
```python
import predictr


```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
