# Available classes
Currently, there is only one class available in the predictr package. I will continue to add new classes in the near future.
## Analysis
Analysis contains all necessary methods for the Weibull analysis.  
### Default arguments and values
This table provides information on alle arguments that are passed to the Analysis class.

|  Parameter  | default value |      type      |                                                description                                               |
|:-----------:|:-------------:|:--------------:|:--------------------------------------------------------------------------------------------------------:|
| df          |      None     | list of floats | List of failures                                                                                         |
| ds          |      None     | list of floats | List of suspensions<br>(right-censored only)                                                             |
| bounds      |      None     |       str      | Confidence bounce method to<br>be used in mle() or mrr()                                                 |
| bounds_type |      None     |       str      | Setting for the bounds:<br>either two-sided or one-sided                                                 |
| show        |      None     |      bool      | If True, the Weibull probability<br>plot will be plotted                                                 |
| bcm         |      None     |       str      | Defines the bias-correction<br>method in mle()                                                           |
| cl          |      0.9      |      float     | Sets the confidence level<br>when bounds are used                                                        |
| bs_size     |      5000     |       int      | Number of bootstrap samples                                                                              |
| est_type    |    'median'   |       str      | Sets the statistic to compute <br>from the bootstrap samples                                             |
| plot_style  |    'ggplot'   |       str      | Choose a style according to your needs.<br>See matplotlib style references for more<br>available styles. |
| unit        |      '-'      |       str      | Unit of failures and suspensions, e.g.<br>'s', 'ms', 'no. of cycle' etc.                                 |

**Important**:

- df = None will raise an error. There has to be at least one failure.

### Parameter estimation methods
One can use either the Maximum Likelihoof Estimation or Median Rank Regression.

**Maximum likelihood estimation (MLE):** 
```python
from predictr import Analysis
prototype_a = Analysis(...) # create an instance
prototype_a.mle() # use instance methods
```
**Median Rank Regression (MRR)**
```python
from predictr import Analysis
prototype_a = Analysis(...) # create an instance
prototype_a.mrr() # use instance methods
```
### Bias-correction methods
Since parameter estimation methods are only asymptotically unbiased (sample sizes -> "infinity"), bias-correction methods are useful when you have only a few failures. These methods correct the Weibull shape and scale parameter.
The following table provides possible configurations. Bias-corrections for mrr() are not supported, yet.<br>

| Bias-correction method              | mle() | mrr() | argument value | config. |             statistic            |
|-------------------------------------|:-----:|:-----:|:--------------:|:-------:|:--------------------------------:|
| C4 aka 'reduced bias adjustment     |   x   |   -   |      'c4'      |    -    |                 -                |
| Hirose and Ross method              |   x   |   -   |     'hrbu'     |    -    |                 -                |
| Non-parametric Bootstrap correction |   x   |   -   |     'np_bs'    | bs_size | 'mean', 'median', 'trimmed_mean' |
| Parametric Bootstrap correction     |   x   |   -   |     'p_bs'     | bs_size | 'mean', 'median', 'trimmed_mean' |

### Confidence bounds methods
Analysis supports nearly all state of the art confidence bounds methods.

| confidence bounds               | mle() | mrr() | uncensored data | censored data |    bounds_type     | argument value |
|---------------------------------|:-----:|:-----:|:---------------:|:-------------:|:------------------:|:--------------:|
| Beta-Binomial Bounds            |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |      'bbb'     |
| Monte-Carlo Pivotal Bounds      |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |     'mcpb'     |
| Non-Parametric Bootstrap Bounds |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |     'npbb'     |
| Parametric Bootstrap Bounds     |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |      'pbb'     |
| Fisher Bounds                   |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'fb'      |
| Likelihood Ratio Bounds         |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'lrb'     |

**Important**:
- mle() and mrr() support only specific confidence bounds methods. For instance, you can't use Beta-Binomial Bounds with mle(). This will also raise an error. Use the table below to check, whether a combination of parameter estimation and confidence bounds method is supported.
- '2s': two-sided confidence bounds, '1su': upper confidence bounds, '1sl': lower confidence bounds. If Beta-Binomial Bounds are used, the lower bound represents the lower percentile bound at a specific time ((pctl) is added in the plot legend). If Fisher Bounds are used, the lower bound represents the lower time bound at a specific percentile.

