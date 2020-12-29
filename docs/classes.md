# Available classes
Currently, there is only one class available in the predictr package. I will continue to add new classes in the near future.
## 1. Analysis
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
- mle() and mrr() support only specific confidence bounds methods. For instance, you can't use Beta-Binomial Bounds with mle(). This will also raise an error. Use the table below to check, whether a combination of parameter estimation and confidence bounds method is supported.

| confidence bounds               | mle() | mrr() | uncensored data | censored data |        type        | argument value |
|---------------------------------|:-----:|:-----:|:---------------:|:-------------:|:------------------:|:--------------:|
| Beta-Binomial Bounds            |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |      'bbb'     |
| Monte-Carlo Pivotal Bounds      |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |     'mcpb'     |
| Non-Parametric Bootstrap Bounds |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |     'npbb'     |
| Parametric Bootstrap Bounds     |   x   |   x   |        x        |       -       | '2s', '1sl', '1su' |      'pbb'     |
| Fisher Bounds                   |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'fb'      |
| Likelihood Ratio Bounds         |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'lrb'     |

Following parameter estimation methods are supported:
**Maximum likelihood estimation (MLE):** 
```python
from predictr import Analysis
prototype_a = Analysis(...) # create an instance
prototype_a.mle() # use instance methods
```
**Median Rank Regression**
```python
from predictr import Analysis
prototype_a = Analysis(...) # create an instance
prototype_a.mrr() # use instance methods
```
