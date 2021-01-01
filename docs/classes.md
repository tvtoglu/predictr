# Available classes
Currently, there is only one class available in the predictr package. I will continue to add new classes in near future.
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

- mle() and mrr() support only specific confidence bounds methods. For instance, you can't use Beta-Binomial Bounds with mle(). This will also raise an error. Use the table above to check, whether a combination of parameter estimation and confidence bounds method is supported.
- '2s': two-sided confidence bounds, '1su': upper confidence bounds, '1sl': lower confidence bounds. If Beta-Binomial Bounds are used, the lower bound represents the lower percentile bound at a specific time ((pctl) is added in the plot legend). If Fisher Bounds are used, the lower bound represents the lower time bound at a specific percentile.

### Examples
#### Maximum Likelihood Estimation (MLE)
##### Uncensored sample
Example: 
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fb',show=True)
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored.png" height="500" />

##### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='lrb',show=True)
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_LRB_censored.png" height="500" />

#### Median Rank Regression (MRR)
##### Uncensored sample
Example: 
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='bbb',show=True)
prototype_a.mrr()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_BBB_uncensored.png" height="500" />

##### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='mcpb',show=True)
prototype_a.mrr()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_MCPB_censored.png" height="500" />

#### Bias-corrections
As already mentioned, only mle() support bias-corrections. The samples in these examples are drawn from a two-parameter Weibull distribution with a shape parameter of 2.0 and a scale parameter of 1.0.

##### Uncensored sample
It is appearent that the estimates of beta and eta are now closer to the ground truth values. The dotted grey line in the plot is the "biased" MLE line, the bia-corrected line is blue. The legend contains all needed information.

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fb', show=True, bcm='c4')
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored_c4.png" height="500" />

The estimates can for the Weibull parameters can be compared directly, since they are available as attributes
```python
print(f'biased beta: {prototype_a.beta:4f} --> bias-corrected beta: {prototype_a.beta_c4:4f}')
```

##### Censored sample
The data is type II right-censored.
```python
failures = [0.38760099164906514, 0.5867052007217437, 0.5878056753744406, 0.602290402929083, 0.6754829518358306, 0.7520219855697948]
suspensions = [0.7520219855697948, 0.7520219855697948]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='lrb', show=True, bcm='hrbu')
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_LRB_censored_hrbu.png" height="500" />
