# Available classes
Currently, there are two classes (Analysis and PlotAll) available in the predictr package. I will continue to add new classes in near future.
## Analysis
Analysis contains all necessary methods for the Weibull analysis.  
### Default arguments and values
This table provides information on alle arguments that are passed to the Analysis class.

|  Parameter  | default value |      type      |                                                description                                               |
|:-----------:|:-------------:|:--------------:|:--------------------------------------------------------------------------------------------------------:|
| df          |      None     | list of floats | List of failures                                                                                         |
| ds          |      None     | list of floats | List of suspensions (right-censored only)                                                             |
| bounds      |      None     |       str      | Confidence bounce method to be used in mle() or mrr()                                                 |
| bounds_type |      None     |       str      | Setting for the bounds: either two-sided or one-sided                                                 |
| show        |      False    |      bool      | If True, the Weibull probability plot will be plotted                                                 |
| bcm         |      None     |       str      | Defines the bias-correction method in mle()                                                           |
| cl          |      0.9      |      float     | Sets the confidence level when bounds are used                                                        |
| bs_size     |      5000     |       int      | Number of bootstrap samples                                                                              |
| est_type    |    'median'   |       str      | Sets the statistic to compute from the bootstrap samples                                             |
| plot_style  |    'ggplot'   |       str      | Choose a style according to your needs. See matplotlib style references for more available styles. |
| unit        |      '-'      |       str      | Unit of failures and suspensions, e.g. 's', 'ms', 'no. of cycle' etc.                                 |

**Important**:

- df = None will raise an error. There has to be at least one failure.

### Parameter estimation methods
One can either use the Maximum Likelihood Estimation or Median Rank Regression.

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
| C4 aka 'reduced bias adjustment'    |   x   |   -   |      'c4'      |    -    |                 -                |
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
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored.png){: width="500" }

##### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='lrb',show=True)
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_LRB_censored.png){: width="500" }

#### Median Rank Regression (MRR)
##### Uncensored sample
Example: 
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='bbb',show=True)
prototype_a.mrr()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_BBB_uncensored.png){: width="500" }

##### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='mcpb',show=True)
prototype_a.mrr()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_MCPB_censored.png){: width="500" }

#### Bias-corrections
As already mentioned, only mle() support bias-corrections. The samples in these examples are drawn from a two-parameter Weibull distribution with a shape parameter of 2.0 and a scale parameter of 1.0.

##### Uncensored sample
It is appearent that the estimates of beta and eta are now closer to the ground truth values. The dotted grey line in the plot is the "biased" MLE line, the bia-corrected line is blue. The legend contains all needed information.

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fb', show=True, bcm='c4')
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored_c4.png){: width="500" }

The estimates can for the Weibull parameters can be compared directly, since they are available as attributes
```python
print(f'biased beta: {prototype_a.beta:4f} --> bias-corrected beta: {prototype_a.beta_c4:4f}')
>>> biased beta: 2.511134 --> bias-corrected beta: 2.108248
```

##### Censored sample
The data is type II right-censored.
```python
failures = [0.38760099164906514, 0.5867052007217437, 0.5878056753744406, 0.602290402929083, 0.6754829518358306, 0.7520219855697948]
suspensions = [0.7520219855697948, 0.7520219855697948]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='lrb', show=True, bcm='hrbu')
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_LRB_censored_hrbu.png){: width="500" }

## PlotAll
PlotAll plots class objects from Analysis in one figure. Currently, only data from mle() is supported.
Theoretically, you can plot as many objects as you like -> provide a list of colors as a kwarg in PlotAll(objects, **kwargs).mult_weibull(). <b>
For now, six colors are supported by default, but you can pass an infinit amount of colors to the mult_weibull() method.

**Available methods**:

| Methods        	| Description                                                           	|
|----------------	|-----------------------------------------------------------------------	|
| mult_weibull() 	| Plots multiple Analysis class instances in one Weibull plot           	|
| contour_plot() 	| Plots contour plots when likelihood ratio bounds are used in Analysis 	|

### mult_weibull()
#### Both with two-sided bounds - default colors
```python
from predictr import Analysis, PlotAll

# Create new objects, e.g. name them prototype_a and prototype_b
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='lrb', bounds_type='2s')
prototype_a.mle()

failures_b = [1.8506941739639076, 2.2685555679846954, 2.380993183650987, 2.642404955035375,
              2.777082863078587, 2.89527127055147, 2.9099992138728927, 3.1425481097241,
              3.3758727398694406, 3.8274990886889997]
prototype_b = Analysis(df=failures_b, bounds='pbb', bounds_type='2s')
prototype_b.mle()

# Create dictionary with Analysis objects
# Keys will be used in figure legend. Name them as you please.
objects = {'proto_a': prototype_a, 'proto_b': prototype_b}

# Use mult_weibull() method
PlotAll(objects).mult_weibull()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_MLE_2s.png){: width="500" }

#### One object with a one-sided lower bound, the other one has two-sided bounds - default colors
You can plot every bounds_type ('2s', '1sl', '1su') and combine them.
```python
from predictr import Analysis, PlotAll

# Create new objects, e.g. name them prototype_a and prototype_b
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='fb', bounds_type='1sl')
prototype_a.mle()

failures_b = [1.8506941739639076, 2.2685555679846954, 2.380993183650987, 2.642404955035375,
              2.777082863078587, 2.89527127055147, 2.9099992138728927, 3.1425481097241,
              3.3758727398694406, 3.8274990886889997]
prototype_b = Analysis(df=failures_b, bounds='npbb', bounds_type='2s')
prototype_b.mle()

# Create dictionary with Analysis objects
# Keys will be used in figure legend. Name them as you please.
objects = {'proto_a': prototype_a, 'proto_b': prototype_b}

# Use mult_weibull() method
PlotAll(objects).mult_weibull()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_MLE_1sl_2s.png){: width="500" }

#### Three objects - custom colors
```python
from predictr import Analysis, PlotAll

# Create new objects
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='fb', bounds_type='2s')
prototype_a.mle()

failures_b = [1.8506941739639076, 2.2685555679846954, 2.380993183650987, 2.642404955035375,
              2.777082863078587, 2.89527127055147, 2.9099992138728927, 3.1425481097241,
              3.3758727398694406, 3.8274990886889997]
prototype_b = Analysis(df=failures_b, bounds='npbb', bounds_type='2s')
prototype_b.mle()

failures_c = [0.04675399107295282, 0.31260891592041457, 0.32121232576015757, 0.6013488316204837,
              0.7755159796641791, 0.8994041575114923, 0.956417788622185, 1.1967354178170764,
              1.6115311492838604, 2.1120891587523793]
prototype_c = Analysis(df=failures_c, bounds='pbb', bounds_type='2s')
prototype_c.mle()

objects = {'proto_a': prototype_a, 'proto_b': prototype_b, 'secret': prototype_c}
# Create list with custom colors and pass to the instance method with the kwarg: set_cmap = [...]
colors = ['green', 'red', 'blue']
PlotAll(objects, set_cmap = colors).mult_weibull()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_MLE_2s_custom_colors.png){: width="500" }

### contour_plot()
contour_plot() only works for likelihood ratio bounds. Hence, you have to use bounds='lrb' in the Analysis class. This method supports all bounds types and all confidence levels. You can pass as many objects as you want to.

#### Plot a single Analysis object
```python
from predictr import Analysis, PlotAll

# Create new objects
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='lrb', bounds_type='2s')
prototype_a.mle()

objects = {'initial design': prototype_a}
PlotAll(objects).contour_plot()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Contour_plot_LRB.png){: width="500" }

#### Plot a multiple Analysis objects
```python
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='lrb', bounds_type='2s')
prototype_a.mle()

failures_c = [0.04675399107295282, 0.31260891592041457, 0.32121232576015757, 0.6013488316204837,
              0.7755159796641791, 0.8994041575114923, 0.956417788622185, 1.1967354178170764,
              1.6115311492838604, 2.1120891587523793]
prototype_c = Analysis(df=failures_c, bounds='lrb', bcm = 'hrbu', bounds_type='2s')
prototype_c.mle()

# Create dictionary with Analysis objects
# Keys will be used in figure legend. Name them as you please.
objects = {'initial design': prototype_a, 'final design': prototype_c}
PlotAll(objects).contour_plot()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Contour_plot_LRB_multiple.png){: width="500" }
