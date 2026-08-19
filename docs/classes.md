# Available classes
Currently, there are two classes (Analysis and PlotAll) available in the predictr package. I will continue to add new classes in near future.
## Analysis
Analysis contains all necessary methods for the Weibull analysis. Since version 0.1.34, it also supports dist='normal', dist='lognormal' and dist='exponential' - see [Distributions](#distributions) below.
### Default arguments and values
This table provides information on alle arguments that are passed to the Analysis class.

| Parameter           | default value              | type            | description                                                                                        |
|---------------------|----------------------------|-----------------|----------------------------------------------------------------------------------------------------|
| df                  | None                       | list of floats  | List of failures                                                                                   |
| ds                  | None                       | list of floats  | List of suspensions (right-censored only)                                                          |
| dist                | 'weibull'                  | str             | Distribution to fit: 'weibull', 'normal', 'lognormal' or 'exponential'                             |
| bounds              | None                       | str             | Confidence bounce method to be used in mle() or mrr()                                              |
| bounds_type         | None                       | str             | Setting for the bounds: either two-sided or one-sided                                              |
| show                | False                      | bool            | If True, the Weibull probability plot will be plotted                                              |
| bcm                 | None                       | str             | Defines the bias-correction method in mle()                                                        |
| cl                  | 0.9                        | float           | Sets the confidence level when bounds are used                                                     |
| bs_size             | 5000                       | int             | Number of bootstrap samples                                                                        |
| est_type            | 'median'                   | str             | Sets the statistic to compute from the bootstrap samples                                           |
| plot_style          | 'predictr'                 | str             | Choose a style according to your needs. 'predictr' is predictr's own built-in style (no setup required); see matplotlib style references for other available styles. |
| unit                | '-'                        | str             | Unit of failures and suspensions, e.g. 's', 'ms', 'no. of cycle' etc.                              |
| x_label             | 'Time to Failure'          | string          | Label for the x-axis                                                                               |
| y_label             | 'Unreliability'            | string          | Label for the y-axis                                                                               |
| xy_fontsize         | 12                         | float           | Fontsize for the axes label                                                                        |
| tick_fontsize       | 10                         | float           | Fontsize for the tick labels (the numbers on the axes)                                             |
| legend_fontsize     | 9                          | float           | Fontsize for the legend                                                                            |
| plot_title          | 'Weibull Probability Plot' | string          | Title for the plot                                                                                 |
| plot_title_fontsize | 14                         | float           | Fontsize of the plot title                                                                         |
| fig_size            | (6, 7)                     | tuple of floats | Sets figure width and height in inches: (width, height)                                            |
| save                | False                      | boolean         | the beta and eta length of lists.                                                                  |
| plot_ranks          | True                       | boolean         | If True, median ranks will be plotted.                                                             |
| show_legend         | True                       | boolean         | If True, the legend will be plotted                                                                |
| kwarg: path         |                            | string          | Path defines the directory and format of the figure E.g. r'var/user/.../test.pdf'                  |


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
Analysis supports nearly all state of the art confidence bounds methods. The table below applies to dist='weibull' (the default). For the other distributions, see [Distributions](#distributions).

| confidence bounds               | mle() | mrr() | uncensored data | censored data |    bounds_type     | argument value |
|---------------------------------|:-----:|:-----:|:---------------:|:-------------:|:------------------:|:--------------:|
| Beta-Binomial Bounds            |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |      'bbb'     |
| Monte-Carlo Pivotal Bounds      |   -   |   x   |        x        |       x       | '2s', '1sl', '1su' |     'mcpb'     |
| Non-Parametric Bootstrap Bounds |   x   |   x   |        x        |       x       | '2s', '1sl', '1su' |     'npbb'     |
| Parametric Bootstrap Bounds     |   x   |   x   |        x        |       x       | '2s', '1sl', '1su' |      'pbb'     |
| Fisher Bounds                   |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'fb'      |
| Likelihood Ratio Bounds         |   x   |   -   |        x        |       x       | '2s', '1sl', '1su' |      'lrb'     |

**Important**:

- mle() and mrr() support only specific confidence bounds methods. For instance, you can't use Beta-Binomial Bounds with mle(). This will also raise an error. Use the table above to check, whether a combination of parameter estimation and confidence bounds method is supported.
- '2s': two-sided confidence bounds, '1su': upper confidence bounds, '1sl': lower confidence bounds. If Beta-Binomial Bounds are used, the lower bound represents the lower percentile bound at a specific time ((pctl) is added in the plot legend). If Fisher Bounds are used, the lower bound represents the lower time bound at a specific percentile.

### Distributions
Since version 0.1.34, dist='normal', dist='lognormal' and dist='exponential' are supported alongside the default dist='weibull'. bcm is not supported for these three (bias-correction stays Weibull-only). Confidence bounds are more limited too:

| dist          | mle() bounds               | mrr() bounds |
|---------------|-----------------------------|:------------:|
| 'weibull'     | 'fb', 'lrb'                  | see table above |
| 'normal'      | 'fb', 'lrb'                  | not supported |
| 'lognormal'   | 'fb', 'lrb'                  | not supported |
| 'exponential' | 'fb', 'chi2'                  | not supported |

'chi2' is an exact chi-square pivotal confidence interval, only available for dist='exponential' (its single-parameter model has a closed-form pivot, so likelihood-ratio bounds aren't needed there). 'lrb' is available for 'normal'/'lognormal' but not 'exponential' for the same reason, the other way round.

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

#### Modifying the Weibull plot
##### Axes labels and title
You can modify the axes label, plot title and the fontsize. Also, you can save the plot by setting save=True and path='your/own/directory/example.pdf'.
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fb',show=True, plot_title='New Project', x_label='No. of Cycles', unit='10^3', y_label='Unreliability: 1-R', xy_fontsize=12, save=True, path=r'var/user/test.pdf')
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Analysis_Plot_Modification.png){: width="500" }

##### Figure size, plot legend and median rank markers
You can customize the fontsize that is being used in the plot legend. If you don't want a legend set show_legend=False.
By default, the markers for the median ranks will be plotted. Set plot_ranks=False if you don't want median rank markers in your plot.
The figure size can be modified with fig_size=(width, height). Width and height set the figure size in inches.
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fb',show=True, show_legend=True, legend_fontsize=10, plot_ranks=False, fig_size=(7, 7))
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Analysis_Plot_Modification2.png){: width="500" }

#### Normal, LogNormal and Exponential
Set dist='normal', dist='lognormal' or dist='exponential' to fit that distribution instead of Weibull. bcm is not supported for these three; see the [Distributions](#distributions) table above for which bounds each one accepts.
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]

normal_fit = Analysis(df=failures, dist='normal', bounds='lrb', show=True)
normal_fit.mle()

lognormal_fit = Analysis(df=failures, dist='lognormal', bounds='fb', show=True)
lognormal_fit.mle()

exp_fit = Analysis(df=failures, dist='exponential', bounds='fb', show=True)
exp_fit.mle()
```

##### Exponential: exact chi-square bounds
dist='exponential' additionally supports bounds='chi2', an exact chi-square pivotal confidence interval (as opposed to the asymptotic Fisher bounds bounds='fb' also available for it).
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
exp_fit = Analysis(df=failures, dist='exponential', bounds='chi2', bounds_type='2s', show=True)
exp_fit.mle()
```

## PlotAll
PlotAll plots class objects from Analysis in one figure. Currently, only data from mle() is supported.
Theoretically, you can plot as many objects as you like -> provide a list of colors (and, for mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential(), optionally a matching list of linestyles) as a kwarg in PlotAll(objects, **kwargs).mult_weibull() / .contour_plot(). <b>
By default, predictr uses its own 6-color categorical palette. If you plot more than 6 datasets without passing your own `color`, the palette repeats, but the linestyle automatically advances (solid -> dashed -> dotted -> dash-dot) with every full pass through the palette, so up to 24 datasets stay visually distinguishable by color+shape before anything repeats outright.

**Available methods**:

| Methods        	| Description                                                           	|
|----------------	|-----------------------------------------------------------------------	|
| mult_weibull() 	| Plots multiple Analysis class instances (dist='weibull') in one Weibull plot           	|
| mult_normal() 	| Plots multiple Analysis class instances (dist='normal') in one Normal probability plot           	|
| mult_lognormal() 	| Plots multiple Analysis class instances (dist='lognormal') in one LogNormal probability plot           	|
| mult_exponential() 	| Plots multiple Analysis class instances (dist='exponential') in one Exponential probability plot (drawn on Weibull paper, since Exponential is Weibull's beta=1 special case)           	|
| contour_plot() 	| Plots contour plots when likelihood ratio bounds are used in Analysis 	|
| weibull_pdf()   | Plots one or more Weibull probability density functions. Axes are completely customizable.|
| simple_weibull()| Plots the Weibull probability plot for a given pair of beta and eta. If failures and/or suspensions are given, the median ranks are plotted as well.|
| compare()       | Fits every distribution predictr supports to one dataset and plots a probability-plot grid ranked by AIC (or Anderson-Darling), optionally with a separate PDF comparison figure.|

Note: mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential() each only accept Analysis objects that share their own dist - a ValueError is raised if you mix, e.g., a dist='weibull' object into mult_normal(). To compare fits across distributions, use compare() instead.

### Default Arguments of each method
Most of the arguments are either self explanatory or already defined in [default arguments and values](https://tvtoglu.github.io/predictr/classes/#default-arguments-and-values)

| Methods          | Default arguments                                                                                                                                                                                                                                         |
|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mult_weibull()   | x_label='Time To Failure', y_label='Unreliability', plot_title='Weibull Probability Plot', xy_fontsize=12, plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7), x_bounds=None, plot_ranks=True, save=False, color=None, linestyle=None, y_min=0.01, y_max=0.99, show=True, **kwargs    |
| mult_normal()   | x_label='Time To Failure', y_label='Unreliability', plot_title='Normal Probability Plot', xy_fontsize=12, plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7), x_bounds=None, plot_ranks=True, save=False, color=None, linestyle=None, y_min=0.01, y_max=0.99, show=True, **kwargs    |
| mult_lognormal()   | x_label='Time To Failure', y_label='Unreliability', plot_title='LogNormal Probability Plot', xy_fontsize=12, plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7), x_bounds=None, plot_ranks=True, save=False, color=None, linestyle=None, y_min=0.01, y_max=0.99, show=True, **kwargs    |
| mult_exponential()   | x_label='Time To Failure', y_label='Unreliability', plot_title='Exponential Probability Plot', xy_fontsize=12, plot_title_fontsize=14, legend_fontsize=9, fig_size=(6, 7), x_bounds=None, plot_ranks=True, save=False, color=None, linestyle=None, y_min=0.01, y_max=0.99, show=True, **kwargs    |
| contour_plot()   | show=True, style='hull', show_weibull=False, show_legend=True, color=None, x_label=r'$\widehat\beta$', y_label=None, plot_title='Contour Plot', xy_fontsize=12, plot_title_fontsize=14, legend_fontsize=9, fig_size=(6.4, 4.8), save=False, scale_mode='auto', log_ratio_threshold=10, cl_set=None, curve_fill=True, fill_alpha=0.25, **kwargs |
| weibull_pdf()    | beta=None, eta=None, linestyle=['-', '--', ':', '-.'], labels=None, x_label=None, y_label=None, xy_fontsize=12, tick_fontsize=10, legend_fontsize=9, plot_title='Weibull PDF', plot_title_fontsize=14, x_bounds=None, fig_size=None, color=None, save=False, plot_style='predictr', **kwargs |
| simple_weibull() | beta, eta, unit='-', x_label = 'Time to Failure', y_label = 'Unreliability', xy_fontsize=12, tick_fontsize=10, plot_title_fontsize=14, plot_title='Weibull Probability Plot', fig_size=(6, 7), show_legend=True, legend_fontsize=9, save=False, df=None, ds=None, **kwargs |
| compare()        | df, ds=None, bounds=None, bounds_type='2s', cl=0.9, x_label='Time to Failure', y_label='Unreliability', fig_size=(7.7, 7), y_min=0.01, y_max=0.99, plot_ranks=False, criteria='aic', plot_pdf=True, pdf_xy_fontsize=12, pdf_tick_fontsize=10, pdf_legend_fontsize=9, pdf_plot_title_fontsize=14, show=True, save=False, plot_style='predictr', **kwargs |


| Parameter(s)        | default value              | type            | description                                                                                        |
|---------------------|----------------------------|-----------------|----------------------------------------------------------------------------------------------------|
| df                  | None                       | list of floats  | List of failures                                                                                   |
| ds                  | None                       | list of floats  | List of suspensions (right-censored only)                                                          |
| plot_style          | 'predictr'                 | str             | Choose a style according to your needs. 'predictr' is predictr's own built-in style (no setup required); see matplotlib style references for other available styles. Only weibull_pdf() exposes this as its own argument - mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential()/contour_plot()/simple_weibull() inherit it from the Analysis object(s) passed in. |
| unit                | '-'                        | str             | Unit of failures and suspensions, e.g. 's', 'ms', 'no. of cycle' etc.                              |
| x_label             | depends on method          | string          | Label for the x-axis                                                                               |
| y_label             | depends on method          | string          | Label for the y-axis                                                                               |
| labels              |                            | string          | List containing the labels for the plot legend in weibull_pdf()                                    |
| xy_fontsize         | 12                         | float           | Fontsize for the axes label                                                                        |
| tick_fontsize       | 10                         | float           | Fontsize for the tick labels (the numbers on the axes). weibull_pdf() and simple_weibull() only.    |
| legend_fontsize     | 9                          | float           | Fontsize for the legend                                                                            |
| plot_title          | 'Weibull Probability Plot' | string          | Title for the plot                                                                                 |
| plot_title_fontsize | 14                         | float           | Fontsize of the plot title                                                                         |
| fig_size            | (6, 7)                     | tuple of floats | Sets figure width and height in inches: (width, height)                                            |
| save                | False                      | boolean         | If True, the plot is saved according to the path (kwargs)                                          |
| style               | 'hull'                     | string          | contour_plot() only. Defines how each dataset's confidence region is drawn: 'hull' (convex hull outline, optionally filled) or 'scatter' (raw sampled points)                       |
| show_weibull        | False                      | boolean         | contour_plot() only. If True, returns the matplotlib Figure object instead of just showing/saving it        |
| scale_mode          | 'auto'                     | string          | contour_plot() only. Scaling of the eta (y) axis: 'auto' switches to a log scale when eta spans more than log_ratio_threshold across datasets, 'linear'/'log' force the respective scale |
| log_ratio_threshold | 10                         | float           | contour_plot() only. max(eta)/min(eta) ratio above which scale_mode='auto' switches to a log scale  |
| cl_set              | None                       | list of floats  | contour_plot() only. Confidence levels to draw per dataset, e.g. [0.95, 0.9, 0.8]. If None, each object's own cl attribute is used (one curve per object, as before) |
| curve_fill          | True                       | boolean         | contour_plot() only. If True, the area enclosed by each confidence-region curve is filled          |
| fill_alpha          | 0.25                       | float           | contour_plot() only. Opacity used for curve_fill                                                    |
| plot_ranks          | True                       | boolean         | If True, median ranks will be plotted.                                                             |
| show_legend         | True                       | boolean         | If True, the legend will be plotted                                                                |
| weibull_pdf: beta, eta| None, None               | list of floats or None | Attributes from Analysis object. Pairs of beta and eta values to be plotted. Each parameter pair must have the same index value.|
| linestyle         |    ['-', '--', ':', '-.']   | list of strings      | weibull_pdf(): required, must match the length of beta/eta. mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential(): optional, must match the number of objects if given.                 |
|color        |             None               | list of strings         | List containing the colors for the plotted lines/datasets. If not given, predictr's built-in 6-color palette is used (see note above the "Available methods" table for what happens with more than 6 datasets). If given, length must match the beta/eta length (weibull_pdf()) or the number of Analysis objects (mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential(), contour_plot()).  |
| y_min, y_max        |             0.01, 0.99               | float         | mult_weibull()/mult_normal()/mult_lognormal()/mult_exponential() only. Y-axis limits (unreliability, as a fraction). Must satisfy 0 < y_min < y_max < 1.  |
| x_bounds    |                            | list of floats          | Sets x-axis boundaries: [start, stop] or [start, end, steps inbetween], respectively.|
| simple_weibull:beta, eta    |                            | float          | Weibull parameter pair which will be plotted|
| criteria            | 'aic'                      | string          | compare() only. Ranks/labels the panels by 'aic' or 'ad' (Anderson-Darling)                        |
| plot_pdf            | True                       | boolean         | compare() only. If True, also produces a separate figure overlaying every fitted distribution's PDF |
| pdf_xy_fontsize, pdf_tick_fontsize, pdf_legend_fontsize, pdf_plot_title_fontsize | 12, 10, 9, 14 | float | compare() only. Font sizes for the separate PDF figure (only used when plot_pdf=True)      |
| kwarg: path         |                            | string          | Path defines the directory and format of the figure E.g. r'var/user/.../test.pdf'                  |

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
# Set plot_ranks=True, if you want to plot the median rank markers
PlotAll(objects).mult_weibull(plot_ranks=False)
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
# Create list with custom colors and pass to the instance method
colors = ['green', 'red', 'blue']
PlotAll(objects).mult_weibull(plot_ranks=False, color=colors)
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_MLE_2s_custom_colors.png){: width="500" }

### mult_normal() / mult_lognormal() / mult_exponential()
Just like mult_weibull(), but for dist='normal'/'lognormal'/'exponential' objects respectively - each draws all given Analysis instances on that distribution's own probability paper. All objects passed to one call must share the same dist; mixing distributions raises a ValueError (use compare() if you want to compare across distributions instead).

```python
from predictr import Analysis, PlotAll

failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, dist='normal', bounds='fb', bounds_type='2s')
prototype_a.mle()

failures_b = [1.8506941739639076, 2.2685555679846954, 2.380993183650987, 2.642404955035375,
              2.777082863078587, 2.89527127055147, 2.9099992138728927, 3.1425481097241,
              3.3758727398694406, 3.8274990886889997]
prototype_b = Analysis(df=failures_b, dist='normal', bounds='fb', bounds_type='2s')
prototype_b.mle()

objects = {'proto_a': prototype_a, 'proto_b': prototype_b}
PlotAll(objects).mult_normal()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_Normal_2s.png){: width="500" }

The same objects fitted with dist='lognormal' instead, plotted with mult_lognormal():

![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_LogNormal_2s.png){: width="500" }

...and with dist='exponential', plotted with mult_exponential() (drawn on Weibull paper - see [Distributions](#distributions) for why):

![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/PlotAll_Exponential_2s.png){: width="500" }

### contour_plot()
contour_plot() only works for likelihood ratio bounds. Hence, you have to use bounds='lrb' in the Analysis class. This method supports all bounds types and all confidence levels. You can pass as many objects as you want to.

Each dataset's confidence region is drawn as a filled hull by default (curve_fill=True), with its confidence level labeled directly on the curve and its point estimate marked as a dot - see the "Default Arguments of each method" and parameter table above for style, cl_set, curve_fill, scale_mode and log_ratio_threshold.

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

#### Multiple confidence levels for one object (cl_set)
Instead of plotting each object's own cl once, pass cl_set to draw several confidence-level curves per dataset, e.g. to compare 80%, 90% and 95% confidence regions at a glance.
```python
failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
              0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
              1.3201853093496474, 1.492172669340363]
prototype_a = Analysis(df=failures_a, bounds='lrb', bounds_type='2s')
prototype_a.mle()

objects = {'initial design': prototype_a}
PlotAll(objects).contour_plot(cl_set=[0.8, 0.9, 0.95])
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Contour_plot_cl_set.png){: width="500" }

#### Normal and LogNormal fits
contour_plot() also works for dist='normal'/'lognormal' objects fitted with bounds='lrb', with axis labels matching each distribution's own parameters. Objects with different dist can't be mixed on one contour_plot() call - they don't share the same axes meaning.
```python
from predictr import Analysis, PlotAll

failures = [0.30481336314657737, 0.5793918872111126, 0.633217732127894, 0.7576700925659532,
            0.8394342818048925, 0.9118100898948334, 1.0110147142055477, 1.0180126386295232,
            1.3201853093496474, 1.492172669340363]

normal_a = Analysis(df=failures, dist='normal', bounds='lrb', bounds_type='2s')
normal_a.mle()

objects = {'sample': normal_a}
PlotAll(objects).contour_plot()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Contour_plot_normal.png){: width="500" }

Note: contour_plot() is not available for dist='exponential', since it has only one parameter (no likelihood-ratio contour to draw) and uses bounds='chi2'/'fb' instead - see [Distributions](#distributions).

### weibull_pdf()
This method plots one or more Weibull probability density functions. Axes are completely customizable.

Arguments:
weibull_pdf(self, beta=None, eta=None, linestyle=['-', '--', ':', '-.'], labels=None,
                    x_label=None, y_label=None, xy_fontsize=12, tick_fontsize=10,
                    legend_fontsize=9, plot_title='Weibull PDF', plot_title_fontsize=14,
                    x_bounds=None, fig_size=None, color=None, save=False,
                    plot_style='predictr', **kwargs)
```python
from predictr import Analysis, PlotAll

# Use analysis for the parameter estimation
failures1 = [3, 3, 3, 3, 3, 3, 4, 4, 9]
failures2 = [3, 3, 5, 6, 6, 4, 9]
failures3 = [5, 6, 6, 6, 7, 9]

a = Analysis(df=failures1, bounds='lrb', bounds_type='2s', show = False, unit= 'min')
a.mle()

b = Analysis(df=failures1, ds = failures2, bounds='fb', bounds_type='2s', show = False, unit= 'min')
b.mle()

c = Analysis(df=failures3, bounds='lrb', bcm='hrbu', bounds_type='2s', show = False, unit= 'min')
c.mle()

# Use weibull_pdf method in PlotAll to plot the Weibull pdfs
# beta contains the Weibull shape parameters, which were estimated using Analysis class. Do the same for the Weibull scale parameter eta.
# Cusomize the path directory in order to use this code
PlotAll().weibull_pdf(beta = [a.beta, b.beta, c.beta], eta = [a.eta, b.eta, c.eta],
                      linestyle=['-', '--', ':'], labels = ['A', 'B', 'C'],
                x_bounds=[0, 20, 100], plot_title = 'Comparison of three Prototypes',
                x_label='Time to Failure', y_label='Density Function',
                save=True, color=['black', 'black', 'black'], path=r'/your/custom/path/test.pdf')
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Weibull_PDF.png){: width="500" }

### simple_weibull()
This method plots the Weibull probability plot for a given pair of beta and eta. If failures and/or suspenions are given, the median ranks are plotted as well.

```python
from predictr import Analysis, PlotAll

# If save=True, you must set the path argument, e.g. path=r'/your/custom/path/test.pdf'
PlotAll().simple_weibull(beta =2.0, eta=1, show_legend=True, x_label='Cycles until failure', plot_title='Simple Weibull')

```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Simple_Weibull.png){: width="500" }

### compare()
This method fits every distribution predictr supports (Weibull, Normal, LogNormal, Exponential) to one dataset and plots a probability-plot grid, ranked by AIC (or Anderson-Darling via criteria='ad'), each on its own native paper. With plot_pdf=True (the default), it also produces a separate figure overlaying every fitted PDF on shared, linear axes.

```python
from predictr import PlotAll

failures = [93.34, 100.87, 96.41, 99.02, 108.9, 95.64, 102.31]
PlotAll().compare(df=failures, criteria='aic', plot_pdf=True)
```

| Ranked by AIC | PDF comparison |
|:---:|:---:|
| <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Compare_Normal.png" alt="PlotAll().compare() ranked by AIC" width="260"> | <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Compare_Normal_pdf.png" alt="PlotAll().compare() PDF comparison figure" width="260"> |
