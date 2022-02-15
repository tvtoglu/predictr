# Available classes
Currently, there are two classes (Analysis and PlotAll) available in the predictr package. I will continue to add new classes in near future.
## Analysis
Analysis contains all necessary methods for the Weibull analysis.  
### Default arguments and values
This table provides information on alle arguments that are passed to the Analysis class.

| Parameter           | default value              | type            | description                                                                                        |
|---------------------|----------------------------|-----------------|----------------------------------------------------------------------------------------------------|
| df                  | None                       | list of floats  | List of failures                                                                                   |
| ds                  | None                       | list of floats  | List of suspensions (right-censored only)                                                          |
| bounds              | None                       | str             | Confidence bounce method to be used in mle() or mrr()                                              |
| bounds_type         | None                       | str             | Setting for the bounds: either two-sided or one-sided                                              |
| show                | False                      | bool            | If True, the Weibull probability plot will be plotted                                              |
| bcm                 | None                       | str             | Defines the bias-correction method in mle()                                                        |
| cl                  | 0.9                        | float           | Sets the confidence level when bounds are used                                                     |
| bs_size             | 5000                       | int             | Number of bootstrap samples                                                                        |
| est_type            | 'median'                   | str             | Sets the statistic to compute from the bootstrap samples                                           |
| plot_style          | 'ggplot'                   | str             | Choose a style according to your needs. See matplotlib style references for more available styles. |
| unit                | '-'                        | str             | Unit of failures and suspensions, e.g. 's', 'ms', 'no. of cycle' etc.                              |
| x_label             | 'Time to Failure'          | string          | Label for the x-axis                                                                               |
| y_label             | 'Unreliability'            | string          | Label for the y-axis                                                                               |
| xy_fontsize         | 12                         | float           | fontsize for the axes label and ticks                                                              |
| legend_fontsize     | 9                          | float           | Fontsize for the legend                                                                            |
| plot_title          | 'Weibull Probability Plot' | string          | Title for the plot                                                                                 |
| plot_title_fontsize | 12                         | float           | Fontsize of the plot title                                                                         |
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
Analysis supports nearly all state of the art confidence bounds methods.

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
prototype_a = Analysis(df=failures, bounds='fb',show=True, show_legend=True, legend_fontsize=10, show_ranks=False, fig_size=(7, 7))
prototype_a.mle()
```
![!Backup Text](https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Analysis_Plot_Modification2.png){: width="500" }

## PlotAll
PlotAll plots class objects from Analysis in one figure. Currently, only data from mle() is supported.
Theoretically, you can plot as many objects as you like -> provide a list of colors as a kwarg in PlotAll(objects, **kwargs).mult_weibull(). <b>
For now, six colors are supported by default, but you can pass an infinit amount of colors to the mult_weibull() method.

**Available methods**:

| Methods        	| Description                                                           	|
|----------------	|-----------------------------------------------------------------------	|
| mult_weibull() 	| Plots multiple Analysis class instances in one Weibull plot           	|
| contour_plot() 	| Plots contour plots when likelihood ratio bounds are used in Analysis 	|
| weibull_pdf()   | Plots one or more Weibull probability density functions. Axes are completely customizable.|
| simple_weibull()| Plots the Weibull probability plot for a given pair of beta and eta. If failures and/or suspensions are given, the median ranks are plotted as well.|

### Default Arguments of each method
Most of the arguments are either self explanatory or already defined in [default arguments and values](https://tvtoglu.github.io/predictr/classes/#default-arguments-and-values)

| Methods          | Default arguments                                                                                                                                                                                                                                         |
|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mult_weibull()   | x_label='Time To Failure', y_label='Unreliability', plot_title='Weibull Probability Plot', xy_fontsize=12, plot_title_fontsize=12, fig_size=(6, 7), x_bounds=None, plot_ranks=True, save=False, linestyle=None, legend_fontsize=9, **kwargs                                                                |
| contour_plot()   | show_legend=True, save=False, **kwargs                  |
| weibull_pdf()    | beta=None, eta=None, linestyle=None, labels = None,x_label = None, y_label=None, xy_fontsize=10, legend_fontsize=8, plot_title='Weibull PDF', plot_title_fontsize=12, x_bounds=None, fig_size=None, color=None, save=False, plot_style='ggplot', **kwargs |
| simple_weibull() | beta, eta, unit='-', x_label = 'Time to Failure', y_label = 'Unreliability', xy_fontsize=12, plot_title_fontsize=12, plot_title='Weibull Probability Plot', fig_size=(6, 7), show_legend=True, legend_fontsize=9, save=False, df=None, ds=None, **kwargs  |


| Parameter(s)           | default value              | type            | description                                                                                        |
|---------------------|----------------------------|-----------------|----------------------------------------------------------------------------------------------------|
| df                  | None                       | list of floats  | List of failures                                                                                   |
| ds                  | None                       | list of floats  | List of suspensions (right-censored only)                                                          |
| plot_style          | 'ggplot'                   | str             | Choose a style according to your needs. See matplotlib style references for more available styles. |
| unit                | '-'                        | str             | Unit of failures and suspensions, e.g. 's', 'ms', 'no. of cycle' etc.                              |
| x_label             | depends on method          | string          | Label for the x-axis                                                                               |
| y_label             | depends on method          | string          | Label for the y-axis                                                                               |
| labels              |                            | string          | List containing the labels for the plot legend in weibull_pdf()                                    |
| xy_fontsize         | 12                         | float           | fontsize for the axes label and ticks                                                              |
| legend_fontsize     | 9                          | float           | Fontsize for the legend                                                                            |
| plot_title          | 'Weibull Probability Plot' | string          | Title for the plot                                                                                 |
| plot_title_fontsize | 12                         | float           | Fontsize of the plot title                                                                         |
| fig_size            | (6, 7)                     | tuple of floats | Sets figure width and height in inches: (width, height)                                            |
| save                | False                      | boolean         | If True, the plot is saved according to the path (kwargs)                                          |
| plot_ranks          | True                       | boolean         | If True, median ranks will be plotted.                                                             |
| show_legend         | True                       | boolean         | If True, the legend will be plotted                                                                |
| weibull_pdf: beta, eta| None, None               | list of floats or None | Attributes from Analysis object. Pairs of beta and eta values to be plotted. Each parameter pair must have the same index value.|
| linestyle         |    ['-', '--', ':', '-.']   | list of strings      | Defines the linestyle(s) in the plot. Must be greater or equal to the length of beta ans eta lists                 |
|color        |             None               | list of strings         | List containing the colormap for the plotted lines. Length of list must be equal to the beta and eta length of lists or the number of Analysis objects.  |
| x_bounds    |                            | list of floats          | Sets x-axis boundaries: [start, stop] or [start, end, steps inbetween], respectively.|
| simple_weibull:beta, eta    |                            | float          | Weibull parameter pair which will be plotted|
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

### weibull_pdf()
This method plots one or more Weibull probability density functions. Axes are completely customizable.

Arguments:
weibull_pdf(self, beta=None, eta=None, linestyle=None, labels = None,
                    x_label = None, y_label=None, xy_fontsize=10, legend_fontsize=8,
                    plot_title='Weibull PDF', plot_title_fontsize=12, x_bounds=None,
                    fig_size=None, color=None, save=False, plot_style='ggplot', kwargs)
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
