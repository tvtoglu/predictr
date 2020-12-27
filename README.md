# predictr - Predict the Reliability

predictr: predict + reliability, in other words: A tool to predict the reliability.<br>
The aim of this package is to provide state of the art tools for all kinds of Weibull analyses. <br>
predictr already includes many methods (see tables below). A guideline when to use which method will be added soon.

![](https://img.shields.io/pypi/v/predictr?color=blue&style=flat&label=pypi)
[![Downloads](https://pepy.tech/badge/predictr)](https://pepy.tech/project/predictr)
![](https://img.shields.io/pypi/pyversions/predictr)
![](https://img.shields.io/pypi/l/predictr)
![](https://img.shields.io/github/stars/tvtoglu/predictr?style=flat)

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install predictr.

```bash
pip install predictr
```

## Usage
### Import predictr in python
```python
from predictr import Analysis
```
### How to use the Analysis class after importing
If both df and ds are None, an error will be raised. Configure all needed arguments that fits your use case.<br>
Arguments to pass in Analysis class:<br>
Analysis(df: list = None, ds: list = None, show: bool = False, plot_style='ggplot', bounds=None, bounds_type='2s',
                 cl=0.9, bcm=None, bs_size=5000, est_type='median', unit='-'))<br>

**Step 1: Create an instance of the class Analysis** <br>
The class object is some prototype data, so we will call it "prototype_test". You can name it whatever you want.<br>
prototype_test = Analysis(df = [200, 300, 400], show=True)<br>

**Step 2: Use one of the instance methods (mle() or mrr()) of class Analysis**<br>
prototype_test.mle() -> MLE or prototype_test.mrr() -> MRR<br>

### Default Parameter values
df: list = None -> failures in seconds, days, no. of cycles etc., e.g. df = [100, 120, 80, 300]<br>
ds: list = None -> suspensions (right-censored) in seconds, days, no. of cycles etc., e.g. ds = [300, 400, 400]<br>
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
cl=0.9 -> configure the confidence level in the intervall (0, 1.0)<br>
bcm=None -> Define the bias-correction method when the MLE is being used. Bootstrap bias-corrections are dependent on the number of bootstrap replication and the chosen statistic, e.g. if bcm = 'np_bs': bs_size = 5000 and est_type = 'median'.<br>
bs_size = 5000 -> Resampling/Bootstrap sample size (number of replication). bs_size should be greater than or equal to 2000 for accurate results. The higher the number of replication, the longer it takes to compute the bias-correction.<br>
est_type = 'median' -> When using bootstrap bias-corrections, this argument decides which statistic to compute from the bootstrap samples.<br>
The following table provides possible configurations. Bias-corrections for mrr() are not supported, yet.<br>

| Bias-correction method              | mle() | mrr() | argument value | config. |             statistic            |
|-------------------------------------|:-----:|:-----:|:--------------:|:-------:|:--------------------------------:|
| C4                                  |   x   |   -   |      'c4'      |    -    |                 -                |
| hrbu                                |   x   |   -   |     'hrbu'     |    -    |                 -                |
| non-parametric Bootstrap correction |   x   |   -   |     'np_bs'    | bs_size | 'mean', 'median', 'trimmed_mean' |
| Parametric Bootstrap correction     |   x   |   -   |     'p_bs'     | bs_size | 'mean', 'median', 'trimmed_mean' |

unit = '-' -> Unit of the elements in df and ds, e.g. unit = 'seconds', unit = 'days', unit = 'ms' etc.

### How to use the Maximum Likelihood Estimation (MLE)
object = Analysis()
object.**.mle()**

#### Uncensored sample
Example: 
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='fisher',show=True)
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_Fisher_uncensored.png" height="700" />

#### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='lrb',show=True)
prototype_a.mle()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MLE_LRB_censored.png" height="700" />

### How to use the Median Rank Regression (MRR)
object = Analysis()
object**.mrr()**

#### Uncensored sample
Example: 
```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
prototype_a = Analysis(df=failures, bounds='bbb',show=True)
prototype_a.mrr()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_BBB_uncensored.png" height="700" />

#### Censored sample
Example: 

```python
failures = [0.4508831,  0.68564703, 0.76826143, 0.88231395, 1.48287253, 1.62876357]
suspensions = [1.9, 2.0, 2.0]
prototype_a = Analysis(df=failures, ds=suspensions, bounds='mcpb',show=True)
prototype_a.mrr()
```
<img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/MRR_MCPB_censored.png" height="700" />

## To Do
I will add a homepage with more detailed examples and guidelines for non-experts in the field of reliability engineering.
