
# Change Log - predictr
All notable changes to this project will be documented in this file.

## Planned updates
### Added
 - .summary() method
 - Weibull AFT

## [0.1.34] - 2026-08-18
### Added
 - Added support for following distributions: **normal**, **lognormal**, and **exponential**
   - predictr had the arg dist='weibull' since the first release. Now, I have updated the classes with dist='normal', dist='lognormal' and dist='exponential'. Check the documentation for use cases
 - New class method: PlotAll().compare() for goodness-of-fit tests for a given dataset:
   - plots pdfs and cdf
   - AIC and AD criterion available
   - optional: plot median ranks and confidence bounds for each distribution

   | Ranked by AIC | PDF comparison |
   |:---:|:---:|
   | <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Compare_Normal.png" alt="PlotAll().compare() ranked by AIC for a random normal sample of 7" width="260"> | <img src="https://raw.githubusercontent.com/tvtoglu/predictr/main/docs/images/Compare_Normal_pdf.png" alt="PlotAll().compare() PDF comparison figure" width="260"> |
 - contour_plot() now supports dist='normal'/'lognormal' objects fitted with bounds='lrb', with axis labels matching each distribution's own parameters ($\widehat\mu$/$\widehat\sigma$ for Normal, $\widehat\mu_{\ln(t)}$/$\widehat\sigma_{\ln(t)}$ for LogNormal). Objects with different dist can no longer be plotted together on one contour_plot() (they'd share axes that don't mean the same thing), and a clear error is raised if an object's LRB bounds weren't computed yet

### Improved
There are many more smart quirks in predictr now, like for instance (not gonna name them all her, check the documentary! :D):
- contour_plot()'s scale_mode ('auto', 'linear', 'log') - previously Weibull-only - now works for Normal and LogNormal fits too, automatically switching to a log scale whenever a dataset's confidence region is spread very widely, so it stays clearly readable instead of being squashed flat.
- y-axis lower limit: new standard is 0.01 (1 %). This lower value is more practical in everyday use than 0.001. However, the limits are still customizable


## [0.1.33] - 2025-08-13
### Improved
 - Say hello to the brand-new 'predictr' style, now replacing 'ggplot' as the default across all classes! Update to the latest version and see for yourself, or check out the docs for a sneak peek
 - Unified categorical color palette used across mult_weibull(), contour_plot() and weibull_pdf() (previously three different, uncoordinated color sources). For more than 6 datasets, the palette now cycles color together with linestyle (solid/dashed/dotted/dash-dot), so datasets stay visually distinguishable instead of silently repeating a color
### Added
 - contour_plot() no support the args: curve_fill, scale_mode and cl_set  
    - curve_fill: if True, the contour plot will be filled
    - scale_mode: standard: 'auto' ('linear', 'log') -> if a linear scale would not properly show the curves due to big differences in x and y values of the datasets, eta will be show as log(eta).  
    - cl_set: is a list containing confidence levels to be drawn for each dataset
 - The confidence level each object is now shown as a contour line label (clabel). Hence, the legend is now decluttered  
 - Point estimates of each dataset are shown as a scatter point
 - New built-in matplotlib style 'predictr' (the new default plot_style), with tuned defaults for presentation/documentation use: colors, tick padding/size, title spacing, and a new tick_fontsize argument to size axis-label and tick-number text independently
### Fixed
 - MLE/MRR plots no longer clip the legend at the right figure edge when its position overflows the canvas (widened automatically instead)
 - contour_plot(save=True) no longer crashes with an AttributeError
 - weibull_pdf() no longer reuses a previous plot's leftover axes (log scale, probability-paper formatting) when no fig_size is given, and its line widths no longer default to an overly thick 4pt
 - The plot border is no longer thicker and lighter than the grid lines, which could make the outermost gridline invisible

## [0.1.32] - 2025-08-12
### Fixed
 - Fixed import handling for get_cmap in order to ensure compability throughout different matplotlib version

## [0.1.31] - 2025-08-11
### Improved
 Long time no see! I’m back! Expect many exciting updates as we make Predictr the best reliability tool out there!  
 - Vectorized MCPB for MRR: up to 300x faster  
 - Vectorized parametric and non-parametric bootstrap bias correction: up to 900x faster  
 - Vectorized bootstrap bounds: up to x900 faster  
 - Vectorized LRB: up to 200x faster  
### Fixed
 - PlotAll.contour_plot() now searches for beta solutions for a fixed eta value, hence, adding more solutions to the contour plot. The curve now looks smoother.

## [0.1.30] - 2025-05-22
### Added
 - Added new arguments for PlotAll.contour_plot(): x_label=r'$\widehat\beta$', y_label=r'$\widehat\eta$', plot_title='Contour Plot', xy_fontsize=12, plot_title_fontsize=12, legend_fontsize=9, fig_size=(6.4, 4.8). You can now fully customize this plot.

### Fixed
 - Showing the legend in PlotAll.contour_plot() now works again. 

## [0.1.29] - 2025-05-21
### Added
 - Added two new styles to the contour plot in PlotAll.contour_plot(): 'spline' and 'angular_line. These are additional styles to the scatter plot. New default is 'spline' for smoother curves through the data points.

## [0.1.28] - 2025-03-08
### Fixed
 - '1sl' attribute inside the plot is now being displayed correctly. Bug was found by bobbolous! Thank you.


## [0.1.27] - 2023-03-01
### Fixed
 - MRR based NPBB are now computed correctly. A typo in the code lead to a bug that computed solely MLE based bounds, no matter if 'mrr' was set as an attribute. A big THANK YOU to William Gandler who reported this bug and various spelling mistakes!
 - Minor other spelling mistakes pointed out by Antonia, Sarahi and Edna. ¡Viva México, Cabrones!

## [0.1.26] - 2022-03-28
### Added
 - added new style argument plot in PlotAll.contour_plot(style='...'): Currently scatter and angular line plots are supported

### Fixed
 - Boosted LRB computation by fully vectorizing each computation step and using algeabric tricks (as a tip: don't use np.power(), because it is way too slow for big computations). Now, it takes half the amount of time to have more precise results. 


## [0.1.25] - 2022-02-15
### Added
 - added show argument in PlotAll.contour_plot()
 - added customizable x-axis limits in PlotAll.mult_weibull: x_bounds=[start, end]


## [0.1.24] - 2022-01-15
 
### Added
 - In previous versions of predictr, saving Weibull plots required to plot and show them. Now, saving plots and actually showing them are independent from each other and can be called separately.
 
## [0.1.23] - 2022-01-02
 
### Fixed
 - Computation of one-sided Likelihood Ratio Bounds had an error, which always resulted in a 95% one-sided bound; no matter what was set by the user.

### Added
 - Legend fontsize can now be set in mult_weibull(legend_fontsize=9)
 - Bootstrap bounds now fully support censored data
 - Non-parametric boostrap bias-correction now fully supports censored data (parametric already did)

## [0.1.22] - 2021-11-02
 
### Fixed
 - Units in axis labels are now according to the International System of Units (SI):e.g. '[%]' -> 'in %'


## [0.1.21] - 2021-08-23
 
### Fixed
 - Custom plot title was not updating in PlotAll
 - MRR plots: pvalue was not shown when adj. ranks were being used
 - RBA now works correctly in PlotAll
 
### Added
 - Linestyles can now be customized  in PlotAll -> mult_weibull(linestyle=['-', '-.']) etc.


## [0.1.20] - 2021-04-15
 
### Fixed
 - Fixed bug that prevented custom axis labels and titles in mult_weibull()

## [0.1.19] - 2021-03-23
 
### Changed
 - Removed set_cmaps argument in PlotAll. You can now customize the colors for each PlotAll method individually using the color argument.
 - Improved consistency of method arguments throughout code: fontsizes, nomenclature etc.

### Added
 - PlotAll methods can plot median rank markers
 - Possibility to save plots for all methods

### Fixed
 - Fixed bug that resulted in empty saved plots

## [0.1.18] - 2021-03-17
 
### Fixed
 - fig_size didn't work in mult_weibull() due to a bug in the code

## [0.1.17] - 2021-03-13
 
### Added
 - Better documentation in the code

### Change
 - Minor changes to have more consistent nomenclature in the code

## [0.1.16] - 2021-03-03
 
### Fixed
 - fixed last bug in the median rank computation. Bug occured when multiple failure times were identical

### Added
 - PlotAll().simple_weibull(): Plots Weibull probability plot according to input (Weibull parameters) without calculations
 - PlotAll().weibull_pdf: Plots one or more Weibull probability density functions with fully customizable figure attributes (color, labels, width, height, fontsize etc.)
 - Added ability to customize axis labels, title, fontsize etc for Analysis() and PlotAll
 - Added ability to hide legend in Weibull plot

## [0.1.15] - 2021-02-27
 
### Fixed
 - fixed wrong percentile values for adjusted ranks when suspenions have lower time to failure values than actual failures
 - permanent fix for the disappearing Weibull probability line (see changelog for version 0.1.14)

### Added
 - raise ValueError when np_bs and p_bs bias-correction methods are applied to data that has suspensions

## [0.1.14] - 2021-02-16
 
### Fixed
 - temporary fix for instances in Analysis where the Weibull probability line disappears when two-sided bounds are used
 
## [0.1.13] - 2021-01-23
 
### Fixed
 - hrbu was misspelled in the code and raised an error when calling this bias-correction method

## [0.1.12] - 2021-01-14
 
### Added
 - contour_plot() method in PlotAll class: Plots likelihood ratio contours for Analysis class instances

## [0.1.11] - 2021-01-01
 
### Fixed
 - When using npbb and pbb as bounds, bounds types 1sl and 1su would return the same bounds limits. Fixed the percentile method.

## [0.1.10] - 2021-01-01
 
### Fixed
 - Minor fix and code restructuring

## [0.1.9] - 2021-01-01
 
### Added
 - new GithubPage: https://tvtoglu.github.io/predictr/
 - new Class: PlotAll -> plot multiple Weibull plots in one figure

## [0.1.7] - 2020-12-29
  
### Changed
 - changed get_bx_percentile() to a classmethod

## [0.1.6] - 2020-12-29
 
### Added
- Official changelog
- New static method get_bx_percentile() for Analysis class. You can now get the time values for given BX-lives

### Changed
  
- Argument "bounds" for Fisher bounds: 'fisher' -> 'fb'. This is more in line with the other confidence bounds
