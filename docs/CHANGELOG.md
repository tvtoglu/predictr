
# Change Log - predictr
All notable changes to this project will be documented in this file.
 
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
