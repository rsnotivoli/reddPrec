# reddPrec
**reddPrec: Reconstruction of Daily Data - Precipitation** is an R package that computes quality control for daily precipitation datasets, reconstructs the original series by estimating precipitation in missing values, creates new series in a specified pair of coordinates and creates grids.

## References
If you use the **reddPrec** package in your scientific research you can cite it as:

Serrano-Notivoli, R., De Luis, M. and Beguería, S. 2017. An R package for daily precipitation climate series reconstruction. *Environmental Modelling and Software* **89**: 190-195. [http://dx.doi.org/10.1016/j.envsoft.2016.11.005](http://dx.doi.org/10.1016/j.envsoft.2016.11.005).

A detailed description of the methodology can be found in:

Serrano-Notivoli, R., De Luis, M., Saz, M.A. and Beguería, S. 2017. Spatially-based reconstruction of daily precipitation instrumental data series. *Climate Research* **73**(3): 167-186. [https://doi.org/10.3354/cr01476](https://doi.org/10.3354/cr01476).

## Version history

**0.5** 
* Added functions for observations-estimates comparison, with validation purposes
* Improved the speed of gridding process (now in two steps: destance calculation and gridding)

## Pending tasks
* Add functions to plot the results of the QC process

## Issues and improvements
Please, address any issue, comment or suggestion in the [issues section](https://github.com/rsnotivoli/reddPrec/issues)

