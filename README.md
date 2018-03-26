# reddPrec
**reddPrec: Reconstruction of Daily Data - Precipitation** is an R package that computes quality control for daily precipitation datasets, reconstructs the original series by estimating precipitation in missing values, creates new series in a specified pair of coordinates and creates grids.

## References
If you use the **reddPrec** package in your scientific research you can cite it as:

Serrano-Notivoli, R., De Luis, M. and Beguería, S. 2017. An R package for daily precipitation climate series reconstruction. *Environmental Modelling and Software* **89**: 190-195. [http://dx.doi.org/10.1016/j.envsoft.2016.11.005](http://dx.doi.org/10.1016/j.envsoft.2016.11.005).

A detailed description of the methodology can be found in:

Serrano-Notivoli, R., De Luis, M., Saz, M.A. and Beguería, S. 2017. Spatially-based reconstruction of daily precipitation instrumental data series. *Climate Research* **73**(3): 167-186. [https://doi.org/10.3354/cr01476](https://doi.org/10.3354/cr01476).

## Version history

**1.0 - [in progress]**
* Added functions for observations-estimates comparison, with validation purposes
* Improved the speed of gridding process (now in two steps: distance calculation and gridding)
* Improved the speed of quality control (QC)
* Now you can set the number of neighbours to use in QC
* New dependencies: *hydroGOF* (for scores calculation); *multiApply* (parallelization); *future* (multicore)
* Already not depending on *snowfall*

**0.4 - 2017-10**
* Improved the speed in gapFilling function
* Fixed minor issues with rounded decimals

## Pending tasks
* Add functions to plot the results of the QC process
* Add functionality to read, manage and write NetCDF files
* Add functionality to work by chunks (handle Big Data files)
* Change distance calculation to work with coordinates in degrees

## Issues and improvements
Please, address any issue, comment or suggestion in the [issues section](https://github.com/rsnotivoli/reddPrec/issues)

