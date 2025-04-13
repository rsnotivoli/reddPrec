## Introduction to reddPrec

The goal of **reddPrec** is to perform a complete reconstruction of daily precipitation. The process follows 3 steps:

1)  quality control of daily raw precipitation observations;
2)  gap filling of missing values in data series;
3)  creation of gridded datasets.

While here is provided a simple explanation about how to use the functions, a detailed explanation of the methodology can be found in the corresponding scientific publications.

To cite the original **reddPrec** package and the methodology:

Serrano-Notivoli R, de Luis M, Beguería S (2017) An R package for daily precipitation climate series reconstruction, *Environmental Modelling & Software* 89:190-195,<https://doi.org/10.1016/j.envsoft.2016.11.005>

To cite the improved and most recent version **reddPrec 3.0** and the new functions:

Huerta A, Serrano-Notivoli R, & Brönnimann S (2024) SC-PREC4SA: A serially complete daily precipitation dataset for South America. <https://doi.org/10.31223/X57D8R>

Each one of the three steps can be applied independently through the three available functions:

-   **qcPrec()**: applies several threshold-based criteria to filter original observations of daily precipitation.
-   **gapFilling()**: estimates new values for missing data in daily precipitation data series.
-   **gridPcp()**: creates a gridded precipitation dataset from a station-based dataset of observations.

While the whole package is designed to work with daily precipitation, monthly and annual data can be addressed.

## Installation

Installation is straightforward since the package is available on CRAN:

``` r
install.packages("reddPrec")
```

However, if you want to have the last stable developing version with bugs fixes and other improvements, it can be installed directly from GitHub

``` r
library(devtools)
install_github("rsnotivoli/reddPrec")
```

## Preparation of data

### Daily precipitation observations

We will use a set of daily precipitation observations from the Spanish Meterological Agency (AEMET) as example.

The package [climaemet](https://ropenspain.github.io/climaemet/) facilitates the process of downloading the data, but you will need an API key from AEMET that can be freely obtained [here](https://opendata.aemet.es/centrodedescargas/altaUsuario?).

``` r
library(climaemet)

# just replace "MY_API_KEY" by your personal API key and run this:
# aemet_api_key("MY_API_KEY", install=TRUE)

data_daily <- aemet_daily_clim(station = "all",
                               start = "2025-03-01", 
                               end = "2025-03-31",
                               return_sf = TRUE)
```

We only want to analyze the Iberian Peninsula and Balearic Islands, so we crop the spatial domain of the data with package [sf](https://r-spatial.github.io/sf/) using the boundaries of the Spanish regions (except Canary Islands) retrieved with package [giscoR](https://ropengov.github.io/giscoR/).

``` r
library(sf)
library(giscoR)
spain <- gisco_get_nuts(country = "Spain", nuts_level = "2")
spain <- spain[-which(spain$NAME_LATN=='Canarias'),]
data_daily <- st_crop(data_daily, spain)
```

The daily data must be organized in a single matrix with columns (stations) and rows (days). To do that, we use the [reshape](https://cran.r-project.org/web/packages/reshape/index.html) package which facilitates the task of *casting* data.

``` r
library(reshape)
dd <- cbind(as.data.frame(data_daily),st_coordinates(data_daily))
obs_pr <- cast(dd[,c('fecha','indicativo','prec')], fecha~indicativo)
```

AEMET codes low precipitation data (*PCP\<0.1*) as "Ip". As the function only accepts numeric values, we ensure that no strings remain.

``` r
obs_pr <- apply(obs_pr, 2, function(x){
  x <- gsub(',','.',x)
  x[x=='Ip'] <- NA
  as.numeric(x)
})
```

Lastly, we create a data.frame with the information of the stations

``` r
stations <- data.frame(ID=dd$indicativo, alt = dd$altitud, lon = dd$X, lat = dd$Y)
stations <- stations[-which(duplicated(stations$ID)),]
```

### Creation of geospatial (raster) data

Te estimation of precipitation is used in all stages of the reconstruction process, and it uses environmental data as predictors. As we don't have available this information for each station's location, it can be extracted from raster data, but must create it first.

In this case, we will use seven predictors: 1) elevation, 2) latitude, 3) longitude, 4) the first four Pincipal Components of a collection of topographic variables.

First, we use the [elevatr](https://github.com/jhollist/elevatr) package to derive a raster of elevations. Then, the [terra](https://rspatial.org/pkg/index.html) package will allow for calculating the required environmental predictors derived from elevations.

``` r
# topographic variables

dem <- get_elev_raster(data_daily, z = 5)
dem <- rast(dem)
dem <- crop(dem, spain)
dem <- mask(dem, spain)
dem[dem<0] <- 0

asp <- terra::terrain(dem, v ='aspect')
aspectcosine <- cos(asp)
aspectsine <- sin(asp)
dist2coast <- costDist(dem)
lon <- rast(cbind(crds(dem),crds(dem)[,1]),type='xyz',crs='EPSG:4326')
lat <- rast(cbind(crds(dem),crds(dem)[,2]),type='xyz',crs='EPSG:4326')
roughness <- terra::terrain(dem, v ='roughness')
slope <- terra::terrain(dem, v ='slope')
tpi <- terra::terrain(dem, v ='TPI')
tri <- terra::terrain(dem, v ='TRI')

covs <- c(aspectcosine, aspectsine, dist2coast, dem, lon, lat, roughness, slope, tpi, tri)
names(covs) <- c('aspectcosine', 'aspectsine', 'dist2coast', 
                         'elevation', 'longitude', 'latitude', 'roughness', 
                         'slope', 'tpi', 'tri')

# pca of topo variables
pca_cvs <- terra::princomp(covs[[-c(4:6)]])
pca_cvs <- terra::predict(covs[[-c(4:6)]], pca_cvs, index=1:4)
names(pca_cvs) <- c("pc1", "pc2", "pc3", "pc4")
plot(pca_cvs)
```

Lastly, we extract the values of de PCs to the stations

``` r
stations <- vect(stations, geom=c('lon','lat'),crs = 'EPSG:4326',keepgeom=TRUE)
e <- terra::extract(pca_cvs, stations, ID=F)
stations <- cbind(stations, e)
stations <- as.data.frame(stations)
```

Those stations with no overlapping with raster don't have data, so we remove them

``` r
stations <- stations[complete.cases(stations),]
obs_pr <- obs_pr[, match(stations$ID, colnames(obs_pr))]
```

At his point, we have available the dataset with the original data.

``` r
st <- vect(stations, geom = c('lon', 'lat'), crs = 'EPSG:4326', keepgeom = TRUE)
st$ndata <- colSums(!is.na(obs_pr)) * 100 / nrow(obs_pr)
dem_df <- as.data.frame(dem, xy = TRUE)
colnames(dem_df)[3] <- "elevation"
st_sf <- st_as_sf(st)
st_sf$prec_acum <- colSums(obs_pr, na.rm = TRUE)

ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradient(name = "Elevation (m)", low = "lightgreen", high = "darkgreen") +
  geom_sf(data = st_transform(st_sf, crs = st_crs(spain)), aes(size = prec_acum), color = "black", alpha = 0.7) +
  scale_size_continuous(name = "Cumulated PCP (mm)") +
  geom_sf(data = spain, fill = NA, color = "black") +
  coord_sf(crs = st_crs(spain)) +
  theme_minimal() +
  theme(legend.position = "right")
```

And the environmental pedictors

``` r
env <- c(dem, lon, lat, pca_cvs)
names(env)[1:3] <- c("alt","lon","lat")
plot(env)
```

## Quality control

The quality control (QC) function allows for a customization of the thresholds applied to each criteria. In this case, we will flag and remove the observations based on the following conditionals:

-   the coordinates (lon, lat) will be used as predictors
-   the 15 nearest observations to estimate precipitation
-   no maximum radius of searching nearest observations
-   a type of model must be chosen from: *learner_glm* (Generalized Linear Model), *learner_rf* (Random Forest), *learner_svm* (Support Vector Machines), *learner_xgboost* (eXtreme Gradient Boosting)
-   we will apply the 5 reference QC criteria: suspect value, suspect zero, suspect outlier, suspect wet day and suspect dry day.
-   a threshold of 10 times higher or lower observation than estimate to detect outliers
-   a threshold of 0.99 wet probability and an observed magnitude higher than 5 mm to detect suspect zeros.
-   a threshold of 0.01 wet probability and an observed magnitude lower than 0.1 mm to detect suspect values higher than 5 mm.
-   two processor cores to compute results in parallel (usually the higher the faster, depending on your hardware capabilities)

``` r
library(reddPrec)
qcdata <- qcPrec(prec = obs_pr, 
                 sts = stations, 
                 model_fun = learner_glm,
                 crs = 'EPSG:4326', 
                 coords = c('lon','lat'),
                 coords_as_preds = TRUE, 
                 neibs = 15, 
                 thres = NA,
                 qc = 'all', 
                 qc3 = 10, 
                 qc4 = c(0.99, 5), 
                 qc5 = c(0.01, 0.1, 5),
                 ncpu = 8)
```

Depending on the velocity of your processor(s), the job will be done quick or slow, but this particular task (with our setting) should take a few minutes. If you use a larger dataset of observations, the computing time rises.

The result is a list of two elements:

-   **cleaned**: a matrix (stations x days) with the filtered observations (cleaned data).
-   **codes**: a matrix (stations x days) with the codes corresponding to the reasons of values removal
    -   "1" (suspect value): obs==0 & all(neibs\>0)
    -   "2" (suspect zero): obs\>0 & all(neibs==0)
    -   "3" (suspect outlier): obs is "qc3" times higher or lower than the estimate
    -   "4" (suspect wet): obs==0 & wet probability \> "qc4[1]" & estimate \> "qc4[2]"
    -   "5" (suspect dry): obs\>"qc5[3]" & dry probability \< "qc5[1]" & estimate \< "qc5[2]")

In our example, the outliers (QC3) were the most flagged values (3.74%), followed by suspect zeros (QC2, 0.69%), suspect values (QC1, 0.33%), suspect wet (QC4, 0.16%) and suspect dry (QC5, 0.11%).

``` r
allcodes <- as.numeric(as.matrix(qcdata$codes))
flagged <- round(table(allcodes)*100/length(allcodes),2)
flagged
```

```         
## allcodes
##    1    2    3    4    5 
## 0.33 0.69 3.74 0.16 0.11 
```

## Gap filling

Missing values are common in raw original data series of observations. After the QC process, the number of these missing data is increased, and the resulting gaps can affect to further analyses at coarser scales (seasonal, annual, etc.). To solve that, a large collection of infilling methods exist, most of them based on regression algorithms.

The gap filling process in **reddPrec** uses the nearest observations (a number defined by the user) and their associated environmental information (which we added to the *stations* data.frame in the "data preparation" section) to compute a Reference Value (RV). The RV is computed through a multivariate linear regression that uses all those data.

The function returns a data.frame with different estimates for all days of every station:

-   **wd_pred**: probability (0 to 1) of wet day
-   **raw_pred**: rainfall prediction (in the original units) obtaind from the model
-   **mod_pred**: modified prediction based on wet day probability (if wd_pred\<0.5 -\> mod_pred = 0)
-   **st_pred**: standardized prediction (based on mod_pred) from selected standardization method (ratio or quantile)
-   **err**: standard error of the model (in the original units)
-   **neibs**: number of nearest stations (neighbors) used in the model. They can vary if a searching threshold (thres) is set.

``` r
gf_res <- gapFilling(prec = qcdata$cleaned, 
                     sts = stations,
                     model_fun = learner_glm,
                     dates = seq.Date(as.Date('2025-03-01'), 
                                      as.Date('2025-03-31'),
                                      by ='day'), 
                     stmethod = 'ratio', 
                     ncpu = 8, 
                     thres = NA, 
                     neibs = 15,
                     coords = c('lon','lat'),
                     crs = 'EPSG:4326',
                     coords_as_preds = TRUE,
                     window = 11)
```

Our example shows slight differences between **mod_pred** and **st_pred** since we used a window of 11 days, but the length of the dataset is just 31 days.

We can compare, for example, the original data series and their reconstructions.

``` r
library(ggplot2)
library(patchwork)
library(scales)
library(hexbin)

# Daily and monthly data
pearson_daily <- round(cor(gf_res$obs, gf_res$st_pred, use = "pairwise.complete.obs"), 2)
obs <- gf_res[complete.cases(gf_res),]
o <- aggregate(obs$obs, by = list(obs$ID), FUN = sum)
p <- aggregate(obs$st_pred, by = list(obs$ID), FUN = sum)
pearson_monthly <- round(cor(o[,2], p[,2], use = "pairwise.complete.obs"), 2)

daily_df <- gf_res
monthly_df <- data.frame(obs_sum = o[,2], pred_sum = p[,2])

# Axes & legend
all_values <- c(daily_df$obs, daily_df$st_pred, monthly_df$obs_sum, monthly_df$pred_sum)
axis_min <- floor(min(all_values, na.rm = TRUE))
axis_max <- ceiling(max(all_values, na.rm = TRUE))
hex_daily <- hexbin(daily_df$obs, daily_df$st_pred, xbins = 60)
hex_monthly <- hexbin(monthly_df$obs_sum, monthly_df$pred_sum, xbins = 60)
max_count <- max(c(hex_daily@count, hex_monthly@count))

# 
p1 <- ggplot(daily_df, aes(x = obs, y = st_pred)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log",
    limits = c(1, max_count),
    name = "Frequency",
    labels = label_number(accuracy = 1)
  ) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  annotate("text", x = axis_min + 5, y = axis_max - 10,
           label = paste0("Pearson = ", pearson_daily), 
           hjust = 0, color = "black") +
  labs(
    x = "Observed (mm)",
    y = "Predicted (mm)",
    title = "Daily PCP comparison (all stations)"
  ) +
  coord_fixed() +
  xlim(axis_min, axis_max) +
  ylim(axis_min, axis_max) +
  theme_minimal()

# 
p2 <- ggplot(monthly_df, aes(x = obs_sum, y = pred_sum)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log",
    limits = c(1, max_count),
    name = "Frequency",
    labels = label_number(accuracy = 1) 
  ) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  annotate("text", x = axis_min + 5, y = axis_max - 10,
           label = paste0("Pearson = ", pearson_monthly), hjust = 0, color = "black") +
  labs(
    x = "Observed (mm)",
    y = "Predicted (mm)",
    title = "Monthly PCP comparison (all stations)"
  ) +
  coord_fixed() +
  xlim(axis_min, axis_max) +
  ylim(axis_min, axis_max) +
  theme_minimal()

# 
p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "right")
```

## Gridding

The final step creates a gridded product based on the environmental variables in raster format. Different inputs can be used as observations: (i) the original observations or (ii) the reconstructed series. Reconstructed series are recommended due to the neighboring stations for all pixels will remain the same in all days of the period. Otherwise, some inconsistencies (inhomogeneities) could be imputed to the grid and spatially propagated.

We will use reconstructed series for our example, 15 neighbors with no radius limitation. Please consider that, although here the three steps (QC, gap filling and gridding) are presented as a workflow, they can be run separately, meaning that the options used in gridding can be different than gap filling, for example.

The gridding process will take a long time depending on multiple factors: the grid resolution, the number of days, the number of neighbors, and the number of used CPUs, mainly. In our example, we will aggregate the grid at a coarser resolution and will reduce the time period to two days just to reduce the computing time. (This example takes about 3 minutes each day)

A proper gridding of value should use the original observations when available and estimates for those days with missing data. You decide what estimates you will use from the output of gap fillling process. In this case, we will use the *st_pred* estimate (the one with standardization). Be careful with this choice because sometimes, depending on the available data and other climatic factors, the standardization process could create unrealistic estimates.

As gridding is a separate process from the previous steps, we can choose different covariables as predictors. In this case, we will only use *alt*, *lat* and *lon* since they yield more coherent spatial patterns. (Note that we select only the first 4 columns from stations [ID, lat, lon, lat]).

``` r
recs <- gf_res$obs
recs[is.na(recs)] <- gf_res$st_pred[is.na(recs)]
rec <- data.frame(date = gf_res$date, ID = gf_res$ID, pred = recs)
rec <- cast(rec, date~ID)
rec <- rec[,-1]

day21 <- as.numeric(rec[21,])

gridPcp(prec = day21,
        grid = env,
        dyncovars = NULL,
        sts = stations[,1:4],
        model_fun = learner_glm,
        dates = as.Date('2025-03-21'),
        ncpu = 8,
        thres = NA,
        neibs = 15,
        coords = c('lon','lat'),
        crs = 'EPSG:4326',
        coords_as_preds = TRUE,
        dir_name = "grid_test")
```

The function creates 2 folders in the working directory, one containing the daily grids of precipitation estimates and one containing the daily grids of uncertainties (errors of the model)

``` r
pre <- terra::rast('./pred_grid_test/20250321.tif')
title <- "Daily Precipitation 2025-03-21"
vals <- c(-1, 0, 1, 2, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 1000)
  legnd <- paste(vals[1:21],vals[2:22],sep='-')
  legnd[1] <- '0'
  legnd[21] <- '>110'
  
  cols <- c('#ffffff', '#ffd576', '#ffe685', '#fff0b0',
               '#fffecf', '#b2fdad', '#a2eda6', '#61cf87',
               '#00c56a','#b2f1fe', '#72f0fe', '#1cdefe',
               '#31c3fe','#6e90fe','#b77dff','#ca8eff',
               '#dc9ae7','#f4a1f3','#fec1ff','#ffdcfe', '#e7e7e7')
  
  m <- cbind(vals[1:21],vals[2:22],1:21)
  d <- classify(pre, m)
  fadd <- function() plot(vect(spain),add=T)
  plot(d,  type="interval", breaks = 1:22, col = cols, plg=list(legend=legnd),
       main = title, fun = fadd)
```
