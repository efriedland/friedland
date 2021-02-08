Making Functions Do Your Econometrics For You
============================================
Currently being updated to use methods and custom classes. 

Overview
--------

friedland contains the work flow tools I use regularly for time series econometrics. A repository of the functions and datasets I use. Some inspiration from code from Professor Roger Keener at Illinois (found here http://www.econ.uiuc.edu/~econ508/e-ta.html) while the work flow methodologies are attributable to my classes with Professor Robert Kaufmann at Boston University. 

Installation
------------

``` r
# The easiest way to get the development version from GitHub is to use remotes:
install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true") # if the loading fails due to outdated dependencies
remotes::install_github("efriedland/friedland")
```

Getting Started
---------------

``` r
# Loading the package and data
library(friedland)
# "MacKinnon" is the workhouse table for critical values used to assess significance of unit root tests
# "Asymmetry" is an example dataset containing price of gasoline, crude oil, utilization and stocks 
# "NhemShem" is an example dataset containing temperature and radiative forcing data across
# the Northern and Southern Hemispheres
dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12) # converted to a time series object

# identify if variables are stationary or nonstationary
UnitRoot(dat.ts[,"Pgasoline"], drift = T, trend = T) # augmented dickey-fuller testing a single variable
UnitRoot(dat.ts, drift = T, trend = T) # vectorized to accepted the entire dataset

# test for cointegration
mu <- residuals(dynlm(Pgasoline ~ Pcrude, dat.ts))
UnitRoot(mu, drift = T, trend = T, cointvariables = 2)

# build a Dynamic Ordinary Least Squares model 
buildDOLS(Pgasoline ~ Pcrude, dat.ts, robusterrors = T)

# build a Vector Error Correction model
buildVECM(Pgasoline ~ Pcrude, dat.ts, stationary_vars = ~ Stocks, robusterrors = T, SplitError = T)

# speed up the process with the lazier "lazy" functions
# test for cointegration across all possibilities of independent variable combinations
results <- lazyCoint("Pgasoline", dat.ts)
results

# evaluate the best by comparing the HAC estimated significance levels of all DOLS coefficients
lazybuildDOLS(results, dat.ts)
```

Getting help
------------

Feel free to email me at evan.friedland@gmail.com.
