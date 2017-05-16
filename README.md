Making Functions Do Your Econometrics For You

<!-- README.md is generated from README.Rmd. Please edit that file -->
============================================

Overview
--------

friedland contains the work flow tools I use regularly for time series econometrics. A repository of the functions and datasets I use. Some inspiration from code from Professor Roger Keener at Illinois (found here http://www.econ.uiuc.edu/~econ508/e-ta.html) while the work flow methodologies are attributable to my classes with Professor Robert Kaufmann at Boston University. 
It is currently being updated to use methods and custom classes. 

Installation
------------

``` r
# The easiest way to get the development version from GitHub is to use remotes:
install.packages("remotes")
remotes::install_github("efriedland/friedland")
```

Usage
-----

``` r
library(friedland)
data(MacKinnon) # the workhouse table for critical values used to assess significance of unit root tests
data(Asymmetry) # example dataset containing price of gasoline, crude oil, utilization and stocks 
dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)

UnitRoot(dat.ts[,"Pgasoline"], drift = T, trend = T) # augmented dickey-fuller testing single variables
UnitRoot(dat.ts, drift = T, trend = T) # the entire dataset
```

Getting help
------------

Feel free to email me at evan.friedland@gmail.com.
