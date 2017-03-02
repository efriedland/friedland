\name{friedland-package}
\alias{friedland-package}
\alias{friedland}
\docType{package}
\title{
\packageTitle{friedland}
}
\description{
\packageDescription{friedland}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{friedland}
\packageIndices{friedland}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{friedland}

Maintainer: \packageMaintainer{friedland}
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{ \code{\link{adf}} \code{\link{UnitRoot}} \code{\link{Asymmetry}} \code{\link{NhemShem}} \code{\link{MacKinnon}} \code{\link{aic}} \code{\link{gatelyPriceAsymm}} }
\Walkthrough{
  \code{
  ## remove.packages("friedland") # may be needed for reinstalling, may have to reopen Rstudio

  ## load package ##
  # install.packages("devtools")
  # devtools::install_github("efriedland/friedland")
  # library(friedland)
  # ?friedland

  ## load data ##
  data(Asymmetry) # dataset
  data(MacKinnon) # needed for computing pvalue of the Augmented Dickey-Fuller statistic (adf)
  # package also includes data(NhemShem)

  ## turn data into time series class ##
  dat.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)

  ## choose max lag with ad hoc method ##
  floor(nrow(dat.ts)^(1/3))

  ## determine if variable stationary or nonstationary, where k is maxlag ##
  args(UnitRoot) # default
  UnitRoot(variable = dat.ts[,"Pcrude"], k = 5, trend = F)

  ## do the above process quickly across all variables in the dataset with the UnitRootApply() function ##
  UnitRootApply(dat.ts, k = 5, trend = F) # drift = T and cointvariables = 1 is the default

  ## We can decompose the Pcrude value into its vectors of its 
  ## maximum values, cumulative upward movements, and cumulative 
  ## downward movements with the gatelyPriceAsymm() function ##
  str(gatelyPriceAsymm(dat.ts[,"Pcrude"])) # output is a list class
  ## putting all the data in the same table despite not having the same start and end dates ##
  newdat.ts.list <- c(as.list(dat.ts), gatelyPriceAsymm(dat.ts[,"Pcrude"]) )
  newdat.ts <- do.call(cbind, newdat.ts.list) # will keep NA
  # do.call(ts.intersect,newdat.ts.list) ## would not keep NA row, don't worry about this

  ## double check unit roots now that we have new variables ##
  UnitRootApply(newdat.ts, k = 5, trend = F)

  ## Test for cointegration ##
  # We do a simple regression for the residual where we use the decomposed prices of crude instead
  mu <- residuals(dynlm(Pgasoline ~ Pmax + Prec + Pcut, data = newdat.ts)) # notices uses dynlm() to keep residuals in ts format, not lm()
  class(mu)
  UnitRoot(mu, k=5, trend = F, cointvariables = 4)$reason
  # the residual is stationary, there is cointegration


  # Use the buildDOLS function which does a number of steps
  ?buildDOLS
  DOLS.list <- buildDOLS(Pgasoline ~ Pmax + Prec + Pcut, data = newdat.ts)

  DOLS.list$data.names # returns the column names of the dataset you put in 

  DOLS.list$call # returns the dynlm formula used to make the DOLS

  DOLS.list$selection # returns a table showing the lags/leads and the SBC, fixed observations, for choosing the number of lags/leads to use

  DOLS.list$leadslags # returns the # of lags/leads identified for the model with the smallest SBC

  fit <- DOLS.list$model # returns the final model, where you can also find the dates used

  summary(fit)
  start(fit)
  end(fit)

  DOLS.list$robusterrors # returns the robust errors using Newey West


  # Vector Error Correction Model is as follows
  ?buildVECM
  buildVECM(Pgasoline ~ Pmax + Prec + Pcut, newdat.ts,  
       stationary.vars = ~ Utilization + Stocks) # add the stationary variables like so
  }
}
