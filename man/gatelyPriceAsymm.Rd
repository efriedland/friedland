\name{gatelyPriceAsymm}
\alias{gatelyPriceAsymm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Decompose Prices
}
\description{
Produces a matrix of 3 variables decomposed from the input time series variable. The historical maximum at each date, 
  the cumulative positive changes, and the cumulative negative changes. 
}
\usage{
gatelyPriceAsymm(x)
}
\arguments{
  \item{x}{ The Time Series variable to be decomposed. }
  \item(prefix){ The string that leads "max" "rec" and "cut" which is by default "P"}
}
\value{
%%  A time series matrix containing "Pmax" "Prec" and "Pcut" where the "P" can be changed using the prefix argument. 
}
\references{
Gately (1991) "Imperfect Price-Reversibility of U.S. Gasoline Demand: Asymmetric Responses to Price Increases and Declines.
}
\seealso{
\code{\link{Asymmetry}}
}
\examples{
# Data
Asymmetry.ts <- ts(Asymmetry[,-1], start = c(2003,6), frequency = 12)
decomposed <- gatelyPriceAsymm(Asymmetry.ts[,"Pcrude"], prefix = "Pcrude_")
newdat.ts <- ts.intersect(Asymmetry.ts, decomposed)
colnames(newdat.ts) <- c(colnames(Asymmetry.ts), colnames(decomposed))
head(newdat.ts)

# 2 models
lessrestricted <- dynlm(Pgasoline ~ Pcrude_max + Pcrude_rec + Pcrude_cut, newdat.ts)
summary(lessrestricted)
restricted <- dynlm(Pgasoline ~ Pcrude, newdat.ts, 
                    start = start(lessrestricted)) # fix the dates to be equal
summary(restricted)

# Choosing between them
anova(dynlm(Pgasoline ~ Pcrude_max + Pcrude_rec + Pcrude_cut, newdat.ts), 
      dynlm(Pgasoline ~ Pcrude, newdat.ts, start = c(2003,7)))
# F test shows significant difference. Use less restricted with smaller RSS. 
}
