\name{gatelyPriceAsymm}
\alias{gatelyPriceAsymm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Decompose Prices
}
\description{
Produces a matrix of 3 variables decomposed from the input time series variable. The historical maximum at each date, the cumulative positive changes, and the cumulative negative changes. Named after the paper Gately (1991) "Imperfect Price-Reversibility of U.S. Gasoline Demand: Asymmetric Responses to Price Increases and Declines.
}
\usage{
gatelyPriceAsymm(x)
}
\arguments{
  \item{x}{ The Time Series variable to be decomposed. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gately (1991) "Imperfect Price-Reversibility of U.S. Gasoline Demand: Asymmetric Responses to Price Increases and Declines.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{Asymmetry}}, ~~~
}
\examples{
## None right now
function (x) 
{
    stopifnot(is.ts(x))
    Pmax <- 0
    for (i in 1:length(x)) {
        Pmax[i] <- max(x[1:i])
    }
    Pchange <- diff(Pmax - x)
    Prec <- cumsum((Pchange > 0) * Pchange)
    Pcut <- cumsum((Pchange <= 0) * Pchange)
    Pmax <- ts(Pmax)
    tsp(Pmax) <- tsp(x)
    Prec <- ts(Prec)
    tsp(Prec) <- tsp(diff(x))
    Pcut <- ts(Pcut)
    tsp(Pcut) <- tsp(diff(x))
    list(Pmax = Pmax, Prec = Prec, Pcut = Pcut)
  }
}
