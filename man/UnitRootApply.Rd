\name{UnitRootApply}
\alias{UnitRootApply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Unit Root Work Flow on a Dataset}
\description{
Puts the results of the UnitRoot() function into a matrix.
}
\usage{ UnitRootApply(data, k, ...) }
\arguments{
  \item{data}{ The Time Series dataset, in mts, ts, matrix format. }
  \item{k}{ The maximum lags desired, Use only a single positive integer. }
  \item{\dots}{ All arguments passed to \code{UnitRoot} function. }
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## UnitRootApply(data.ts, k = 4, drift = T, trend = F)
## UnitRootApply(data.ts[,c("Pcrude","Pgasoline")], k = 4, drift = T, trend = F)
}