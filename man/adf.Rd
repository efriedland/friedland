\name{adf}
\alias{adf}
\title{ Augmented Dickey-Fuller Test }
\description{ The Augmented Dickey-Fuller test is for unit root testing. Remember that the P-Value from traditional summary() of the regression of is often spurious. Uses MacKinnon (1994) Table 3 "“Approximate asymptotic distribution functions for unitroot and cointegration tests” to determine significance of the ADF test statistic result. }
\usage{ adf(variable, k = 1, drift = TRUE, trend = TRUE, startfrom = 0, ADFstatOnly = FALSE) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{ A Time Series object. }
  \item{k}{ The maximum number of lags desired. Use only a single positive integer. }
  \item{drift}{ Logic statement. If TRUE, includes an intercept/constant in the model. }
  \item{trend}{ Logic statement. If TRUE, includes a time trend in the model. }
  \item{startfrom}{ The start date of the sample in the test, can be used to limit observations }
  \item{ADFstatOnly}{ Logic statement. If TRUE, the function only returns the Augmented Dickey Fuller statistic. }
}
\details{ 
regular text
\emph{italicized text}
\bold{bold text}

A new paragraph
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ MacKinnon (1994) }
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{UnitRoot}}, \code{\link{UnitRootApply}}
}
\examples{
GermanInterest.ts <- as.ts(dynlm::M1Germany[,"interest"])
x <- adf(GermanInterest.ts, k = 4, drift = T, trend = T)
summary(x)
}
