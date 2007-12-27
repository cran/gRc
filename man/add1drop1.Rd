\name{add1drop1}
\alias{add1.rcox}
\alias{drop1.rcox}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Add or drop colour classes to RCOX models}
\description{
  Make a test for adding/dropping all colour classes in scope for an
  RCOX model.
}
\usage{
\method{add1}{rcox}(object, scope, details = 1, trace = 0, ...)
\method{drop1}{rcox}(object, scope, details = 1, trace = 0, stat = "wald", ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model, an object of class 'rcox'}
  \item{scope}{A set of edge colour classes to be considered for
    addition or deletion, see 'details'.}
  \item{details}{Control the amount of output created.}
  \item{trace}{For debugging purposes}
  \item{stat}{Either "wald" for a Wald statistic or "dev" for a
    deviance statistic.}
  \item{\dots}{Additional arguments, currently unused.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with entries:
  \item{tab}{A data frame with the test results}
  \item{cc}{A list of colour classes}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
\note{
  Note that the keyword 'stat' is not available for add1 because this function expands
  the current model and hence the Wald statistic is not available. 
}

\seealso{
  \code{\link{comparecc}},
  \code{\link{stepadd1}},
  \code{\link{stepdrop1}}
}
\examples{
data(math)
gc.sat <- ~me:ve:al:st:an
gc.1   <- ~me+ve+al+st+an

m.sat <- rcox(gm=gc.sat, data=math)
m.1   <- rcox(gm=gc.1,   data=math)

t.sat <- drop1(m.sat)
t.sat$tab
t.sat$cc

t.1   <- add1(m.1)
t.1$tab
t.1$cc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
