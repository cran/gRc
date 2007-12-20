\name{fit}
\alias{fit}
\alias{fit.RCOX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit RCOX models}
\description{
  This is a general function for fitting RCOX models (i.e.\ RCON and
  RCOR models) using different estimation algorithms.
}
\usage{
fit(object, method = object$method, control = object$control,
details = object$details, trace = object$trace, returnModel = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model object (an object of class RCOX)}
  \item{method}{The specific estimation method. Can be either "scoring",
    (a modified Fisher scoring algorithm),
    "ipm" (iterative partial maximization), "matching" (score matching)
    or "user" (currently not used)}
  \item{control}{A list controlling the fitting algorithms. See the
    'details' section.}
  \item{details}{The amount of details printed on the screen. 0 means no
    details at all.}
  \item{trace}{Controls various diagnostics print outs. A debugging
    feature not intended for the user.}
  \item{returnModel}{If TRUE the model object m is returned with fitting
    info added to it. If FALSE only the fitting info is returned.}
}
\details{
  The fitted parameters etc. can be extracted using 'fitInfo(m)'.

  The control argument is a list with named entries. Most important are
  the entries 'maxouter' and 'maxinner' (which both defaults to 25) for
  controlling the estimation algorithms. For other components please
  refer to the code. 
}
\value{
  An RCOX model object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
  \code{\link{rcox}},
  \code{\link{update.RCOX}}
}
\examples{

data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)

fit(m1, method="matching")
fit(m1, method="scoring")
fit(m1, method="ipm")

## MISSING
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
