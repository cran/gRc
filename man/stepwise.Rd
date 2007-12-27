\name{stepwise}
\alias{stepadd1}
\alias{stepdrop1}
\alias{stepjoin1}
\alias{stepsplit1}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepwise model selection in RCOX models}
\description{
  These allow for stepwise model selection in RCOX models by. Model
  expansion (i.e. forward selection) is obtained by adding edge colour
  classes and by splitting edge/vertex colour classes.
  Model reduction (i.e. backward selection) is obtained by dropping edge
  colour classes and by joining edge/vertex colour classes.
}
\usage{
stepadd1  (object, criterion = "aic", steps = 1000, k = 2, alpha = 0.05, details=1, trace=0)
stepdrop1 (object, criterion = "aic", steps = 1000, k = 2, alpha = 0.05, stat = "wald", details=1, trace=0) 
stepjoin1 (object, type = "ecc", criterion = "aic", steps = 1000, k = 2, alpha = 0.05, stat = "wald",   details = 1, trace=0) 
stepsplit1(object, type = "ecc", criterion = "aic", steps = 1000, k = 2, alpha = 0.05, stat = "wald",   details = 1, trace=0) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model, an object of class RCOX}
  \item{criterion}{Either "aic" (the default), "bic" or "test" (for
    significance test)}
  \item{type}{Either "ecc" for edge colour classes or "vcc" for vertex
    colour classes.}
  \item{k}{The multiple of the number of degrees of freedom used for the
    penalty when criterion is "aic". Ignored when criterion is "bic" or
    "test". Only k = 2 gives the genuine AIC.}
  \item{steps}{The maximum number of steps to be considered. The default
    is 1000 (essentially as many as required). It is typically used to
    stop the process early} 
  \item{stat}{Either "wald" for a Wald statistic or "dev" for a
    deviance statistic.}
  \item{alpha}{Critical value if 'criterion' is "test". If criterion is
    "aic" or "bic", the critical value is 0.}
  \item{details}{Control the amount of output created.}
  \item{trace}{For debugging purposes}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Either NULL or a new RCOX model.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}

%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\note{
  Note that the keyword 'stat' is not available for stepadd1  and
  stepsplit1 because these functions expand
  the current model and hence the Wald statistic is not available. 
}

\seealso{
  \code{\link{split1}}
  \code{\link{join1}}
  \code{\link{add1.rcox}}
  \code{\link{drop1.rcox}}
  \code{\link{comparecc}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
