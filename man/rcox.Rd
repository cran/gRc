\name{rcox}
\alias{rcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main function for specifying RCON/RCOR models}
\description{
  This is the main function for specifying and fitting RCON/RCOR models in the
  package along with certain utility functions. 
}
\usage{
rcox(gm = NULL, vcc = NULL, ecc = NULL, type = c("rcon", "rcor"),
method = "ipm",
fit = TRUE, data = NULL, S = NULL, n = NULL, Kstart, control = list(),
details=1, trace=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{Generating class for a grapical Gaussian model, see
    'Examples' for an illustration}
  \item{vcc}{List of vertex colour classes for the model}
  \item{ecc}{List of edge colour classes for the model}
  \item{type}{Type of model. Default is RCON}
  \item{method}{Estimation method; see 'Details' below. }
  \item{fit}{Should the model be fitted}
  \item{data}{A dataframe}
  \item{S}{An empirical covariance matrix (as alternative to giving data
  as a dataframe)}
  \item{n}{The number of observations (which is needed if data is
    specified as an empirical covariance matrix)}
  \item{Kstart}{An initial value for K. Can be omitted.}
  \item{control}{Controlling the fitting algorithms}
  \item{details}{Controls the amount of output}
  \item{trace}{Debugging info}
  %\item{object}{An RCOX model object.}
  %\item{slot}{A name of a slot. }
  %\item{complement}{If FALSE, the edges of the model is returned. If
  %  TRUE, the edges not in the model is returned}
}
 \details{
   Estimation methods:

   'ipm' (default) is iterative partial maximization
   which when finished calculates the information matrix so that
   approximate variances of the parameters can be obtained using vcov().

   'ipms' is iterative partial maximization without
   calculating the information matrix. This is the fastest method.
   
   'scoring' is stabilised Fisher scoring.
   
   'matching' is score matching followed by one step with Fisher
   scoring.  

   'hybrid1' is for
   internal use and should not be called directly
 }
\value{
  A model object of type 'RCOX'.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
\note{
  demo("gRc-JSS") gives a more comprehensive demo. 
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching')
m2 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='scoring')
m3 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='ipm')

m1
m2
m3

summary(m1)
summary(m2)
summary(m3)

coef(m1)
coef(m2)
coef(m3)

vcov(m1)
vcov(m2)
vcov(m3)
}
\keyword{models}% at least one, from doc/KEYWORDS
