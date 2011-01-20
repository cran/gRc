\name{fitInteral}
\alias{scoring}
\alias{scoring.rcox}
\alias{ipm}
\alias{ipm.rcon}
\alias{ipm.rcor}
\alias{matching}
\alias{matching.rcon}
\alias{matching.rcor}

\alias{rconIPM}
\alias{rcorIPM}

\alias{rconScoreMatch}
\alias{rcorScoreMatch}
\alias{rconScoreTheta}
\alias{rcorScoreTheta}

\alias{fitIPSedge}
\alias{fitIPSset}
\alias{fitNR2}
\alias{fitNR}
\alias{modNewt}
\alias{refitA}
%%\alias{}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions used in connection with fitting of RCOX models}
\description{
  This is an overview over core functions used in connection with
  fitting of RCOX models. The functions described here will typically
  not be directly called by the user. 
}
\usage{
scoring (object, K0, control = object$control, maxit=control$maxouter, trace = object$trace)
ipm     (object, K0, control = object$control, trace=object$trace)
matching(object, control=object$control, trace=object$trace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model object (an object of class rcox.Internal)}
  \item{K0}{An initial value for K}
  \item{control}{A list controlling the fitting algorithms.}
  \item{maxit}{Maximal number of scoring iterations}
  \item{trace}{Controls various diagnostics print outs. A debugging
    feature not intended for the user.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with fitted values (K, logL etc).
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ %
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{\code{\link{fit}}, \code{\link{rcox}}}
\examples{
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)
data(math)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)
f1 <- matching(m1)

## Use f1$K as starting value
scoring(m1, K0=f1$K)
ipm(m1, K0=f1$K)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
