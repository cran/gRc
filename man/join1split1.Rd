\name{join1split1}
\alias{join1}
\alias{split1}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Joining and splitting of colour classes in RCOX models}
\description{
  Test for joining of two colour classes  (of a specific type) by testing if their
  corresponding parameters are not significantly different.
  Split a colour class and test how much this changes the fit of the
  model.
  }
\usage{
join1(object, scope=NULL, type = "ecc", details = 1,stat = "wald")
split1(object, scope=NULL, type = "ecc", details = 1) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model, an object of class RCOX}
  \item{scope}{A specification of colour classes which should be
    considered for joining/splitting. If NULL, then all colour classes
    are considered.}
  \item{type}{Either "ecc" for edge colour classes or "vcc" for vertex
    colour classes.}
  \item{stat}{Either "wald" for a Wald statistic or "dev" for deviance
    statistic.}
  \item{details}{Control the amount of output}
}
%\details{
  
  
%}
\value{
  A list with entries:
  \item{tab}{A data frame with the test results}
  \item{cc}{A list of colour classes}
}
\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
\note{
  Note that the keyword 'stat' is not available for split1 because this function expands
  the current and hence the Wald statistic is not available. 
  Note also that join1 is simply a wrapper for comparecc applied to edge
  colour classes.  
  }
\seealso{\code{\link{rcox}}, \code{\link{update}}, \code{\link{comparecc}}}
\examples{
data(math)
g1     <- ~me:ve:al+al:st:an
m1     <- rcox(gm=g1, data=math)
join1(m1)

gm  = ~al:an:st
vcc = list(~me+st, ~ve+an)
ecc = list(~me:ve+me:al, ~ve:al+al:st)
m2 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, type="rcon")
split1(m2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
