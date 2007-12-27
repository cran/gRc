\name{comparecc}
\alias{comparecc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare colour classes of an RCOX model}
\description{
  A general function for  pairwise comparisons of
  colour classes in an RCOX model, i.e. for testing whether the corresponding parameters
  are significantly different
}
\usage{
comparecc(object, cc1 = NULL, cc2 = NULL, type = "ecc", stat = "wald", details = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model, an object of class 'rcox'}
  \item{cc1, cc2}{Lists of colour classes of type 'type', see 'details'
    for an explanation of the defaults.}
  \item{type}{Either "ecc" for edge colour classes or "vcc" for vertex
    colour classes}
  \item{stat}{Base the comparison on either "wald" for a Wald statistic
    or "dev" for a deviance statistic}
  \item{details}{Control the amount of output created.}
}
\details{
All colour classes specified in \code{cc1} are
compared with all those
given in \code{cc2} (duplicate entries are not compared). If
\code{cc2=NULL} (the default) then all colour classes specified in
\code{cc1} are compared with all colour classes in the model except
those specified in \code{cc1}. If \code{cc1=NULL} (the default) and
\code{cc2=NULL} then all pairwise comparisons are made.
}
\value{
  A list with entries:
  \item{tab}{A data frame with the test results}
  \item{cc1, cc2}{Lists of colour classes}
}
%%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
  \code{\link{add1.rcox}},
  \code{\link{drop1.rcox}},
  \code{\link{stepadd1}},
  \code{\link{stepdrop1}},
  \code{\link{join1}},
  \code{\link{split1}},
  \code{\link{stepjoin1}},
  \code{\link{stepsplit1}}
}
\examples{

data (math)

gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math)
m1

comparecc(m1, type="vcc")
comparecc(m1, type="ecc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
