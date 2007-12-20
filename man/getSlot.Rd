\name{getSlot}
\alias{getSlot}
\alias{getSlot}
\alias{fitInfo}
\alias{dataRep}
\alias{intRep}
\alias{getedges}
\alias{getecc}
\alias{getvcc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessing RCOX model objects}
\description{
  Accessing RCOX model objects
}
\usage{
## Accessor functions
getSlot(object,slot)
fitInfo(object,slot)
intRep(object,slot)
dataRep(object,slot)
getecc(object)
getvcc(object)
getedges(object, complement=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model object.}
  \item{slot}{A name of a slot. }
  \item{complement}{If FALSE, the edges of the model is returned. If
    TRUE, the edges not in the model is returned}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
% \value{
%   A list. 
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
% \references{ ~put references to the literature/web site here ~ }

\author{Søren Højsgaard, sorenh@agrsci.dk}

%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ %
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{\code{\link{rcox}}}
\examples{
data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math)
getecc(m1)

getSlot(m1,"type")
fitInfo(m1)
fitInfo(m1,"K")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
