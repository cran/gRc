\name{update.rcox}
\alias{update.rcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update an RCOX model}
\description{
  update will update and (by default) re-fit an RCOX model.
  %It does this by extracting the call stored in the object, updating the
  %call and (by default) evaluating that call. Sometimes it is useful to
  %call update with only one argument, for example if the data frame has been corrected
}
\usage{
\method{update}{rcox}(object, vcc = NULL, ecc = NULL, splitecc = NULL,
splitvcc = NULL, joinvcc = NULL, joinecc = NULL, addecc = NULL, dropecc
= NULL, Kstart = NULL, fit = TRUE, control=NULL, trace = object$trace, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An RCOX model, an object of class RCOX}
  \item{vcc}{Specification of the vertex colour classes in the model}
  \item{ecc}{Specification of the edge colour classes in the model}
  \item{splitvcc}{Existing vertex colour class to be split}
  \item{splitecc}{Existing edge colour class to be split}
  \item{joinvcc}{Existing vertex colour classes to be joined}
  \item{joinecc}{Existing vertex colour classes to be joined}
  \item{addecc}{New edge colour classes to be added}
  \item{dropecc}{Existing vertex color classes to be dropped (deleted)}
  \item{Kstart}{A start value for K}
  \item{fit}{Should the updated model be fitted.}
  \item{control}{A list of control parameters.}
  \item{trace}{For debugging purposes}
  \item{...}{Additional arguments, currently not used. }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A new model object of class 'rcox'.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}

\section{Warning}{
  Only one of the arguments pertaining to edge colour classes (i.e. ecc,
  splitecc, joinecc, dropecc, addecc) should be applied at the
  time. Likewise for the arguments pertaining to the vertex colour
  classes.

  The result will otherwise be highly unpredictable and is likely to
  cause an error. 
}

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

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching', trace=0)

update(m1, joinvcc=list(~me+st, ~ve+an))
update(m1, joinecc=list(~al:an, ~an:st))

update(m1, splitvcc=~ve+an)
update(m1, splitecc=~me:ve+me:al)


update(m1, dropecc=list(~me:st+st:an,~al:an,~st:al))
update(m1, addecc=list(~an:me+st:ve))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
