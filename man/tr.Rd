\name{tr}
\alias{trA}
\alias{trAW}
\alias{trAWB}
\alias{trAWBW}
\alias{trAWBV}



\title{Calculate trace of various matrix products}
\description{
  Calculate trace of various matrix products.
}
\usage{
trA(A)
trAW(A,W)
trAWB(A,W,B)
trAWBW(A,W,B)
trAWBV(A,W,B,V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A,B}{Square matrices represented as matrices or lists (see examples below).}
  \item{W,V}{Square matrices}
}
%\details{
%
%}
\value{
  A number
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{
%  The storage mode of the matrices must be 'double', see examples below
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
d <- 5
W <- matrix(rnorm(d*d),nr=d,nc=d); 
V <- W <- W+t(W)

## Turn list into matrix
##
tomat <- function(x){
  ans <- do.call("rbind", x)
  storage.mode(ans)<-"double"
  return(ans)
}

A1 <- tomat(list(c(1,2),c(1,3)))
A2 <- tomat(list(1,3,5))


## Just for checking the calculations
##
symMat <- function(A,d){
  ans <- matrix(0,nr=d,nc=d)
  for (i in 1:length(A)){
    e <- A[[i]]
    if (length(e)==1){
      ans[e,e] <- 1
    } else { 
      ans[e[1],e[2]] <-   ans[e[2],e[1]] <- 1 
    }
  }

  return(ans)
}

trAW(A1, W)
#sum(diag(symMat(A1,d=d) %*% W))

trAW(A2, W)
#sum(diag(symMat(A2,d=d) %*% W))

trAWB(A1, W, A2)
#sum(diag(symMat(A1,d=d) %*% W %*% symMat(A2,d=d)))

trAWBV(A1, W, A2, V)
#sum(diag(symMat(A1,d=d) %*% W %*% symMat(A2,d=d) %*% V))
}
\keyword{utilities}% at least one, from doc/KEYWORDS
