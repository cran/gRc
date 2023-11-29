######################################################################
#' @title Calculate trace of various matrix products.
#' @description Calculate trace of various matrix products.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name clever_trace
######################################################################

#' @param rA,rB Square matrices represented as "sparse" matrices. The
#'   dimension is inferred from dimension of V,W.
#' @param Alist,Blist Lists of such "sparse" matrices.
#' @param A A matrix
#' @param mode 0 or 1.
#' @param rW,rV,W Square matrices; must be of same dimension.

#' @return A number
#' @keywords utilities
#'

#' @examples
#' 
#' ## Turn list into "sparse" matrix
#' to_sp <- function(x){
#'   ans <- do.call("rbind", x)
#'   storage.mode(ans)<-"double"
#'   return(ans)
#' }
#' 
#' ## Turn "sparse" matrix into dense matrix
#' ##
#' to_de <- function(A,d){
#'   if (nrow(A) < 1) stop("need non-empty matrix\n")
#'   ans <- matrix(0,nr=d,nc=d)
#'   for (i in 1:nrow(A)){
#'     e <- A[i,]
#'     if (length(e) == 1){
#'       ans[e,e] <- 1
#'     } else { 
#'       ans[e[1],e[2]] <-   ans[e[2],e[1]] <- 1 
#'     }
#'   }  
#'   return(ans)
#' }
#' 
#' d <- 5
#' W <- matrix(rnorm(d * d), nr=d, nc=d);
#' V <- matrix(rnorm(d * d), nr=d, nc=d); 
#' W <- W + t(W)
#' V <- V + t(V)
#' 
#' 
#' A1.lst <- list(c(1,2),c(1,3))
#' A2.lst <- list(1,3,5)
#' 
#' A1.sp <- to_sp(A1.lst)
#' A2.sp <- to_sp(A2.lst)
#' 
#' A1.de <- to_de(A1.sp, d)
#' A2.de <- to_de(A2.sp, d)
#' 
#' trAW(A1.sp, W)
#' sum(diag(A1.de %*% W))
#' 
#' trAW(A2.sp, W)
#' sum(diag(A2.de %*% W))
#' 
#' trAWB(A1.sp, W, A2.sp)
#' trAWB(A2.sp, W, A1.sp)
#' sum(diag(A1.de %*% W %*% A2.de))
#' 
#' trAWBW(A1.sp, W, A2.sp)
#' trAWBW(A2.sp, W, A1.sp)
#' sum(diag(A1.de %*% W %*% A2.de %*% W))
#' 
#' trAWBV(A1.sp, W, A2.sp, V)
#' trAWBV(A2.sp, W, A1.sp, V)
#' sum(diag(A1.de %*% W %*% A2.de %*% V))
#' 

#' @rdname clever_trace
#' @export
trA <- function(A){
  if (ncol(A)==1) return(nrow(A))
  else return(0)
}


#' @rdname clever_trace
#' @export
trAW <- trAW_

#' @rdname clever_trace
#' @export
trAWB <- trAWB_

#' @rdname clever_trace
#' @export
trAWBW <- trAWBW_

#' @rdname clever_trace
#' @export
trAWBV <- trAWBV_

#' @rdname clever_trace
#' @export
trAWBlist <- trAWBlist_

#' @rdname clever_trace
#' @export
trAWBWlist <- trAWBWlist_


## #' @rdname clever_trace
## #' @export
## trAW <- function(A, W){
##   .Call("trAW", A, W, PACKAGE="gRc")
## }

## #' @rdname clever_trace
## #' @export
## trAWB <- function(A,W,B){
##   .Call("trAWB", A, W, B, PACKAGE="gRc")
## }

## #' @rdname clever_trace
## #' @export
## trAWBW <- function(A,W,B){
##    .Call("trAWBW", A, W, B, PACKAGE="gRc")
## }

## #' @rdname clever_trace
## #' @export
## trAWBV <- function(A,W,B,V){
##   .Call("trAWBV", A, W, B, V, PACKAGE="gRc")
## }

## #' @rdname clever_trace
## #' @export
## trAWBlist <- function(Alist,W,Blist,mode){
##   .Call("trAWBlist", Alist, W, Blist,mode, PACKAGE="gRc")
## }

## #' @rdname clever_trace
## #' @export
## trAWBWlist <- function(Alist,W,Blist,mode){
##   .Call("trAWBWlist", Alist, W, Blist, mode,PACKAGE="gRc")
## }


