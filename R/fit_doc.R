###################################################################### 
#' @title Functions used in connection with fitting of RCOX models
#' @description This is an overview over core functions used in
#'   connection with fitting of RCOX models. The functions described
#'   here will typically not be directly called by the user.
#' @author  Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name fit-doc
###################################################################### 
##' 
##' @aliases ipm.rcon ipm.rcor
##'   matching.rcon matching.rcor rconIPM_c rconIPM_r rcorIPM_c
##'   rcorIPM_r rconScoreMatch rcorScoreMatch rconScoreTheta
##'   rcorScoreTheta fitIPSedge fitIPSset fitNR2 fitNR modNewt refitA
##'   refitA_old
##'
##' @param object An RCOX model object (an object of class rcox.Internal)
##' @param K0 An initial value for K
##' @param control A list controlling the fitting algorithms.
##' @param maxit Maximal number of scoring iterations
##' @param trace Controls various diagnostics print outs. A debugging feature
##' not intended for the user.
##'
##' @return A list with fitted values (K, logL etc).
##' @seealso \code{\link[gRbase]{fit}}, \code{\link{rcox}}
##' @keywords internal
##' @examples
##' 
##' gm  = ~al:an:st
##' vcc = list(~me+st, ~ve+an, ~al)
##' ecc = list(~me:ve+me:al, ~ve:al+al:st)
##' data(math)
##' 
##' m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)
##' f1 <- matching(m1)
##' 
##' ## Use f1$K as starting value
##' scoring(m1, K0=f1$K)
##' ipm(m1, K0=f1$K)
##'
NULL
