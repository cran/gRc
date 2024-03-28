#' @useDynLib gRc

## Vanilla R imports (and exports)
## -------------------------------

#' @importFrom stats as.formula cov cov2cor delete.response deviance
#'     formula runif terms xtabs simulate xtabs runif terms addmargins
#'     as.formula cov.wt fitted formula ftable getCall logLik loglin
#'     na.omit pchisq pf pnorm r2dtable terms update update.formula
#'     coef vcov
#'
#' @importFrom grDevices heat.colors topo.colors
#'
#' @importFrom MASS mvrnorm
#' 
#' @importFrom utils combn str install.packages
#'
#' @importMethodsFrom stats4 plot
#' @exportMethod plot
#' 

## Miscellaneous
## -------------
#' @importFrom Rcpp evalCpp
#'
#' @import methods
#' @import gRbase
#'
#' @importFrom igraph 
#'     get.adjacency V "V<-" E "E<-" is.directed layout.lgl
#'     layout.graphopt plot.igraph as_edgelist get.edge.ids

 
NULL

