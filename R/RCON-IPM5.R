
rconIPM <- function(object, K0,
                     control=object$control, trace=object$trace){


  S       <- object$dataRep$S
  nobs    <- object$dataRep$n
  Kstart  <- object$Kstart
  
  ctrl <- object$control
  logL      = 0
  converged = 1

  deltaeps  = ctrl$deltaeps
  maxouter  = ctrl$maxouter
  maxinner  = ctrl$maxinner
  
  eccfit <- control$eccfit
  vccfit <- control$vccfit

  if (vccfit)
    vcc <- object$intRep$vcc
  else
    vcc <- NULL

  if (eccfit)
    ecc <- object$intRep$ecc
  else
    ecc <- NULL

  glist   <- c(vcc,ecc)

  logL0       <- prevlogL <- ellK(K0,S,nobs-1)
  logLeps     <- ctrl$logLeps * abs(logL0)

  Kwork<-.Call("rconipm", S=S, nobs=nobs-1, K0, Glist=glist, 
               maxouter=maxouter, maxinner=maxinner, 
               logL=logL, logLeps=logLeps, deltaeps=deltaeps,
               converged=converged, trace=0,
               PACKAGE="gRc")

  coef <- K2theta(object,Kwork, scale='original')
  vn   <- unlist(lapply(getcc(object),names))
  names(coef) <- vn

  if (object$method=="ipms"){ ## IPM without finalizing with finding score
    J <- NULL
  } else {
    if (!is.null(control$vcov)){
      J  <- getScore(object,K=Kwork)$J
      dimnames(J) <- list(vn, vn)
    } else {
      J <- NULL
    }
  }
  
  ans <- list(K=Kwork, logL=logL, coef=coef, J=J)
  return(ans)
}


