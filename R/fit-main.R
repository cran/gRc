###################################################################### 
#' @title Fit rcox model
#' @description This is a general function for fitting RCOX models
#'   (i.e.\ RCON and RCOR models) using different estimation
#'   algorithms.
#' @author  Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name fit-main
###################################################################### 
#'
#' @param object An RCOX model object (an object of class 'rcox')
#' @param Kstart An initial value for the concentration matrix.
#' @param method The specific estimation method. Can be either
#'   "scoring", (a modified Fisher scoring algorithm), "ipm"
#'   (iterative partial maximization), "matching" (score matching) or
#'   "user" (currently not used).
#' @param K0 Initial value for concentration matrix.
#' @param control A list controlling the fitting algorithms. See the
#'   'details' section.
#' @param details The amount of details printed on the screen. 0 means
#'   no details at all.
#' @param trace Controls various diagnostics print outs. A debugging
#'   feature not intended for the user.
#' @param returnModel If TRUE the model object m is returned with
#'   fitting info added to it. If FALSE only the fitting info is
#'   returned.
#' @param ... Additional arguments; currently not used.
#' 
#' @details The fitted parameters etc. can be extracted using
#'   'fitInfo(m)'.  The control argument is a list with named
#'   entries. Most important are the entries 'maxouter' and 'maxinner'
#'   (which both defaults to 25) for controlling the estimation
#'   algorithms. For other components please  refer to the code. 
#' 
#' @keywords models 
#' @return An RCOX model object.
#' @seealso \code{\link{rcox}}, \code{\link{update.rcox}}
#' @examples
#' 
#' data(math)
#' gm  = ~al:an:st
#' vcc = list(~me + st, ~ve + an, ~al)
#' ecc = list(~me:ve + me:al, ~ve:al + al:st)
#' 
#' m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)
#' 
#' fit(m1, method="matching")
#' fit(m1, method="scoring")
#' fit(m1, method="ipm")
#'

#' @rdname fit-main
#' @export
fit.rcox <- function(object,
                     Kstart  = object$Kstart,
                     method  = object$method,
                     control = object$control,
                     details = object$details,
                     trace   = object$trace,
                     returnModel = TRUE,...){

  ##  cat("fit.rcox\n")
  if (is.null(Kstart)){
    #cat("Finding Kstart\n")
    Kstart    <- matching(object, trace=trace)$K
    #print(Kstart)
  }
  
  ##cat("Kstart:\n"); print(Kstart)
  ##Kstart    <- findKinModel(object, KS=object$Kstart,type=object$type, regularize=TRUE)
  
  ##print(method); print(trace)
  tstart <- proc.time()
  ans <- switch(method,
    "matching"=
      {
        scoring(object, K0=Kstart, control=control, maxit=1, trace=trace)
      },
    "scoring"=,
    "ipm" =,
    "ipms"=
      {
        switch(method,
          "scoring"={
            #print(Kstart)
            scoring(object, K0=Kstart, control=control, trace=trace)
          },
          "ipm"=,
          "ipms"={
            #print(Kstart)
            ipm(object, K0=Kstart, control=control, trace=trace)         
          })
      },
    "hybrid1"={
      object2 <- object
      ctrl          <- object$control
      ctrl$maxouter <- ctrl$hybrid1switch
      ctrl$vcov     <- NULL
      KK  <- ipm(object2, K0=Kstart, control=ctrl, trace=trace)$K
      scoring(object, K0=KK, control=control, trace=trace)
    }
  )
  ans$method <- method
  ans$Kstart <- Kstart 
  ans$time   <- (proc.time()-tstart)[3]
  
  if (returnModel){
    object$Kstart   <- ans$Kstart
    ans$Kstart      <- NULL    
    object$fitInfo  <- ans
    object$method   <- method
    return(object)
  } else {
    return(ans)
  }
}


#' @rdname fit-main
#' @export
matching      <- function(object, control=object$control, trace=object$trace){
  if (inherits(object, "rcon"))
    rconScoreMatch(object, control=control, trace=trace)
  else
    rcorScoreMatch(object, control=control, trace=trace)
}

#' @rdname fit-main
#' @export
ipm <- function(object, K0, control=object$control, trace=object$trace){
  if (inherits(object, "rcon"))
    rconIPM(object, K0, control, trace, engine="r")
  else
    rcorIPM(object, K0, control, trace, engine="r")
}

rconIPM <- function(object, K0, control=object$control, trace=object$trace, engine="R"){
  engine=match.arg(tolower(engine), c("r", "cpp"))
  switch(engine,
    "r"  = rconIPM_r(object, K0, control, trace),
    "cpp"= stop("Enigne 'cpp' currently not available.")
    )
}

rcorIPM <- function(object, K0, control=object$control, trace=object$trace, engine="R"){
  engine=match.arg(tolower(engine), c("r", "cpp"))
  switch(engine,
    "r"  = rcorIPM_r(object, K0, control, trace),
    "cpp"= stop("Enigne 'cpp' currently not available.")
    )
}















## matching.rcon <- function(object, control=object$control, trace=object$trace){
##   rconScoreMatch(object, control=control, trace=trace)
## }

## matching.rcor <- function(object, control=object$control, trace=object$trace){
##   rcorScoreMatch(object, control=control, trace=trace)
## }

## ipm.rcon <- function(object, K0, control=object$control, trace=object$trace){
##   rconIPM_r(object, K0, control, trace)
## }

## ipm.rcor <- function(object, K0, control=object$control, trace=object$trace){
##   rcorIPM_r(object, K0, control, trace)
## }










