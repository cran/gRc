fit              <- function(object,
                             method  = object$method,
                             control = object$control,
                             details = object$details,
                             trace   = object$trace,
                             returnModel=TRUE
                             ){
  UseMethod("fit")
}

fit.RCOX <- function(object,
                     method  = object$method,
                     control = object$control,
                     details = object$details,
                     trace   = object$trace,
                     returnModel=TRUE){
  ans        <- .fitit(object, method=method, control=control,trace=trace)

  if (returnModel){
    object$fitInfo  <- ans
    object$method   <- method
    return(object)
  } else {
    return(ans)
  }
}


.fitit <-  function(m,
                    method  = m$method,
                    control = m$control,
                    details = m$details,
                    trace   = m$trace){
  if (trace>=1) cat(".Fitting method: ", method, "\n")
  tstart <- proc.time()
  ans <- switch(method,
                "matching"=
                {
                  ##matching(m, control=control, trace=trace)
                  ctrl      <- m$control;
                  ctrl$vcov <- NULL
                  Kstart    <- matching(m, control=ctrl, trace=trace)$K
                  scoring(m, K0=Kstart, control=control, maxit=1, trace=trace)
                  
                },
                "scoring"=,
                "ipm"=
                {
                  if (is.null(m$Kstart)){
                    ctrl      <- m$control;
                    ctrl$vcov <- NULL
                    Kstart    <- matching(m, control=ctrl, trace=trace)$K
                  } else {
                    Kstart <- m$Kstart
                  }
                  switch(method,
                         "scoring"={
                           scoring(m, K0=Kstart, control=control, trace=trace)
                         },
                         "ipm"={
                           ipm(m, K0=Kstart, control=control, trace=trace)         
                         })
                }
         )
  ans$method <-  method
  ans$time   <- (proc.time()-tstart)[3]
  return(ans)
}




matching      <- function(m, control=m$control, trace=m$trace){
  UseMethod("matching")
}

matching.rcon <- function(m, control=m$control, trace=m$trace){
  rconScoreMatch(m, control=control, trace=trace)
}

matching.rcor <- function(m, control=m$control, trace=m$trace){
  rcorScoreMatch(m, control=control, trace=trace)
}

ipm <- function(m, K0, control=m$control, trace=m$trace){
  UseMethod("ipm")
}

ipm.rcon <- function(m, K0, control=m$control, trace=m$trace){
  ##print("ipm.rcon"); print(trace)
  rconIPM(m, K0, control, trace)
}

ipm.rcor <- function(m, K0, control=m$control, trace=m$trace){
  #print("ipm.rcor"); print(trace)
  rcorIPM(m, K0, control, trace)
}

scoring <- function(m, K0, control=m$control, maxit=control$maxouter,trace=m$trace) {
  UseMethod("scoring")
}

