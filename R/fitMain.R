
fit.rcox <- function(m,
                     Kstart  = m$Kstart,
                     method  = m$method,
                     control = m$control,
                     details = m$details,
                     trace   = m$trace,
                     returnModel=TRUE,...){
  
  if (is.null(Kstart)){
    ## cat("Finding Kstart\n")
    Kstart    <- matching(m, trace=trace)$K
  }
                                        #  else
#    Kstart    <- findKinModel(m, KS=m$Kstart,type=m$type, regularize=TRUE)

  tstart <- proc.time()
  ans <- switch(method,
                "matching"=
                {
                  scoring(m, K0=Kstart, control=control, maxit=1, trace=trace)
                },
                "scoring"=,
                "ipm"=
                {
                  switch(method,
                         "scoring"={
                                        #print(Kstart)
                           scoring(m, K0=Kstart, control=control, trace=trace)
                         },
                         "ipm"={
                                        #print(Kstart)
                           ipm(m, K0=Kstart, control=control, trace=trace)         
                         })
                },
                "hybrid1"={
                  m2 <- m
                  ctrl          <- m$control
                  ctrl$maxouter <- ctrl$hybrid1switch
                  ctrl$vcov     <- NULL
                  KK  <-ipm(m2, K0=Kstart, control=ctrl, trace=trace)$K
                  scoring(m, K0=KK, control=control, trace=trace)
                }
                )
  ans$method <- method
  ans$Kstart <- Kstart 
  ans$time   <- (proc.time()-tstart)[3]
  
  if (returnModel){
    m$Kstart   <- ans$Kstart
    ans$Kstart <- NULL    
    m$fitInfo  <- ans
    m$method   <- method
    return(m)
  } else {
    return(ans)
  }
}

matching      <- function(m, control=m$control, trace=m$trace){
  if (inherits(m,"rcon"))
    rconScoreMatch(m, control=control, trace=trace)
  else
    rcorScoreMatch(m, control=control, trace=trace)
  ##UseMethod("matching")
}

ipm <- function(m, K0, control=m$control, trace=m$trace){
  if (inherits(m,"rcon"))
    rconIPM(m, K0, control, trace)
  else
    rcorIPM(m, K0, control, trace)
  ##UseMethod("ipm")
}


#scoring <- function(m, K0, control=m$control, maxit=control$maxouter,trace=m$trace) {
#  UseMethod("scoring")
#}

matching.rcon <- function(m, control=m$control, trace=m$trace){
  rconScoreMatch(m, control=control, trace=trace)
}

matching.rcor <- function(m, control=m$control, trace=m$trace){
  rcorScoreMatch(m, control=control, trace=trace)
}

ipm.rcon <- function(m, K0, control=m$control, trace=m$trace){
  rconIPM(m, K0, control, trace)
}

ipm.rcor <- function(m, K0, control=m$control, trace=m$trace){
  rcorIPM(m, K0, control, trace)
}









