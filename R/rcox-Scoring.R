## Multivariate Newton
## 

## Start
#   Get start value for theta(0) = fn(K(0))
#   Get f
#   Calc: logL(0)
#

## At iteration p; given theta(p-1) and K(p-1)
#   Set:  stepsize = 1
#   Get (score, information) = f(K,...)
#   Calc: updating 'discrepancy' DISC 

## Linesearch:
#   Calc: theta(p) = theta(p-1) + DISC*stepsize
#   Calc: K(p) = fn(theta(p))
#   Calc: logL(p)

# If logL(p)<logL(p-1)
#   Set: stepsize = stepsize / 2
#   Do linesearch

# Stop doing lineserach when steps have been halved more than 3 times...

scoring.RCOX <- function(m, K0, ## =fitInfo(m, "K"),
                         control=m$control,
                         maxit=control$maxouter,
                         trace=m$trace){


  tstart <- proc.time()
  if (trace>=2)cat("..Fitting model with scoring\n")
  
  f       <- dataRep(m, 'n') - 1
  S       <- dataRep(m,'S')
  V     <- NULL

  #Kdiag <- diag(K0)
  #Krest <- K0-Kdiag
  ##Kmod <- Kdiag + 0.95*Krest
  #K0 <- Kmod

  logL0 <- ellK(K0, S, f)

  theta0      <- K2theta(m, K0, scale='free')
  curr.logL   <- logL0
  curr.theta  <- theta0
  
  curr.K      <- K0
  logLeps     <- 1e-6 * abs(logL0)
  ##print(logLeps)
  logL.vec    <- rep(NA, maxit)

  ## Iterate here...
  if (trace>=3) cat("...Scoring start:", "logL:", curr.logL, "\n")
  itcount <- 0
  repeat {
    x       <- getScore(m, curr.K, scale='free') 
    #print(x)
    Sc      <- x$score
    J       <- x$J
    ## print(J)
    ## stop()
    ##print(solve(J))
    DISC    <- try(qr.solve(J+ (Sc%*%t(Sc)/sqrt(f)), Sc)) ## /f
    if (class(DISC)=='try-error'){
      cat("Error in Fisher scoring, please report...\n")
    }

    stepsize      <- 1
    stephalfcount <- 0
    repeat{
      new.theta   <- curr.theta + DISC*stepsize
      new.K       <- theta2K(m, new.theta, scale='free')

      new.logL    <- ellK(new.K, S, f)
      diff.logL   <- new.logL - curr.logL
      ##print(diff.logL)
      
      if (is.na(new.logL) | (diff.logL < 0 & stephalfcount < 4)){
        stephalfcount <- stephalfcount + 1
        stepsize <- stepsize / 2
        if(trace>=4) cat ("....New stepsize:", stepsize, "\n")
      } else {
        break()
      }
    }

    itcount     <- itcount + 1
    if (trace>=3) cat("...scoring iteration:", itcount, "logL:", new.logL, "\n")
    
    curr.logL   <- new.logL
    curr.theta  <- new.theta
    curr.K      <- new.K
    logL.vec[itcount] <- curr.logL
    
    if (diff.logL < logLeps | itcount >= maxit){
      break()
    }
  }

  if (trace>=3) cat("...Scoring iterations:", itcount, "\n")
  logL.vec <- logL.vec[!is.na(logL.vec)]

  vn <- unlist(lapply(getcc(m),names))
  names(new.theta) <- vn

  ## Back to original scale
  l              <- length(getSlot(m, "vcc"))
  new.theta[1:l] <- exp(new.theta[1:l])  ### CHECK - er det rigtigt for rcor???

  dimnames(new.K) <- dimnames(S)

  #print("GET SCORE AT THE END...")
  J       <- getScore(m, new.K, scale='original')$J

#   print("J"); print(J)
  
#   if (!is.null(control$vcov)){
#     V <- calculateVCOV(m, K=new.K, vcov=control$vcov, nboot=control$nboot)
#   }

  dimnames(J) <- list(vn, vn)

  ans      <- list(K=new.K, logL=new.logL, coef=new.theta, J=J, logL.vec=logL.vec)

  ##cat("scoring time:\n"); print(proc.time()-tstart)
  
  return(ans)
}



