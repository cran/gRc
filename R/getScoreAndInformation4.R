getScore   <- function(m, K, scale='original'){
  UseMethod('getScore')
}


getScore.rcon <- function(m, K, scale='original'){ ### OK !!!

  md    <- m$dataRep
  ir    <- m$intRep
  
  S     <- md$S;  
  Sigma <- solve(K)

  vccTerms <- ir$vccI
  eccTerms <- ir$eccI
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc

  score    <- rep(NA, nparm)
  J        <- matrix(0,nrow=nparm, ncol=nparm)

  f <- md$n-1; f2 <- f/2

  DSigma <- Sigma-S
  
  for (u in 1:lvcc){
    Ku <- vccTerms[[u]]
    val <-  f2 * trAW(Ku, DSigma)
    if (scale=='free')
      score[u] <- val * K[Ku[[1]],Ku[[1]]]
    else
      score[u] <- val 
  }    

  if (lecc>0){
    for (u in 1:lecc){
      Ku  <- eccTerms[[u]]
      val <-  f2 * trAW(Ku, DSigma)
      score[u+lvcc] <- val
    }    
  }  

  for (u in 1:lvcc){    ##print("V,V")
    Ku  <- vccTerms[[u]]
    for (v in u:lvcc){
      Kv  <- vccTerms[[v]]
      val <- f2* trAWBW(Ku, Sigma, Kv)      
      if (scale=='free')
        J[v,u] <- J[u,v] <- val * (K[Ku[[1]],Ku[[1]]]*K[Kv[[1]],Kv[[1]]])
      else
        J[v,u] <- J[u,v] <- val 
    }    
  }

  if (lecc>0){    ##print("V,E")    
    for (u in 1:lvcc){
      Ku  <- vccTerms[[u]]
      for (v in 1:lecc){
        Kv  <- eccTerms[[v]]
        val <- f2* trAWBW(Ku, Sigma, Kv)      
        if (scale=='free')
          J[u,v+lvcc] <- J[v+lvcc,u] <- val * K[Ku[[1]],Ku[[1]]]
        else
          J[u,v+lvcc] <- J[v+lvcc,u] <- val 
      }    
    }

    for (u in 1:lecc){    ##print("E,E")    
      Ku  <- eccTerms[[u]]
      for (v in u:lecc){
        Kv  <- eccTerms[[v]]
        J[v+lvcc,u+lvcc] <-J[u+lvcc,v+lvcc] <- f2* trAWBW(Ku, Sigma, Kv)      
      }    
    }
  }
  
  return(list(score=score, J=J))
}


getScore.rcor <- function(m,K, scale='original'){
  
  md    <- m$dataRep
  ir    <- m$intRep

  S     <- md$S;
  f     <- md$n - 1

  vccTerms <- ir$vccI
  eccTerms <- ir$eccI
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc

  score    <- rep(NA, nparm)
  J        <- matrix(0,nrow=nparm, ncol=nparm)

  C        <- cov2cor(K); 
  Cinv     <- solve(C)             ## Brug IKKE cholSolve(C) - numerisk ustabil

  a        <- sqrt(diag(K))        ## a indeholder eta'erne
  A        <- diag(a)              

  ASA      <- a* t(a*S)
  CASA     <- C %*% ASA
  ACAS     <- (a * C) %*% (a * S)

  ## Run through VCC:
  for (u in 1:lvcc){
    
    ## Score for VCC x VCC :
    Ku  <-  vccTerms[[u]]
    val <-  f*(trA(Ku) - trAW(Ku, ACAS)) 

    if (scale=='free')
      score[u] <- val
    else
      score[u] <- val / (A[Ku[[1]],Ku[[1]]]) ## OK, dec 07


    ## Information for VCC
    for (v in 1:lvcc){
      Kv  <- vccTerms[[v]]      

      if (u==v)
        val <- 2*f* trAWBV(Ku, Cinv, Kv, C) ## OK, sept 07
      else
        val <- f * trAWBV(Ku, Cinv, Kv, C) ## OK, sept 07
      
      if (scale=='original')
        val <- val / (A[Ku[[1]],Ku[[1]]]*A[Kv[[1]],Kv[[1]]])

      J[u,v] <- J[v,u] <-val
    }    
  }

  if (lecc>0){
    ## Score for ECC:
    for (u in 1:lecc){
      Ku <- term.u <- eccTerms[[u]]
      score[u+lvcc] <- (f/2) * (trAW(Ku, Cinv-ASA))
    }    
    
    for (u in 1:lvcc){
      ## Score for VCC x ECC
      Ku <- term.u <- vccTerms[[u]]
      for (v in 1:lecc){
        Kv <- term.v <- eccTerms[[v]]
        val <- f * trAWB(Ku, Cinv, Kv) ## OK, sept 07
        if (scale=='free')
          J[u,v+lvcc] <- J[v+lvcc,u] <- val
        else
          J[u,v+lvcc] <- J[v+lvcc,u] <- val /(A[Ku[[1]],Ku[[1]]]);
      }    
    }
    
    for (u in 1:lecc){
      ## Score for ECC x ECC
      Ku <- eccTerms[[u]]
      for (v in u:lecc){
        Kv  <- eccTerms[[v]]
        J[u+lvcc,v+lvcc] <- J[v+lvcc,u+lvcc] <- (f/2)* trAWBW(Ku, Cinv, Kv) ## OK, sept 07
      }    
    }
    
  }
  ##print(J)
  return(list(score=score, J=J))
}









