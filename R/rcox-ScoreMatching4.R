
## iR <- getSlot(m, "intRep")
## vccTerms <- iR$vcc
## eccTerms <- iR$ecc
## allTerms <- c(vccTerms, eccTerms)

  #ctrl      <- control
  ##ctrl$vcov <- NULL ## We don't want vcov calculated until at the end...
  #m$control <- ctrl

rcorScoreMatch <- function(m, control=m$control, trace=0){
  if (trace>=2)
    cat("..Fitting with score matching\n")
  #theta     <- rcorScoreTheta(m)
  theta        <- rconScoreTheta(m)
  vn           <- unlist(lapply(getcc(m),names))
  names(theta) <- vn

  S        <- getSlot(m, "dataRep")$S
  n        <- getSlot(m, "dataRep")$n
  
  oclass       <- class(m)
  class(m)     <- c("rcon","rcox")
  K            <- theta2K(m, theta, scale='original')
  ## print("As RCON"); print(K)

  ## Ensure diagonals are positive
  if (min(diag(K))<0){  
    vccI <- intRep(m)$vccI
    dS <- diag(S)
    for (i in 1:length(vccI)){
      cc <- as.numeric(vccI[[i]])
      dS[cc] <- mean(dS[cc])
    }
    diag(K) <- 1/dS
  }
  ##print("After redoing diagonals");  print(K)

  if (min(eigen(K)$values)<0)  
    K     <- regularizeK(K)
  ##print("After regularizing diagonals"); print(K)

  dimnames(K)  <- dimnames(S)
  ##print(ellK(K,S,n-1))
  
  class(m) <- oclass
  K <- findKinModel(m,K,type="rcor")

  #print(K)  
  #print(ellK(K,S,n-1))
  #print(1/sqrt(diag(S)))
  
  ## K            <- theta2K(m, theta, scale="original")
  ## dimnames(K)  <- dimnames(S)

  logL   <- ellK(K, S, n-1)
  if (trace>=3)
    cat("..Score matching, logL:", logL, "\n")

  return(list(K=K, logL=logL, coef=theta))
}

rconScoreMatch <- function(m, control=m$control, trace=0){

  tstart <- proc.time()
  if (trace>=2)
    cat("..Fitting with score matching\n")

  theta <- rconScoreTheta(m)
  vn           <- unlist(lapply(getcc(m),names))
  names(theta) <- vn

  S        <- getSlot(m, "dataRep")$S
  n        <- getSlot(m, "dataRep")$n
    
  K            <- theta2K(m, theta, scale='original')
  dimnames(K)  <- dimnames(S)

  ##print("After score matching"); print(ellK(K,S,n-1))
  ## Ensure diagonals are positive
  if (min(diag(K))<0){
    vccI <- intRep(m)$vccI
    dS <- diag(S)
    for (i in 1:length(vccI)){
      cc <- as.numeric(vccI[[i]])
      dS[cc] <- mean(dS[cc])
    }
    diag(K) <- 1/dS
  }

  if (min(eigen(K)$values)<0)
    K     <- regularizeK(K)
  
  logL  <- ellK(K,S,n-1)
  ##print("After redoing diagonals"); print(logL)

  if (trace>=3)
    cat("..Score matching, logL:", logL, "\n")

  ans <- list(K=K, logL=logL, coef=theta)
  return(ans)
}


rconScoreTheta <- function(m){

  iR <- getSlot(m, "intRep")
  vccTerms <- iR$vcc
  eccTerms <- iR$ecc
  allTerms <- c(vccTerms, eccTerms)

  S    <- getSlot(m, "dataRep")$S
  n    <- getSlot(m, "dataRep")$n

  A <- matrix(0, ncol=length(allTerms), nrow=length(allTerms))
  B <- rep(0,length(allTerms))

  #print(A)
  
  for (u in 1:length(allTerms)){
    Ku <- allTerms[[u]]
    #print(Ku)
    bu     <- trA(Ku)
    B[u]   <- bu
    for (v in u:length(allTerms)){
      Kv <- allTerms[[v]]
      auv    <- trAWB(Ku, S, Kv)      
      A[v,u] <- A[u,v] <- auv
    }
  }
  #print(A); print(B)
  theta <- solve(A,B)  
  ##theta <- colSums(cholSolve(A)*B) ????  
  return(theta)
}

rcorScoreTheta <- function(m){


  theta        <- rconScoreTheta(m)
  vn           <- unlist(lapply(getcc(m),names))
  names(theta) <- vn

  S        <- getSlot(m, "dataRep")$S
  n        <- getSlot(m, "dataRep")$n

  oclass <- class(m)
  class(m) <- c("rcon","rcox")
  K            <- theta2K(m, theta, scale='original')
  print(K)
  dimnames(K)  <- dimnames(S)
  ##print(ellK(K,S,n-1))

  ## Ensure diagonals are positive
  vccI <- intRep(m)$vccI
  dS <- diag(S)
  for (i in 1:length(vccI)){
    cc <- as.numeric(vccI[[i]])
    dS[cc] <- mean(dS[cc])
  }
  diag(K) <- 1/dS
  
  K     <- regularizeK(K) 
  logL  <- ellK(K,S,n-1)
  ##print(ellK(K,S,n-1))
  
  class(m) <- oclass
  theta    <- K2theta(m, K, scale="original")

  vn       <- unlist(lapply(getcc(m),names))
  names(theta) <- vn
  return(theta)
}


  

#   ctrl   <- m$control
#   ctrl$vcov <- NULL
#   class(m) <- c("rcon","rcox")
#   KS       <- rconScoreMatch(m, control=ctrl)$K
#   class(m) <- oclass
#   K.curr   <- findKinModel(m, KS=KS, type='rcor', regularize=TRUE)

#   theta    <- K2theta(m, K.curr, scale="original")
#   vn       <- unlist(lapply(getcc(m),names))
#   names(theta) <- vn
#   return(theta)

#   if (!is.null(control$vcov)){
#     V <- calculateVCOV(m, K=K, vcov=control$vcov, nboot=control$nboot)
#   }




#   print(K)  
#   K <<- K
#   S <<- S
#   print(S)
#   K <<- K; S<<-S

#  print(ellK(K,S,n-1))
  
#  if (min(eigen(K)$values)<0){
#    cat("Score matching gave n.d. K; regularizing to make K p.d.\n")
#    K <- regularizeK(K)
#  }   
