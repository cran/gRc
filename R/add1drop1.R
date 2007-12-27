## Calculate test statistic for addition of each atomic ecc not already in the model
## Note: Based on fitting a new model for each new ecc, and may hence be slow...
##



add1.rcox <- function(object, scope, details=1, trace=0, ...){
  n       <- dataRep(object,"n")
  if (missing(scope)){
    eNew <- getedges(object, complement=TRUE)
  } else{
    eNew <- .addccnames(formula2names(scope),type="ecc")
  }
  
  if (is.L(eNew)){
    eNew <- lapply(eNew, list)
    class(eNew)<- c("colourClass", "list")
  }
  
  if (length(eNew)==0)
    return(NULL)

  object$control$vcov <- NULL

  res <- rep(NA, length(eNew))
  for (i in 1:length(eNew)){
    e      <- eNew[i]
    ##mtmp   <- update(object, addecc=list(e))
    mtmp   <- update(object, addecc=e)
    dev    <- -2*(logL(object)-logL(mtmp))
    res[i] <- dev
  }

  eNew <- .addccnames(eNew,"ecc")
  ##ans <- data.frame(cc=sapply(tocc(eNew), cc2str), X2=res, df=1)
  ans <- data.frame(cc=names(eNew), X2=res, df=1)
  ans <- .addStat(ans, n=n, direction="add")
  
  attr(ans,"ccterms") <- eNew
  ans2 <- structure(list(tab=ans, cc=eNew, details=details), class=c("statTable","data.frame"))
  ans2
}

print.statTable <- function(x,...){
  print(x$tab)
  if (x$details>=1){
    if (!is.null(x$cc)){
      cat("\ncc:\n")
      print(x$cc)
    }
    if (!is.null(x$cc1)){
      cat("\ncc1:\n")
      print(x$cc1)
      cat("cc2:\n")
      print(x$cc2)
    }
  }
  cat(paste("\nAvailable components:", paste(setdiff(names(x),"details"),collapse=' ')),"\n")
  return(invisible(x))
}
  

## Calculate test statistic for deletion of each ecc in the model
## Note: By default based on Wald statistics from existing model
##
drop1.rcox <- function(object, scope, details=1, trace=0, stat="wald", ...){
  stat <- match.arg(stat,c("wald","dev"))
  n   <- dataRep(object,"n")

  if (missing(scope))
    ec  <- getSlot(object,'ecc')
  else
    ec  <- .addccnames(formula2names(scope),type="ecc")
  
  if (details>=1){
    cat("Statistic:", stat, "\n")
  }
  
  if (length(ec)==0)
    return(NULL)

  res <- rep(NA, length(ec))

  if (stat=="wald"){
    V   <- vcov(object)
    b   <- coef(object)
    ofs <- length(getSlot(object,"vcc"))
  
    ccidx <- sapply(ec, matchLL2, ec)
    lcc   <- length(ccidx)
    for (i in ccidx){    
      i2   <- i + ofs
      dev  <- b[i2]^2/V[i2,i2]
      res[i] <- dev
    }
  } else {
    for (i in 1:length(ec)){
      e      <- ec[[i]]
      mtmp   <- update(object, dropecc=list(e))
      dev    <- 2*(logL(object)-logL(mtmp))
      res[i] <- dev
    }
  }
  
  ans <- data.frame(cc=names(ec), X2=res, df=1)
  ans <- .addStat(ans, n=n, direction="drop")
  attr(ans,"ccterms") <- ec
  ans
  ec <- .addccnames(ec,"ecc")
  ans2 <- structure(list(tab=ans, cc=ec, details=details), class=c("statTable","data.frame"))
  ans2
}




