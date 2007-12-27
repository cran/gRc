
update.rcox <- function(object,
                        vcc       = NULL,
                        ecc       = NULL,
                        splitecc  = NULL,
                        splitvcc  = NULL,
                        joinvcc   = NULL,
                        joinecc   = NULL,
                        addecc    = NULL,
                        dropecc   = NULL,
                        Kstart    = NULL,
                        fit       = TRUE,
                        control   = NULL,
                        trace=object$trace,
                        ...){

  if (!is.null(joinvcc)){
    old.ccl   <- getSlot(object,"vcc")
    joinvcc   <- .ccl2names(joinvcc, old.ccl)
    new.ccl   <- .joincc(joinvcc, old.ccl)
    new.ccl   <- .addccnames(new.ccl, type="vcc")
    vcc       <- new.ccl
    if (trace>=1)cat(".joining vcc:", toLisp(joinvcc),"\n")
  }

  if (!is.null(joinecc)){
    old.ccl   <- getSlot(object,"ecc")
    joinecc   <- .ccl2names(joinecc, old.ccl)
    new.ccl   <- .joincc(joinecc, old.ccl)
    new.ccl   <- .addccnames(new.ccl, type="ecc")
    ecc       <- new.ccl
    if (trace>=1)cat(".joining ecc:", toLisp(joinecc),"\n")
  }
  
  if (!is.null(splitvcc)){
    old.ccl    <- getSlot(object,"vcc")
    splitvcc   <- .ccl2names(splitvcc, old.ccl)
    new.ccl    <- .splitcc(splitvcc, old.ccl)
    new.ccl    <- .addccnames(new.ccl, type="vcc")
    vcc        <- new.ccl
    if (trace>=1)cat(".splitting vcc:", toLisp(splitvcc),"\n")
  }
  
  if (!is.null(splitecc)){
    old.ccl    <- getSlot(object,"ecc")
    splitecc   <- .ccl2names(splitecc, old.ccl)
    new.ccl    <- .splitcc(splitecc, old.ccl)
    new.ccl    <- .addccnames(new.ccl, type="ecc")
    ecc        <- new.ccl
    if (trace>=1)cat(".splitting ecc:", toLisp(splitecc),"\n")
  }
  
  if (!is.null(addecc)){
    old.ccl    <- getSlot(object,"ecc")
    if (length(old.ccl)>0){
        
      if (is.L(addecc)){
        addecc <- lapply(addecc, list)
        class(addecc)<- c("colourClass", "list")
      }
      
      addecc  <- .ccl2names(addecc, old.ccl)
      idx     <- sapply(addecc, function(e1)
                        any(is.na(sapply(e1, matchLL2, old.ccl))) )
      addecc  <- addecc[idx]
      new.ccl <- unionL2L2(addecc, old.ccl)
      ecc      <- new.ccl
    } else {
      ecc <- addecc
    }    
    ecc    <- .addccnames(ecc, type="ecc")
    if (trace>=1)cat(".add ecc:", toLisp(addecc),"\n")
  }
    
  if (!is.null(dropecc)){
    old.ccl    <- getSlot(object,"ecc")
    dropecc   <- .ccl2names(dropecc, old.ccl)
    idx <- sapply(dropecc, matchLL2, old.ccl)
    idx <- which(!is.na(idx))
    dropecc <- dropecc[idx]    
    idx <- sapply(dropecc, matchLL2, old.ccl)
    new.ccl <- old.ccl[-idx]    
    ecc        <- new.ccl

    if (trace>=1)cat(".drop ecc:", toLisp(dropecc),"\n")
  }

  if (!is.null(vcc)){
    vcc <- .addccnames(vcc, "vcc")
    object$vcc <- vcc
  }

  if (!is.null(ecc)){
    ecc <- .addccnames(ecc, "ecc")
    object$ecc <- ecc
  }
  
  if (trace>=3)
    cat("...(update) Updating internal representation of model object...\n")
  intRep <- .buildInternalRepresentation(vccN      = object$vcc,
                                         eccN      = object$ecc,
                                         dataNames = object$dataRep$dataNames,
                                         trace     = 2)
  object$intRep <- intRep
  object$Kstart <- Kstart
  
  if (!is.null(control)){
    object$control[(namc <- names(control))] <- control
  }
  
  if (fit)# & !is.null(object$fitInfo))
    object$fitInfo <- .fitit(object, trace=trace)
  else
    object$fitInfo <- NULL
  return(object)
}







.splitcc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2,old.ccl)
  new.cc    <- lapply(unlist(cc, recursive=FALSE),list)
  new.ccl   <- unionL2L2(old.ccl[-idx],new.cc)
  new.ccl
}

.joincc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2,old.ccl)
  new.cc    <- list(unlist(cc, recursive=FALSE))
  new.ccl   <- unionL2L2(old.ccl[-idx],new.cc)
  class(new.ccl) <- union(class(old.ccl),class(cc))
  new.ccl
}

.ccl2names <- function(x,y){
  if (class(x)[1]=="formula")
    x <- list(x)
  cc <- lapply(x, function(xx){
    switch(class(xx),
      "list"    = {xx},
      "numeric" = {y[[xx]]},
      "formula" = {formula2names(xx)}
      )
  })
  cc
}



