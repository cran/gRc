#######################################################################
#' @title Update an RCOX model
#' @description update will update and (by default) re-fit an RCOX
#'   model.  It does this by extracting the call stored in the
#'   object, updating the call and (by default) evaluating that
#'   call. Sometimes it is useful to call update with only one
#'   argument, for example if the data frame has been corrected
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name update-doc
#######################################################################
#'
#' @aliases update.rcox
#' 
#' @param object An RCOX model, an object of class RCOX
#' @param vcc Specification of the vertex colour classes in the model
#' @param ecc Specification of the edge colour classes in the model
#' @param splitvcc Existing vertex colour class to be split
#' @param splitecc Existing edge colour class to be split
#' @param joinvcc Existing vertex colour classes to be joined
#' @param joinecc Existing vertex colour classes to be joined
#' @param addecc New edge colour classes to be added
#' @param dropecc Existing vertex color classes to be dropped
#'   (deleted)
#' @param Kstart A start value for K
#' @param fit Should the updated model be fitted.
#' @param control A list of control parameters.
#' @param trace For debugging purposes
#' @param ... Additional arguments, currently not used.
#' @return A new model object of class 'rcox'.
#'
#' @section Warning: Only one of the arguments pertaining to edge
#'   colour classes (i.e. ecc, splitecc, joinecc, dropecc, addecc)
#'   should be applied at the time. Likewise for the arguments
#'   pertaining to the vertex colour classes.
#' 
#' The result will otherwise be highly unpredictable and is likely to cause an
#' error.

#' @seealso \code{\link{rcox}}
#' @keywords models
#' @examples
#' 
#' 
#' data(math)
#' gm  = ~al:an:st
#' vcc = list(~me+st, ~ve+an, ~al)
#' ecc = list(~me:ve+me:al, ~ve:al+al:st)
#' 
#' m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching', trace=0)
#' 
#' update(m1, joinvcc=list(~me+st, ~ve+an))
#' update(m1, joinecc=list(~al:an, ~an:st))
#' 
#' update(m1, splitvcc=~ve+an)
#' update(m1, splitecc=~me:ve+me:al)
#' 
#' 
#' update(m1, dropecc=list(~me:st+st:an,~al:an,~st:al))
#' update(m1, addecc=list(~an:me+st:ve))
#' 
NULL

#' @export
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


  .intersectListList <- function(addccV, oldccV){
    x <- lapply(addccV,
                function(ccnew){
                  lapply(oldccV,
                         function(ccold){ intersect(ccnew,ccold) })})
    x
  }

  nodes <- object$nodes
  
  if (!is.null(joinvcc)){
    oldcc   <- getSlot(object,"vcc")
    joinvcc   <- .ccl2names(joinvcc, oldcc)
    new.ccl   <- join_cc(joinvcc, oldcc)
    new.ccl   <- .addccnames(new.ccl, type="vcc")
    vcc       <- new.ccl
    if (trace>=1)cat(".joining vcc:", toLisp(joinvcc),"\n")
  }

  if (!is.null(joinecc)){
    oldcc   <- getSlot(object,"ecc")
    joinecc   <- .ccl2names(joinecc, oldcc)
    new.ccl   <- join_cc(joinecc, oldcc)
    new.ccl   <- .addccnames(new.ccl, type="ecc")
    ecc       <- new.ccl
    if (trace>=1)cat(".joining ecc:", toLisp(joinecc),"\n")
  }
  
  if (!is.null(splitvcc)){
    oldcc    <- getSlot(object,"vcc")
    splitvcc   <- .ccl2names(splitvcc, oldcc)
    new.ccl    <- split_cc(splitvcc, oldcc)
    new.ccl    <- .addccnames(new.ccl, type="vcc")
    vcc        <- new.ccl
    if (trace>=1)cat(".splitting vcc:", toLisp(splitvcc),"\n")
  }
  
  if (!is.null(splitecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    splitecc   <- .ccl2names(splitecc, oldcc)    
    
    new.ccl    <- split_cc(splitecc, oldcc)
    new.ccl    <- .addccnames(new.ccl, type="ecc")
    ecc        <- new.ccl
    if (trace>=1)cat(".splitting ecc:", toLisp(splitecc),"\n")
  }

  if (!is.null(addecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    addecc  <- formula2names(addecc)
    #cat("addecc:\n"); print(addecc)
    #cat("oldecc:\n"); print(oldcc)
    if (length(oldcc)==0){
      ecc <- addecc
    } else {
      addccV <- indices2vectors(names2indices(formula2names(addecc), nodes, matrix=TRUE))
      x      <- .intersectListList(addccV, oldccV)
      if (length(unlist(x))>0)
        stop("Can not add ecc to model\n");
      ecc <- c(oldcc, addecc)
    } 
    ecc    <- .addccnames(ecc, type="ecc")
    #print(ecc)
    if (trace>=1)cat(".add ecc:", toLisp(addecc),"\n")
  }
  
  if (!is.null(dropecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    dropecc  <- formula2names(dropecc)

    dropccV <- indices2vectors(names2indices(formula2names(dropecc), nodes, matrix=TRUE))
    #print(dropccV)
    #print(oldccV)
    
                                        #    idx       <- sapply(dropecc, matchLL2, oldcc)
    idx <- unlistPrim(lapply(dropecc, matchLL2, oldcc))
    idx       <- which(!is.na(idx))
    dropecc   <- dropecc[idx]    
    #idx       <- sapply(dropecc, matchLL2, oldcc)
    idx       <- unlistPrim(lapply(dropecc, matchLL2, oldcc))
    
    new.ccl   <- oldcc[-idx]    
    ecc       <- new.ccl

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
    object$fitInfo <- fit(object, trace=trace, returnModel=FALSE)
  else
    object$fitInfo <- NULL
  return(object)
}


split_cc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2, old.ccl)
  if ((length(idx)>1) || (!is.na(idx))){
    old.ccl <- old.ccl[-idx]
  }   
  new.cc    <- lapply(unlist(cc, recursive=FALSE),list)
  new.ccl   <- unionL2L2(old.ccl,new.cc)
  new.ccl
}

join_cc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2, old.ccl)
  if ((length(idx)>1) || (!is.na(idx))){
    old.ccl <- old.ccl[-idx]
  }   
  new.cc         <- list(unlist(cc, recursive=FALSE))
  new.ccl        <- unionL2L2(old.ccl, new.cc)
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



