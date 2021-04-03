###################################################################### 
#' @title Get slots from RCOX model object.
#' @description Get slots from RCOX model object.
#' @author  Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name get-slot
###################################################################### 
#'
#' @param object RCOX model object.
#' @param slot slot.
#' @param type Type of colour class.
#' @param complement If FALSE, the edges of the model is returned. If
#'   TRUE, the edges not in the model is returned.
#' @examples
#'
#' 
#' data(math)
#' gm  = ~al:an:st
#' vcc = list(~me+st, ~ve+an, ~al)
#' ecc = list(~me:ve+me:al, ~ve:al+al:st)
#'
#' m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math)
#' getecc(m1)
#'
#' getSlot(m1,"type")
#' fitInfo(m1)
#' fitInfo(m1,"K")
NULL



#' @rdname get-slot
#' @export
getSlot<-function(object, slot){
  object[[slot]]
}

#' @rdname get-slot
#' @export
dimension  <- function(object){
  length(c(getSlot(object,'vcc'), getSlot(object,'ecc')))
}

#' @rdname get-slot
#' @export
logL <- function(object){
  getSlot(object,'fitInfo')$logL
}


#' @rdname get-slot
#' @export
getSlot<-function(object, slot){
  if(is.null(slot))
    return(object)
  
  return(object[[slot]])
}

#' @rdname get-slot
#' @export
dataRep <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"dataRep"))
  getSlot(object,"dataRep")[[slot]]
}

#' @rdname get-slot
#' @export
intRep <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"intRep"))
  getSlot(object,"intRep")[[slot]]
}

#' @rdname get-slot
#' @export
fitInfo <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"fitInfo"))
  getSlot(object,"fitInfo")[[slot]]
}


#' @rdname get-slot
#' @export
getcc <- function(object,type){
  if (missing(type))
    list(vcc=object$vcc, ecc=object$ecc)
  else {
    switch(type,
           "ecc"={object$ecc},
           "vcc"={object$vcc})
  }
}

#' @rdname get-slot
#' @export
getecc <- function(object){
  object$ecc
}

#' @rdname get-slot
#' @export
getvcc <- function(object){
  object$vcc
}

#' @rdname get-slot
#' @export
getedges <- function(object,complement=FALSE){
  ans <- ecc2edges(getecc(object))
  if (complement){
    eAll <- names2pairs(getSlot(object,"nodes"))
    ans  <- setdiffLL(eAll, ans)
  }
  ans
}

print.colourClass <- function(x,...){
  xf <- names2formula(x)
  xs <- formula2string(xf)
  mapply(function(n,xxx) cat(n,xxx,"\n"), names(xs),xs)
  return(invisible(x))
}

#' @export
coef.rcox <- function(object, ...){
  co  <- fitInfo(object,"coef")
  co
}

tocc <- function(v){
  if(length(v)==0)
    return(NULL)
  as.cclist(
  lapply(v, function(x) {
    if (length(x)==1)
      as.cc(as.atom(x))
    else
      as.cc(lapply(x, as.atom))
    })
  )
}


cc2str <- function(cc){
  paste(sapply(cc, toLisp), collapse='')  
}

.addccnames <- function(x, type){
  if (length(x)){
    names(x) <- paste(type,paste(1:length(x)),sep="")
    class(x) <- c("colourClass","list")
    x
  } else {
    NULL
  }
}



cholSolve <- function(ma)
  chol2inv(chol(  ma  ))


#' @export
logLik.rcox <- function(object, ...){
  out <- fitInfo(object)$logL
  attr(out, "df") <- NA
  attr(out, "dimension") <- dimension(object)
  class(out) <- "logLik"
  out
}
