######################################################################
#' @title Conversions of the type xxx2yyy
#' @description Conversions of the type xxx2yyy.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name xxx2yyy
######################################################################
#'
#' @aliases formula2names formula2names.default formula2names.list
#'   formula2names.formula formula2names.NULL names2indices
#'   names2formula names2formula.list names2formula.default
#'   indices2vectors
#'   formula2string formula2string.list formula2string.default
#'   getIndex getIndex.default getIndex.list getIndex.character
#'    ecc2edges cc2formula
NULL


#' @export
formula2string <- function(x){
  UseMethod("formula2string")
}

#' @export
formula2string.list <- function(x){
  lapply(x, formula2string)
}

#' @export
formula2string.default <- function(x){
  paste(paste(x),collapse="")  
}

getIndex           <- function(x, vn){
  UseMethod('getIndex')
}

getIndex.default   <- function(x, vn){
  x
}

getIndex.list      <- function(x, vn){
  lapply(x, function(a) getIndex(a, vn))
}

getIndex.character <- function(x, vn){
  match(x,vn)
}


## formula2names
##

#' @export
formula2names <- function(x){
  UseMethod("formula2names")
}

#' @export
formula2names.default     <- function(x){
  x
}

#' @export
formula2names.NULL     <- function(x){
  
}
#' @export
formula2names.list     <- function(x){
  lapply(x, function(x2){
    if (class(x2)=="list") {
      x2
    } else {
      formula2names(x2)
    }
  })
}

#' @export
formula2names.formula  <- function(x){  
  mf <- paste(deparse(x[[2]]), collapse="")
  mf <- gsub(" +", "", mf)
  strsplit(unlist(strsplit(mf, "\\+")), ":")
}


## names2indices
##

names2indices <- function(x,vn,matrix=TRUE){
  if(!is.null(x)){
    z <-listOrder(getIndex(x,vn))
    names(z)<- NULL
    ans <- if (length(z)>0) z
    if (length(ans)>0 & matrix==TRUE)
      ans <- lapply(ans, function(xx) do.call(rbind, xx))
    return(ans)
  }}


## names2formula
##
names2formula           <- function(x) 
  UseMethod("names2formula")
names2formula.list      <- function(x) 
  lapply(x, names2formula.default)
names2formula.default <- function(x){
  if (is.null(x))
    return(NULL)
  as.formula(paste("~",paste(sapply(x, paste, collapse=':'),collapse='+')))
}


## indices2vectors
##
indices2vectors <- function(x){
  if (length(x)==0)
    return(NULL)
  if (ncol(x[[1]])==1){
    lapply(x, as.numeric)
  } else {
    lapply(x, function(xx) xx[,1]*100000 + xx[,2])  
  }
}





ecc2edges <- function(x){
  if (length(x)==0)
    return(NULL)
  x2<- do.call("rbind", (lapply(x, function(d) do.call("rbind", d))))
  mapply(function(a,b)c(a,b), x2[,1],x2[,2], SIMPLIFY=FALSE,USE.NAMES=FALSE)
}


## Fragile........
cc2formula <- function(cc){
  v<-as.formula(paste("~", paste(unlist(lapply(lapply(cc, unlist), paste, collapse=":")),collapse="+")))
  v
}
