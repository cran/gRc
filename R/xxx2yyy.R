##
## Conversions of the type xxx2yyy
##

formula2names <- function(x)
  UseMethod("formula2names")
formula2names.NULL     <- function(x){}
formula2names.list     <- function(x){
  ##lapply(x,formula2names)
  
  lapply(x, function(x2){
    if (class(x2)=="list") {
      x2
    } else {
      formula2names(x2)
    }
  })
}

formula2names.default     <- function(x){x}
  
formula2names.formula  <- function(x){  
  mf <- paste(deparse(x[[2]]),collapse="")
  mf <- gsub(" +","",mf)
  strsplit(unlist(strsplit(mf,"\\+")),":")
}



formula2string <- function(x)
  UseMethod("formula2string")

formula2string.list <- function(x)
  lapply(x, formula2string)
  
formula2string.default <- function(x){
  paste(paste(x),collapse="")  
}

names2formula           <- function(x) 
  UseMethod("names2formula")
names2formula.list      <- function(x) 
  lapply(x, names2formula.default)
names2formula.default <- function(x){
  if (is.null(x))
    return(NULL)
  as.formula(paste("~",paste(sapply(x, paste, collapse=':'),collapse='+')))
}


## Fragile........
cc2formula <- function(cc){
  v<-as.formula(paste("~", paste(unlist(lapply(lapply(cc, unlist), paste, collapse=":")),collapse="+")))
  v
}


names2indices <- function(x,vn){
  ##print("names2indices"); print(x); print(vn)
  if(!is.null(x)){
    z <-listOrder(getIndex(x,vn))
    return( if (length(z)>0) z)
  }}


getIndex           <- function(x,vn)UseMethod('getIndex')
getIndex.list      <- function(x,vn){ lapply(x, getIndex, vn) }
getIndex.character <- function(x,vn){ match(x,vn) }
##getIndex.character <- function(x,vn){ vn[x] }
getIndex.default   <- function(x,vn){ x }


names2pairs <- function(x,y=x){
  val <- as.list(rep(NA, length(x)))
  for (i in 1:length(x)){
    xi  <- x[i]
    y2  <- setdiff(y,xi)
    val[[i]] <- lapply(lapply(y2,c,xi),sort)
  }

  val <- unique(unlist(val, recursive=FALSE))
  val
}



ecc2edges <- function(x){
  if (length(x)==0)
    return(NULL)
  x2<- do.call("rbind", (lapply(x, function(d) do.call("rbind", d))))
  mapply(function(a,b)c(a,b), x2[,1],x2[,2], SIMPLIFY=FALSE,USE.NAMES=FALSE)
}
