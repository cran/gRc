## ######################################################################
## #' @title Set operations
## #' @description These functions are not intended to be called by the
## #'   user.
## #' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
## #' @name set-op
## ######################################################################
## #'
## #' @aliases cardOrder listOrder listOrder.list listOrder.numeric
## #'   listOrder.default
## #' 
## NULL

## cardOrder <- function(xlist){
##   x <- xlist
##   len <- unlist(lapply(x,length))
##   unlen <- sort(unique(len))
##   x2 <- NULL
##   for (i in seq(along=unlen)){
##     x2  <- c(x2, x[len==unlen[i]])
##   }
##   x2
## }


## Sorting lists
##

listOrder         <- function(x) UseMethod('listOrder')

#' @export
listOrder.numeric <- function(x){ #print("numeric"); print(x); 
  cl <- class(x)
  x <- x[order(x)] 
  class(x)<- cl
  x
}

#' @export
listOrder.list    <- function(x){ #print("list   "); print(x); 
  cl <- class(x)
  
  ##  print("----------------");print(x)
  ##x <- lapply(x,function(v)listOrder(v))
  x <- lapply(x,
    function(v){
      
      #print(class(v))
      #listOrder(v)
      if (is(v, "list"))
        listOrder.list(v)
      else {
        cl <- class(v)
        x <- v[order(v)]
        class(v)<- cl
        v
      }
    }
  ) 
  class(x)<- cl
  x
}

## listOrder.default <- function(x){ #print("default"); print(x); 
##   x 
## }


#' @export
listOrder.cc    <- listOrder.list

## listOrder.atom  <- listOrder.numeric

