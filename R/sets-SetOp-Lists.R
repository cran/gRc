
# ## Throw out elements which are contained in other elements
# ##

# maximalSet <- function(set){
#   set <- unique(cardOrder(lapply(set,sort)))
  
#   idxv <- rep(TRUE, length(set))
#   for (i in 1:length(set)){
#     el <- set[[i]]
#     r  <- set[-c(1:i)]
#     if(length(r)>1){  
#       for (k in 2:length(r)){
#         if (is.subset(el, r[[k]])){
#           idxv[i] <- FALSE
#           #print(r[[k]])
#           break()
#         }  
#       }
#     }
#   }
#   set[idxv]
# }


# ## Order set by the cardinality of the elements
# ##
# cardOrder <- function(x){
#   len <- unlist(lapply(x,length))
#   unlen <- sort(unique(len))
#   x2 <- NULL
#   for (i in seq(along=unlen)){
#     x2  <- c(x2, x[len==unlen[i]])
#   }
#   x2
# }



# is.subset <- function(x,y){
#   setequal(intersect(x,y),x)
# }



## Sorting lists
##
listOrder         <- function(x) UseMethod('listOrder')

listOrder.numeric <- function(x){ #print("numeric"); print(x); 
  cl <- class(x)
  x <- x[order(x)] 
  class(x)<- cl
  x
  }

listOrder.list    <- function(x){ #print("list   "); print(x); 
  cl <- class(x)
  x <- lapply(x,function(v)listOrder(v)) 
  class(x)<- cl
  x
  }

listOrder.default <- function(x){ #print("default"); print(x); 
  x 
  }
