
##############################################
trAWB <- function(A,W,B){
                                        # storage.mode(A) <- storage.mode(B) <- "double"
  .Call("trawb", A, W, B, PACKAGE="gRc")
}
##############################################


##############################################
 trAWBW <- function(A,W,B){
                                        # storage.mode(A) <- storage.mode(B) <- "double"
   .Call("trawbw", A, W, B, PACKAGE="gRc")
 }
##############################################


##############################################
trA <- function(A){
  if (ncol(A)==1)
    return(nrow(A))
  else
    return(0)
}
##############################################


##############################################
trAW <- function(A, W){
  
  if (ncol(A)==2)
    return(2*sum(W[A]))
  else{
    ##stop()
    return(sum(diag(W)[A]))
  }
}
##############################################


##############################################
trAWBV <- function(A,W,B,V){
  
  ncolA <- ncol(A)
  ncolB <- ncol(B)

  J     <- extractIndexMatrix(A,B)
  #print(J)

# if (ncolA==1 && ncolB==1){
#     JJ  <- cbind(rep(A,each=nrow(B)),rep(B,nrow(A)))
#     ans <- sum(W[JJ]*V[JJ[,2:1]])
#     print(JJ)
#     return(ans)
#   }

  scale <- ncolA*ncolB/2
  sum(W[J[,1:2]] * V[J[,3:4]]) * scale
}
##############################################



##############################################

extractIndexMatrix <- function(x, y) {
  storage.mode(x) <- storage.mode(y) <- "double"
  
  .Call("eim", x, y, PACKAGE="gRc")
}

##############################################






# .zeroAB <- function(A,B){
# #  if (ncol(A)>1)
# #    return(TRUE)
# #  print(A); print(B)
#   dim(B) <- NULL
#   d <- rep(A, each=length(B)) - rep(B, length(A))
#   !any(d==0)
# }


# trAWB <- function(A, W, B){
#   UseMethod("trAWB")
# }

# trAWB.list <- function(A, W, B){
#   A <- do.call("rbind", A)
#   B <- do.call("rbind", B)
#   trAWB.matrix(A,W,B)
# }

# trAWB.matrix <- function(A,W,B){
  
#   ncolA <- ncol(A)
#   ncolB <- ncol(B)

# #   if (ncolA==1){
# #     ii <- zeroAB(A,B)
# #     print(ii)
# #   }
  
#   scale <- ncolA*ncolB/4
#   J     <- extractIndexMatrix(A,B)
  
#   Ui  <- as.numeric(J[,1]==J[,2])
#   Zi  <- as.numeric(J[,3]==J[,4])
#   ans <- sum(c(W[J[,1:2]] * Zi, W[J[,3:4]] * Ui))* scale  
#   ans
# }


# .extractIndexMatrix <- function(A,B){
#   ncolA <- ncol(A)
#   ncolB <- ncol(B)

#   B       <- t(B)
#   b2      <- B;
#   dim(b2) <- NULL
  
#   if (ncolB==2){
#     b22 <- B[2:1,];
#     dim(b22) <- NULL
#   } else {
#     b22 <- rep(b2, each=2)
#   }
  
#   lenb <- length(b22)  
#   J <- matrix(0, nc=4, nr=nrow(A)*lenb)

#   if (ncolA==2){    
#     J[,1] <- rep(A[,1], each=lenb)
#     J[,3] <- rep(A[,2], each=lenb)    
#   } else {
#     J[,1] <- J[,3] <- rep(A[,1], each=lenb)
#   }
  
#   if (nrow(B)==2) {
#     J[,2] <- b22
#     J[,4] <- b2         
#   } else {
#     J[,2] <- J[,4] <- b2
#   }
#   return(J)
# }



# trAWBW <- function(A, W, B){
#   UseMethod("trAWBW")
# }

# trAWBW.list <- function(A, W, B){
#   A <- do.call("rbind", A)
#   B <- do.call("rbind", B)
#   trAWBW.matrix(A,W,B)
# }

# trAWBW.matrix <- function(A, W, B){
  
#   ncolA <- ncol(A)
#   ncolB <- ncol(B)

#   if (ncolA==1 && ncolB==1){
#     #print("HHHHH")
#     JJ <- cbind(rep(A,each=nrow(B)),rep(B,nrow(A)))
#     return(sum(W[JJ]^2))
#   }

#   scale <- ncolA*ncolB/2
#   J     <- extractIndexMatrix(A,B)

#   sum(W[J[,1:2]] * W[J[,3:4]]) * scale
# }


# trA2 <- function(A){
#   UseMethod("trA2")
# }

# trA2.list <- function(A){
#   length(which(sapply(A,length)==1))
# }


# trAW2 <- function(A, W){
#   UseMethod("trAW2")
# }

# trAW2.list <- function(A, W){
#   A <- do.call("rbind", A)
#   trAW2.matrix(A,W)
# }



# trAWBV2 <- function(A, W, B, V){
#   UseMethod("trAWBV2")
# }

# trAWBV2.list <- function(A, W, B, V){
#   A <- do.call("rbind", A)
#   B <- do.call("rbind", B)
#   trAWBV2.matrix(A,W,B,V)
# }
