## Based on a start value KS an estimate K in the model is returned.
##
## Entries in K which should be zero according to the model are so.
## Non-zero elements in K are obtained by averaging relevant entries in KS
##

findKinModel <- function(m, KS, type="rcon", regularize=TRUE)
  UseMethod("findKinModel")


## MATRIX VERSION
findKinModel.rcox <- function(m, KS, type="rcon", regularize=TRUE){
  if (is.null(KS))
    return(NULL)
  
  switch(type,
         'rcor'={
           ##C <- cov2cor(KS)
           ##print(C)
           C <- findKinModel(m,cov2cor(KS),   type='rcon')
           A <- findKinModel(m,diag(diag(KS)),type='rcon')
           a <- sqrt(diag(A))
           KKmod <- a * t(a*C) ## Short for KKmod <- A%*%C%*%A
           KKmod           
         },
         
         'rcon'={
           ##print("RCON")
           VCC <- intRep(m,"vccI")
           ECC <- intRep(m,"eccI")
           for (i in seq(along=ECC)){
             x      <- ECC[[i]]
             if (nrow(x)>1){             
               x2    <- x[,2:1,drop=FALSE]
               KS[x] <- KS[x2] <- mean(KS[x])
             }
           }
           
           for (i in seq(along=VCC)){
             x           <-unlist(VCC[[i]])
             if (nrow(x)>1){
               x           <- as.numeric(x)
               diag(KS)[x] <- mean(diag(KS)[x])
             }
           }
           if (regularize && (min(eigen(KS)$values) < 0)){             
             KKmod <- regularizeK(KS)
           } else {
             KKmod <- KS
           }
         }
         )
  return(KKmod)
}

## K is made positive definite
##

regularizeK <- function(K){
  Kdiag <- abs(diag(diag(K)))           
  Krest <- K-Kdiag
  alpha <- 0.95
  repeat{
    Kmod  <- Kdiag + alpha*Krest
    
    if (min(eigen(Kmod)$values) > 0){
      Kmod <- Kdiag + 0.95*alpha*Krest ## Be on the safe side...
      cat("The end", alpha,"\n")
      break()
    }
    alpha <- alpha-0.05    
  }
  Kmod
}


# regularizeK2 <- function(K){
#   Kdiag <- abs(diag(diag(K)))           
#   Krest <- K-Kdiag
#   alpha <- 1.1
#   repeat{
#     Kmod  <- alpha*Kdiag + Krest 
# #    print(alpha)
# #    print(eigen(Kmod)$values)
#     alpha <- alpha+0.1
#     if (min(eigen(Kmod)$values) > 0){
#       Kmod <- alpha * Kdiag + Krest ## Be on the safe side...
#       ##cat("The end", alpha,"\n")
#       break()
#     }
#   }
#   Kmod
# }












# findKinModel.rcox <- function(m, KS, type="rcon", regularize=TRUE){
#   if (is.null(KS))
#     return(NULL)
  
#   switch(type,
#          'rcor'={
#            C <- findKinModel(m,cov2cor(KS),   type='rcon')
#            A <- findKinModel(m,diag(diag(KS)),type='rcon')
#            a <- sqrt(diag(A))
#            #A <- diag(sqrt(diag(A)))
#            KKmod <- a * t(a*C) ## Short for KKmod <- A%*%C%*%A
#            KKmod
#          },
         
#          'rcon'={
#            KK  <- matrix(0, ncol=ncol(KS), nrow=nrow(KS))
#            dimnames(KK) <- dimnames(KS)##print(dimnames(dataRep(m,"S")))
#            VCC <- intRep(m,"vccI")
#            ECC <- intRep(m,"eccI")
           
#            for (i in seq(along=ECC)){
#              x      <- ECC[[i]]
#              id     <- matrix(unlist(x),ncol=2,byrow=TRUE)
#              KK[id] <- mean(KS[id])
#            }
#            KK <- KK + t(KK)
#            for (i in seq(along=VCC)){
#              x           <-unlist(VCC[[i]])
#              diag(KK)[x] <- mean(diag(KS)[x])
#            }
           
#            if (regularize && (min(eigen(KK)$values) < 0)){             
#              KKmod <- regularizeK(KK)
#            } else {
#              KKmod <- KK
#            }
#          }
#        )
#   return(KKmod)
# }








##findKinModelPrimitive <- function(m,KS,type='rcon'){
#              KKdiag <- diag(diag(KK))           
#              ## Regularize K
#              KKrest <- KK-KKdiag
#              alpha  <- 0.95
#              repeat{
#                ##print(alpha)
#                KKmod <- KKdiag + alpha*KKrest
#                alpha <- alpha-0.05
#                if (min(eigen(KKmod)$values) > 0){
#                  KKmod <- KKdiag + alpha*KKrest ## Be on the safe side...
#                  break()
#                }
#              }
