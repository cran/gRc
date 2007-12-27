
rcox <- function(gm=NULL, vcc=NULL, ecc=NULL, 
                 type   = c('rcon','rcor'),
                 method = c("scoring","ipm", "matching", "user"),
                 fit=TRUE, data=NULL, S=NULL, n=NULL,
                 Kstart=NULL,
                 control=list(),
                 details=1,
                 trace=0){

  type           <- match.arg(type)
  method         <- match.arg(method)
  modelSpec      <- list(gm=gm, vcc=vcc, ecc=ecc)

  con <- list(
              vccfit   = TRUE,
              eccfit   = TRUE,
              vcov     = {if(method=="scoring" || method=="ipm") "inf" else NULL},
              nboot    = 100,
              maxouter = 25,
              maxinner = 25,
              logL     = FALSE,
              logLeps  = 1e-6,
              deltaeps = 1e-6,
              short    = FALSE
              )

  con[(namc <- names(control))] <- control

  if (trace>=2) cat("..Building standard representation\n")
  stdRep   <- .buildStandardRepresentation (gm, vcc, ecc, trace)
  vccN     <- stdRep$vccN
  eccN     <- stdRep$eccN
    
  vccN     <- .addccnames(vccN, type="vcc")
  eccN     <- .addccnames(eccN, type="ecc")
  usedVars <- unique(unlist(c(vccN, eccN)))

  if (trace>=2) cat("..Building data representation\n")
  dataRep <- .buildDataRepresentation (data, S, n, usedVars, type, trace)

  nodes <- dataRep$nodes
  
  if (trace>=2) cat("..Building internal representation\n")
  intRep <- .buildInternalRepresentation(vccN=vccN, eccN=eccN,
                                         dataNames=dataRep$dataNames, trace)
  ans <- structure(list(##call    = match.call(),
                        vcc     = vccN,                        
                        ecc     = eccN,
                        nodes   = nodes,
                        intRep  = intRep,
                        dataRep = dataRep,
                        Kstart  = Kstart,
                        type    = type,
                        method  = method,
                        trace   = trace,
                        control = con
                        ),
                   class=c(type, "rcox"))

  if (fit){
    ans$fitInfo <- .fitit(ans, method=method, trace=trace)
  }
  return(ans)

}

print.rcox <- function(x, ...){
  cat(toupper(getSlot(x,"type")), "model: ")
  if (!is.null(x$fitInfo)){
    cat("logL=", x$fitInfo$logL, " dimension=", dimension(x),
        " method=", x$method, " time=", fitInfo(x,"time"),sep='')
  }
  cat("\n")

  if (!(x$control$short)){
    xcc <-getvcc(x)
    xf  <- names2formula(xcc)
    xs  <-formula2string(xf)
    cat("vcc: ",paste(unlist(xs),collapse=", "),"\n")
    
    xcc <-getecc(x)
    if (length(xcc)){
      xf  <- names2formula(xcc)
      xs  <-formula2string(xf)
      cat("ecc: ",paste(unlist(xs),collapse=", "),"\n")
    }
  }

}

.buildInternalRepresentation <- function(vccN, eccN, dataNames, trace=2){

  vccI <- names2indices(vccN, dataNames)
  eccI <- names2indices(eccN, dataNames)
  
  eccI <- lapply(eccI, function(x) do.call("rbind",x))
  vccI <- lapply(vccI, function(x) do.call("rbind",x))

  if (trace>=6){
    cat("...Internal representation:\n")
    print(tocc(vccI))
    print(tocc(eccI))
  }
  ans <- list(vccI=vccI, eccI=eccI)
  return(invisible(ans))
}

.buildStandardRepresentation <- function(gm, vcc, ecc, trace=2){

  ## Get from formulas to lists
  ##
  gmN   <- formula2names(gm)
  vccN  <- formula2names(vcc)
  eccN  <- formula2names(ecc)
  
  ## Get vertices/edges from gm-spec
  ##
  
  idx         <- sapply(gmN, length)
  gmNvertices <- lapply(unique(unlist(gmN)),list) ## lapply(gmN[idx==1], list)
  x           <- unlist(lapply(gmN[idx>1], names2pairs),recursive=FALSE)
  gmNedges    <- lapply(x, list)

  ## Make standard representation
  ##

  eccN <- c(eccN, gmNedges)

  uuu <- unlist(eccN)  
  uuu <- lapply(uuu, as.list)
  vccN <- c(uuu, vccN, gmNvertices)

  ## Remove entries contained in other entries
  ##
  vccN <-maximalSetL2(vccN)
  eccN <-maximalSetL2(eccN)
  
  varNames <- unique(unlist(c(vccN,eccN)))

  if (trace>=3){
    cat("...Standard representation:\n")
    print(tocc(vccN))
    print(tocc(eccN))
  }
  ans <- list(vccN=vccN, eccN=eccN, varNames=varNames)
  return(ans)


}



.buildDataRepresentation <- function(data=NULL, S=NULL, n=NULL, nodes, type="rcon",
                                     trace=2){
  if (is.null(data) & is.null(S)){
    stop("No data given...\n")
  }
  
  #print(nodes)

  if (!is.null(data)){    
    dataNames <- names(data)  
    S <- cov(data)
    n <- nrow(data)##+1
  } else {
    dataNames <- colnames(S)
  }

  nodes     <- dataNames[sort(match(nodes, dataNames))]
  
  S <- S[nodes, nodes]
  if (type=="rcor"){
    ##cat("Rescaling S\n")
    S <- cov2cor(S)
  }
  
  ans <- list(S=S,n=n, dataNames=rownames(S), nodes=nodes)

  return(ans)
}


  
  ##nodes     <- dataNames[sort(match(nodes, dataNames))]




  
  #cat("vcc: ", cc2str(getSlot(x,'vcc')),"\n")
  #cat("ecc: ", cc2str(getSlot(x,'ecc')),"\n")

  # cat("vcc:\n")
#   aa<-lapply(getSlot(x,"vcc"),cc2str)
#   lapply(aa, cat, "\n")

#   cat("ecc:\n")
#   aa<-lapply(getSlot(x,"ecc"),cc2str)
#   lapply(aa, cat, "\n")

