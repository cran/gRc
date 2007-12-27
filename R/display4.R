plot.rcox <- function(x,y,...){

  m2 <- x
#  eccList <- getCC2(m2,'ecc','list')
#  vccList <- getCC2(m2,'vcc','list')

  eccList <- getSlot(m2,'ecc')
  vccList <- getSlot(m2,'vcc')

  coef <- coef(x)
  
  
  #l <- sapply(vccList, length)
  
  o <- order(coef[1:length(vccList)])
  vccColors <- heat.colors(length(vccList))
  vccColors <- vccColors[o]
  
  V <- unlist(vccList)
  V <- V[order(V)]
  vertexColors <- NULL
  for (i in 1:length(vccList)){
    tmp <- vccList[[i]]
    #print(tmp)
    if (length(tmp)==1){
       vcolor <- "white"    
    } else {
       vcolor <- vccColors[i] 
    }
    d <- c(vstr = rep(vcolor, length(tmp)))
    names(d) <- tmp
    vertexColors <- c(vertexColors, d)    
  }

  nAttrs <- list()
  nAttrs$fillcolor <- vertexColors

  G <- new("graphNEL", nodes=V)
  
  eccColors<-topo.colors(length(eccList))

  edgeColors <- NULL
  if (length(eccList)>0){
    o <- order(coef[-(1:length(vccList))])
    eccColors <- eccColors[o]

    for (i in 1:length(eccList)){
      tmp <- eccList[[i]]; ltmp <- length(tmp)
      for (j in 1:ltmp){
        ee <- tmp[[j]]
        ee <- ee[order(ee)]
                                        #print(ee)
        G <- addEdge(ee[1], ee[2], G, weight=1)
        estr <- paste(ee[1],"~",ee[2],sep='')
        if (ltmp > 1){
          ecolor <- eccColors[i]
          d <- c(estr = ecolor)
          names(d) <- estr
          edgeColors <- c(edgeColors, d)                  
      }
      }
    }
  }

  #print(edgeColors)
  if (!is.null(edgeColors))
    eAttrs <- list(color=edgeColors)
  else
    eAttrs <- list()

  #cc <- rep(3, length(edgeColors))
  #names(cc) <- names(edgeColors)
  #eAttrs$lwd <- cc
  
  plot(G, "neato", nodeAttrs = nAttrs, edgeAttrs = eAttrs)

  #return(G)

}
