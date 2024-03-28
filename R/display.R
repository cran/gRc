#' @title Plot rcox object
#' @param x An rcox object
#' @param y Ignored
#' @param ... Currently not used
#'
#' @examples
#' gm  = ~al:an:st
#' vcc = list(~me+st, ~ve+an, ~al)
#' ecc = list(~me:ve+me:al, ~ve:al+al:st)
#' m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching')
#' plot(m1)
#' 
#' @method plot rcox
#' @export
plot.rcox <- function(x, y, ...){
    ecc_lst <- getSlot(x, 'ecc')
    vcc_lst <- getSlot(x, 'vcc')
    
    gen <- c(unlist(vcc_lst, recursive=FALSE), unlist(ecc_lst, recursive=FALSE))
    gen <- remove_redundant(gen)
    gg <- ug(gen)
    ## plot(gg)
    
    
    
    ecc_cols<-topo.colors(length(ecc_lst))
    reps <- sapply(ecc_lst, length)      
    ecc_cols[reps==1] <- "black"
    ecc_cols
    cc <- ecc_cols[rep(1:length(ecc_lst), times=reps)]
    
    el <- as_edgelist(gg)
    idx <- apply(el, 1, FUN=function(r){
        get.edge.ids(gg, r)
    })
    E(gg)[idx]$color <- cc
    
    
    coef_ <- coef(x)
    if (is.null(coef_))
        coef_ <- 1:c(length(ecc_lst) + length(vcc_lst))
    
    o <- order(coef_[1:length(vcc_lst)])
    vcc_cols <- heat.colors(length(vcc_lst))
    vcc_cols <- vcc_cols[o]
    
    nms <- nodes(gg)
    reps <- sapply(vcc_lst, length)  
    vcc_cols[reps==1] <- "white"
    cols <- vcc_cols[rep(1:length(vcc_lst), times=reps)]
    idx <- lapply(vcc_lst, match, nms)
    V(gg)$color[unlist(idx)]     <- cols
    
    plot(gg)
    return(invisible(gg))
    
}
