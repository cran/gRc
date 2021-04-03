#############################################################################
#' Joining and splitting of colour classes in RCOX models
#' 
#' Test for joining of two colour classes (of a specific type) by testing if
#' their corresponding parameters are not significantly different.  Split a
#' colour class and test how much this changes the fit of the model.
#' 
#' @author  Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name join-split
#############################################################################
#'
#' @aliases join1 split1
#'
#' @param object An RCOX model, an object of class RCOX
#' @param scope A specification of colour classes which should be
#'   considered for joining/splitting. If NULL, then all colour
#'   classes are considered.
#' @param type Either "ecc" for edge colour classes or "vcc" for
#'   vertex colour classes.
#' @param stat Either "wald" for a Wald statistic or "dev" for
#'   deviance statistic.
#' @param details Control the amount of output
#'
#'
#' @return A list with entries: \item{tab}{A data frame with the test
#'   results} \item{cc}{A list of colour classes}
#' @note Note that the keyword 'stat' is not available for split1
#'   because this function expands the current and hence the Wald
#'   statistic is not available.  Note also that join1 is simply a
#'   wrapper for comparecc applied to edge colour classes.

#' @seealso \code{\link{rcox}}, \code{\link{update}},
#'   \code{\link{comparecc}}
#' @keywords htest
#' @examples
#' 
#' data(math)
#' g1     <- ~me:ve:al+al:st:an
#' m1     <- rcox(gm=g1, data=math)
#' join1(m1)
#' 
#' gm  = ~al:an:st
#' vcc = list(~me+st, ~ve+an)
#' ecc = list(~me:ve+me:al, ~ve:al+al:st)
#' m2 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, type="rcon")
#' split1(m2)
#' 
NULL

####
#### What do we do here???
####

#' @export
split1 <- function(object, scope=NULL, type='ecc', details=1){
    logL0 <- logL(object)
    
    if (missing(scope)){
        cc.list   <- getSlot(object, type)
        cc.list   <- cc.list[sapply(cc.list, length) > 1]
    } else {
        cc.list   <- .addccnames(formula2names(scope), type=type)
    }
    if (length(cc.list) == 0)  return(NULL) ## sensible return value???

    ## create a list of model objects, compute lr-statistic for each
    ans <- lapply(cc.list, function(cc){
      switch(type, 
        "ecc"={
          mtmp <- update(object, splitecc=list(cc))
        },
        "vcc"={
          mtmp <- update(object, splitvcc=list(cc))
        }
      )
      list(cc2str(cc), 2 * (logL(mtmp) - logL0), length(cc) - 1)
    })

    #print(ans)
    ans <- .nestedList2df(ans)

    names(ans)    <- c("cc", "X2", "df")
    #print(ans)

    cc.list       <- .addccnames(cc.list,type)
    ans$cc        <- names(cc.list)
    rownames(ans) <- 1:nrow(ans)
    ans           <- .addStat(ans, n=dataRep(object,"n"))
    
    #print(ans)
    attr(ans, "ccterms") <- cc.list
    #ans
    
    ans2 <- list(tab=ans, cc=.addccnames(cc.list,type),
      details=details)
    class(ans2) <- c("statTable", "data.frame")
    ans2
}  

#' @export
join1 <- function(object, scope=NULL, type='ecc',details=1, stat='wald'){

  if (missing(scope)){
    cc.list   <- getSlot(object, type)
  } else {
    cc.list   <- .addccnames(formula2names(scope), type=type)
  }

  ans <- comparecc(object, cc1=cc.list, cc2=cc.list, type=type, stat=stat, details=details)
  return(ans)
}

