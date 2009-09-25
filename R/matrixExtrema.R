## Dec 2006
## sort matrix values
##
## based on getLowestPvals in tallyUtilities.R

## N:desired number of entries
## cutoff: cutoff value
## decreasing: the usual sort variable 

matrixExtrema <- function(inMat,N=NULL,cutoff=NULL,decreasing=FALSE){

  matdim <- nrow(inMat)
  labels <- rownames(inMat)

  sortobject <- sort(inMat,index.return=TRUE,decreasing)
  jj <- sortobject$ix
  sortedvals <- sortobject$x
  names(sortedvals) <- NULL

  p1 <- NULL
  p2 <- NULL
  ept <- NULL

  if ( !is.null(N) & is.null(cutoff) ){
    stop.at <- N
  }
  if ( is.null(N) & !is.null(cutoff) ){
    if ( !decreasing){
      ws <- which(sortedvals <= cutoff)
    } else {
      ws <- which(sortedvals >= cutoff)
    }
    if ( length(ws) == 0 ) { return(NULL) }
    stop.at <- max(ws)
  }

  for ( ind in 1:stop.at ){
    ind1 <- ceiling(jj[ind]/matdim)
    ind2 <- jj[ind]-(ind1-1)*matdim
    partner1 <- labels[ind2]
    partner2 <- labels[ind1]
    p1 <- c(p1,partner1)
    p2 <- c(p2,partner2)
    ept <- c(ept, inMat[partner1,partner2])
  }
  pairs <- cbind(p1,p2)
  dimnames(pairs) <- NULL
  return (list(pairs=pairs,vals=ept))
}

