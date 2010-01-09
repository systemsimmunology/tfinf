

                                       
## subset of ratio data above some significance level

sigSlice <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){

  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]
  
}

## ids of ratio data following a binary "significance pattern"
## e.g c(1,0,0,0) 
## Should add check for consistency in dimensions. Include above too.


sigPattern <- function( lambdaCutoff , boole.pattern, ratioMatrix, lambdaMatrix){

  n.conds <- ncol(ratioMatrix)
  n.rows <- nrow(ratioMatrix)
  ids <- row.names(ratioMatrix)
  
  index.pattern <- boole.pattern+1
  desired.logic.pattern <- c(FALSE,TRUE)[index.pattern]

  sig.logic.matrix <- lambdaMatrix>lambdaCutoff
  collection <- character()

  for ( n in 1:n.rows ){
    row.pattern=sig.logic.matrix[n,]
    comp.pattern=(row.pattern == desired.logic.pattern)
    if ( sum(comp.pattern) == n.conds ){
      collection <- append(collection,ids[n]) 
    }
  }

  collection
  
}
