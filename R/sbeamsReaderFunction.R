
### May need to remove AFFX..

getMatricesFromSbeamsFile <- function ( sbeams.file, sbeamsCondOrder,
                                         desiredSuffixes,cond.names.short){

  data.sbeams <- read.table(sbeams.file,header=TRUE,sep='\t')
  pset <- as.character(data.sbeams$reporter_name)

  colheaders <- names(data.sbeams) 
  ##Can get colheaders without reading in all rows 
  ##colheaders <- names(read.table(sbeams.file,header=TRUE,sep='\t',nrow=1))

  matrixList <- list()
  for ( desiredSuffix in desiredSuffixes ){
    slice <- getMatSlice(desiredSuffix,colheaders,data.sbeams,pset,cond.names.short,sbeamsCondOrder)
    matrixList[[desiredSuffix]] <- slice
  }
  matrixList
}

## getSlice for given suffix
getMatSlice <- function(desiredSuffix,colheaders,data.sbeams,pset,
                          cond.names.short,sbeamsCondOrder){ 

  nsuff <- nchar(desiredSuffix)

  ## find columns with this suffix, and re-order from sbeams order
  colSuffs <- mapply(tailEnd,colheaders,rep(nsuff,length(colheaders)))
  needed.conds <- names(which(colSuffs==desiredSuffix)[sbeamsCondOrder])

  ## Slice through data matrix and relabel rows
  matSlice <- apply(as.matrix(data.sbeams)[,needed.conds],2,as.numeric)

  rownames(matSlice) <- pset
  colnames(matSlice) <- cond.names.short

  matSlice

}



