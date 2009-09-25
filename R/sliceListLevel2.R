
## Feb 2007
## Slice list at level two for a given  label
## There is a fancy way to do this in R that I cannot recall


sliceListLevel2 <- function( inlist, label ){

  outList <- list()
  level1Labels <- names(inlist)
  
  for ( i in 1:length(inlist) ) { 
    
    level1label <- level1Labels[i]
    inlistElement <- inlist[[i]]
    
    if( !is.null( inlistElement[[label]] ) ){
      outList[[level1label]] <- inlistElement[[label]]
      
    }
  }
  
  outList
  
}
