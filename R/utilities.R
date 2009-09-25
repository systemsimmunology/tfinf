
##
## General utility functions
##

##
## Utilities for converting names
## map canonical name to common name
"gid2gname" <- function ( gids ){
  gnames <- character(length=length(gids))
  for ( i in 1:length(gids) ) {
    gnames[i] <- gene.names.all[ which(gene.ids.all==gids[i])]
  }
  gnames
}

## map common name to canonical name
"gname2gid" <- function (gnames){

  gids <- character(length=length(gnames))
  for ( i in 1:length(gnames) ) {
    gids[i] <- gene.ids.all[ which(gene.names.all == gnames[i])]
  }
  gids

}


#
# Remove leading and trailing spaces in a string
#
"pruneFlankingSpaces" <- function (inputString){

  ##Count leading spaces
  n.pre <- 0
  foundFirstNonSpace <- FALSE
  i<-1 
  while (i<=nchar(inputString) && foundFirstNonSpace==FALSE ){
    if ( substring(inputString,i,i)==" "){
      n.pre <- n.pre+1
    } else {
      foundFirstNonSpace <- TRUE
    }
    i <-i+1
  }

  ##Count trailing spaces
  n.post <- 0
  foundLastNonSpace <- FALSE  
  i<-nchar(inputString) 
  while (i>=0 && foundLastNonSpace==FALSE ){
    if ( substring(inputString,i,i)==" "){
      n.post <- n.post+1
    } else {
      foundLastNonSpace <- TRUE
    }
    i <-i-1
  }
  outputString <- substring(inputString,n.pre+1,nchar(inputString)-n.post)
  outputString

}
  
## Testing function, to see if pair has been found already
## Takes account of interaction type (but not directionality)
"isPairThere" <- function ( g1, g2, it, e1, e2, itList ){
  pairThere <- FALSE
  if ( length(e1) > 1 ){ ## if length is zero, no need to test
    for ( i in 1:length(e1) ){
      if ( g1==e1[i] && g2==e2[i] && it==itList[i]){
        pairThere <- TRUE
      }
      if ( g2==e1[i] && g1==e2[i] && it==itList[i] ){
        pairThere <- TRUE
      }
    }
  }
  pairThere
}

## Testing function, to see if pair members found, respectively
## in two sets.
"areBothInSets" <- function ( g1, g2, e1, e2){
  pairThere <- FALSE

  if ( length(e1) > 1 && length(e2)>1 ){
    ## if length is zero, no need to test
    if ( (g1 %in% e1) && (g2 %in% e2)){
      pairThere <- TRUE
    }
  }
  pairThere
}



##
## Return list restricted to a subset of arguments
##

"sliceList" <- function ( inList, args ){
  
  returnList <- list()

  for ( arg in args ){
    returnList[[arg]] <- inList[[arg]]
  }

  returnList
}


## standardize rows of matrix
## i.e. transform to mean zero, standard deviation one
"standardizeRows" <- function ( inMatrix ){

  outMatrix <- matrix(nrow=nrow(inMatrix),ncol=ncol(inMatrix))
  rownames(outMatrix) <- rownames(inMatrix)
  colnames(outMatrix) <- colnames(inMatrix)
  
  for ( i in 1:nrow(inMatrix) ){
    tt <- inMatrix[i,]
    outMatrix[i,] <- (tt-mean(tt))/sd(tt)
  }
    
  outMatrix

}
                           

## Quick version of profileplot.R
## For ratio data, the legend
## is assumed to be commonNames
## no title on plot
"ratioProfile" <- function( gids ){

  profileplot(ratios[gids,],gid2gname(gids),"")

}


## elements of a vector having a certain substring
## produces warnings, don't know why
"elementsWithSubString" <- function( desiredSubstring, invec ){
  outvec <- character()
  ssl <- nchar(desiredSubstring)
  for ( st in invec ){
    if ( nchar(st) >= ssl ){
      query=substring(st,1,ssl)
      if( desiredSubstring == query ){
        outvec <- append(outvec,st)
      }
    }
  }
  outvec
}

## 
"stringVec2CommaSeparated" <- function(invec){
  outstring <- invec[1]
  for ( st in invec[2:length(invec)] ){
    outstring <- paste(outstring,",",st,sep="")
  }
  outstring
}
  

                                       
## subset of ratio data above some significance level

"sigSlice" <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){

  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]
  
}

## ids of ratio data following a binary "significance pattern"
## e.g c(1,0,0,0) 
## Should add check for consistency in dimensions. Include above too.


"sigPattern" <- function( lambdaCutoff , boole.pattern, ratioMatrix, lambdaMatrix){

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

###
### sort rows of a matrix by values in a given column 
###
"sortRows" <- function ( inMatrix, colIndex, decreasing= FALSE ){
  o <- order( inMatrix[,colIndex], decreasing=decreasing )
  outMatrix <- inMatrix[o,]
  outMatrix
}

###
### colDiff
### 
### Column differences in a matrix 
### indexing follows m
###

"colDiff" <- function (inMat){

  ## Set up matrix
  n.cols.in  <- ncol(inMat)
  n.cols.out <- n.cols.in*(n.cols.in-1)/2.
  n.rows <- nrow(inMat)
  mat.out <- matrix(data=0,nrow=n.rows,ncol=n.cols.out)
  rownames(mat.out) <- rownames(inMat)
  col.headers.in <- colnames(inMat)
  col.headers.out <- character(length=n.cols.out)
  
  ## Populate with differences
  colcounter <- 1
  for ( i in 1:(n.cols.in-1) ){

    colname.i <- col.headers.in[i]
    
    for ( j in (i+1):n.cols.in ){
      
      mat.out[,colcounter] <- inMat[,i]-inMat[,j]
      
      colname.j <- col.headers.in[j]
      colname.new <- paste(colname.i,"-",colname.j,sep="")
      col.headers.out[colcounter] <- colname.new 
      
      colcounter <- colcounter + 1
      
    }
  }

  colnames(mat.out) <- col.headers.out 
  mat.out
  
}


### oneColSort
###
### Sort matrix rows by values in a single column
### col: The specified column. Can be numeric index or string index.

"oneColSort" <- function ( col, inMat, decreasing=FALSE){
  sortedOrder <- sort(inMat[,col],index=TRUE,decreasing=decreasing)$ix
  outMat <- inMat[sortedOrder,]
  outMat
}


### oneRowSort
###
### Sorts matrix columns by values in a single row
### row: The specified row. Can be numeric index or string index.

"oneRowSort" <- function ( row, inMat, decreasing=FALSE){
  sortedOrder <- sort(inMat[row,],index=TRUE,decreasing=decreasing)$ix
  outMat <- inMat[,sortedOrder]
  outMat
}
