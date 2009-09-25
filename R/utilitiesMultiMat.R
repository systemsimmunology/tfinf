##
## Expand out vector values into a multidimensional matrix
##
## April 2007

## d: a vector of matrix dimensions
## ds: an index of d, denoting which direction the vector lies
## inVec: vector values
## we must have length(inVec)==d[ds]

expandOut <- function( inVec, ds, d){

  if ( length(inVec) != d[ds] ){
    stop("Inconsistency between vector and multidim array")
  }

  ## Consruct empty multidimensional array
  nd <- length(d)
  nel <- 1
  for ( n in 1:nd ){
    nel <- nel* d[n]
  }
  em <- numeric(length=nel)
  dim(em) <- d
  
  for ( i in 1:d[ds] ){

    ## Need to construct string to address the multi-array

    as <- "em["
    ## Blanks for leading dimensions
    if ( ds != 1 ){
      for ( n in 1:(ds-1) ){
        as <- paste(as,",",sep="")
      }
    }
    ## ds
    as <- paste(as,as.character(i),sep="")
    ## Balnks for trailing dimensions
    if ( ds != nd ){
      for ( n in (ds+1):nd ){
        as <- paste(as,",",sep="")
      }
    }
    as <- paste(as,"]",sep="")
    
    ## assignment 
    eval( parse( text=paste(as,"<-",as.character(inVec[i]))))
  
  }    

  em
}

##
## Expand out matrix values into a multidimensional matrix
##
## April 2007

## d: a vector of matrix dimensions
## ds: an index of d, denoting which direction the vector lies
## inVec: vector values

## we must have dim(inVec) == c(d[ds1],d[ds2])

expandOut2 <- function( inMat, ds1, ds2, d){

  d1 <- dim(inMat)[1]
  d2 <- dim(inMat)[2]
  
  if ( (d1!=d[ds1]) | (d2!=d[ds2]) ){
    stop("Inconsistency between matrix and multidim array")
  }

  ## Consruct empty multidimensional array
  nd <- length(d)
  nel <- 1
  for ( n in 1:nd ){
    nel <- nel* d[n]
  }  
  em <- numeric(length=nel)
  dim(em) <- d

  ## for indexing, need to keep track of index rank order
  d.lower <- min(ds1,ds2)
  d.higher <- max(ds1,ds2)
  d.sep <- d.higher-d.lower

  for ( i in 1:d[d.lower] ){
    for ( j in 1:d[d.higher] ){
    
      ## Need to construct string to address the multi-array
      
      as <- "em["
      ## Blanks for leading dimensions
      if ( d.lower != 1 ){
        for ( n in 1:(d.lower-1) ){
          as <- paste(as,",",sep="")
        }
      }
      
      ## d.lower
      as <- paste(as,as.character(i),sep="")
      
      ## Between
      bets <- paste(rep(",",d.sep),collapse="")
      as <- paste(as,bets,sep="")
      
      ## d.higher
      as <- paste(as,as.character(j),sep="")
    
      ## Balnks for trailing dimensions
      if ( d.higher != nd ){
        for ( n in (d.higher+1):nd ){
          as <- paste(as,",",sep="")
        }
      }
      as <- paste(as,"]",sep="")
      
      ## assignment 
      eval( parse( text=paste(as,"<-",as.character(inMat[i,j]))))
  
    }    
  }
  
  em
}
  

  

  

  
  


