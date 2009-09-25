

##
## Selective loading of members of RData object
##

selectiveLoad <- function(objectNameArray,RDataFileLocation){
  jj <- load(RDataFileLocation,.GlobalEnv)
  chuckthese <- setdiff(jj,objectNameArray)
  rm( list=chuckthese, envir=.GlobalEnv)
}

##
## Memory use by all objects
## 
memoryUse <- function(unit="mb"){

  allobj.string <- ls(env=.GlobalEnv)
  n.objs <- length(allobj.string)

  objSizeVec <- integer(length=n.objs)
  names(objSizeVec) <- allobj.string

  for ( as in allobj.string ){
    objSizeVec[as] <- object.size(eval(parse(text=as)))
  }

  totalUse <- as.integer(sum(objSizeVec))
  
  if ( unit == "b" ){
    returnVal <- totalUse
  } else if ( unit == "kb" ){
    returnVal <- round(totalUse/1000)
  } else if ( unit == "mb" ){
    returnVal <- round(totalUse/1e6,1)
  }
  returnVal
}



##
## Get largest objects, according to object.size()
## 
largestObjects <- function(n=15,unit="kb"){

  allobj.string <- ls(env=.GlobalEnv)
  n.objs <- length(allobj.string)

  if ( n > n.objs ){ n=n.objs } ##for any large n, return all  
  
  returnVec <- integer(length=n.objs)
  names(returnVec) <- allobj.string

  for ( as in allobj.string ){
    returnVec[as] <- object.size(eval(parse(text=as)))
  }

  returnVec <- sort(returnVec,decreasing=TRUE)
  returnVec <- returnVec[1:n]

  if ( unit == "b" ){
    returnVec <- returnVec
  }
  
  if ( unit == "kb" ){
    returnVec <- returnVec/1000.
    returnVec <- sapply(returnVec,round)
  }

  if ( unit == "mb" ){
    returnVec <- returnVec/1e6
    returnVec <- sapply(returnVec,round,2)
  }
  
  returnVec
}
