

## returns pairs
grabPairs <- function(soi.tfs,near.tfs){

  ns <- length(soi.tfs)
  nt <- length(near.tfs)
  
  if ( ns ==0 || nt == 0){
    cat("Error: Missing TFs as input to grabPairs function\n")
  }

  if ( (FALSE %in% (soi.tfs %in% colnames(tf.dist.cn))) || (FALSE %in% (near.tfs %in% colnames(tf.dist.cn)))){
      cat("Error: No distance available ( in grabPairs function )\n")
      cat(setdiff(soi.tfs,colnames(tf.dist.cn)),"\n")
      cat(setdiff(near.tfs,colnames(tf.dist.cn)),"\n")
  }

  distMat <- tf.dist.cn[soi.tfs,near.tfs] ## Can be a vector in some cases
  collection <- character()
  ## Unfortunately, matrices collapsing to vectors (ns=1 and/or nt=1)
  ## have posed the usual R headaches with losing names
  ## and nrow/ncol evaluating to NaN and so on
  ## ( using data.frame didn't solve this completely )
  ## Hence all the case statements and cumbersome code
  
  if ( ns==1 && nt==1 ){
    partner1 <- soi.tfs
    partner2 <- near.tfs
    if ( distMat==1 | distMat==2 ){
      collection <- rbind(character(),c(partner1,partner2))
    }
  }## end ns==1, nt==1 case

  if ( ns>1 && nt==1){
    partner1 <- near.tfs
    one.indices <- which(distMat==1)
    two.indices <- which(distMat==2)
    keepers <- c(one.indices,two.indices)
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(near.tfs,length(keepers)),soi.tfs[keepers]) 
      rownames(collection) <- NULL
    }
    #if ( length(one.indices) != 0 ){
    #  for ( one.index  in 1:length(one.indices) ){
    #    partner2 <- soi.tfs[one.index]
    #    collection <- rbind(collection,c(partner1,partner2))
    #  }
    #}
    #if ( length(two.indices) != 0 ){
    #  for ( two.index  in 1:length(two.indices) ){
    #    partner2 <- soi.tfs[two.index]
    #    collection <- rbind(collection,c(partner1,partner2))
    #  }
    #}
  }## end ns>1, nt==1 case

  if ( nt>1 && ns==1){
    partner1 <- soi.tfs
    one.indices <- which(distMat==1)
    two.indices <- which(distMat==2)
    keepers <- c(one.indices,two.indices)
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(soi.tfs,length(keepers)),near.tfs[keepers]) 
      rownames(collection) <- NULL
    }
    ##if ( length(one.indices) != 0 ){
    ##  for ( one.index  in 1:length(one.indices) ){
    ##    partner2 <- near.tfs[one.index]
    ##    collection <- rbind(collection,c(partner1,partner2))
    ##  }
    ##}
    ##if ( length(two.indices) != 0 ){
    ##  for ( two.index  in 1:length(two.indices) ){
    ##    partner2 <- near.tfs[two.index]
    ##    collection <- rbind(collection,c(partner1,partner2))
    ##  }
    ##}
  }## end nt>1, ns==1 case

  if ( nt>1 &&  ns> 1 ){
    one.indices <- which(distMat==1,arr.ind=TRUE) ## 
    two.indices <- which(distMat==2,arr.ind=TRUE) 

    if ( length(one.indices) != 0 ){
      for ( one.index.row  in 1:nrow(one.indices) ){
        partner1 <- soi.tfs[one.indices[one.index.row,1]]
        partner2 <- near.tfs[one.indices[one.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
  
    if ( length(two.indices) != 0 ){
      for ( two.index.row  in 1:nrow(two.indices) ){
        partner1 <- soi.tfs[two.indices[two.index.row,1]]
        partner2 <- near.tfs[two.indices[two.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
    
  }## end nt>1, ns>1 case

  colnames(collection) <- NULL
  return (collection) 
}
## returns pairs


## Grab pairs, nearest neighbors only
grabPairsNN <- function(soi.tfs,near.tfs){

  ns <- length(soi.tfs)
  nt <- length(near.tfs)
  
  if ( ns ==0 || nt == 0){
    cat("Error: Missing TFs as input to grabPairs function\n")
  }

  if ( (FALSE %in% (soi.tfs %in% colnames(tf.dist.cn))) || (FALSE %in% (near.tfs %in% colnames(tf.dist.cn)))){
      cat("Error: No distance available ( in grabPairs function )\n")
  }

  distMat <- tf.dist.cn[soi.tfs,near.tfs] ## Can be a vector in some cases
  collection <- character()
  ## Unfortunately, matrices collapsing to vectors (ns=1 and/or nt=1)
  ## have posed the usual R headaches with losing names
  ## and nrow/ncol evaluating to NaN and so on
  ## ( using data.frame didn't solve this completely )
  ## Hence all the case statements and cumbersome code
  
  if ( ns==1 && nt==1 ){
    partner1 <- soi.tfs
    partner2 <- near.tfs
    if ( distMat==1 ){
      collection <- rbind(character(),c(partner1,partner2))
    }
  }## end ns==1, nt==1 case

  if ( ns>1 && nt==1){
    partner1 <- near.tfs
    one.indices <- which(distMat==1)
    keepers <- one.indices
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(near.tfs,length(keepers)),soi.tfs[keepers]) 
      rownames(collection) <- NULL
    }
  }## end ns>1, nt==1 case

  if ( nt>1 && ns==1){
    partner1 <- soi.tfs
    one.indices <- which(distMat==1)
    keepers <- one.indices
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(soi.tfs,length(keepers)),near.tfs[keepers]) 
      rownames(collection) <- NULL
    }
  }## end nt>1, ns==1 case

  if ( nt>1 &&  ns> 1 ){
    one.indices <- which(distMat==1,arr.ind=TRUE)

    if ( length(one.indices) != 0 ){
      for ( one.index.row  in 1:nrow(one.indices) ){
        partner1 <- soi.tfs[one.indices[one.index.row,1]]
        partner2 <- near.tfs[one.indices[one.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
  
  }## end nt>1, ns>1 case

  colnames(collection) <- NULL
  return (collection) 
}

