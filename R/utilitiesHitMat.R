
## From a tag in as subset of  extract hm meeting threshold and containing expressedMats
getFilteredHitMat <- function ( ensid, treoutSubset, matrixThreshold ) {
  hitMat <- treoutSubset[[ensid]]$tre.hits
  keepers <- which(hitMat[,"tre"] %in% expressedMats)
  hitMat.filt1 <- hitMat[keepers,]
  ## Filter away the matrix hist below threshold
  nf1 <- dim(hitMat.filt1)[1]
  keepers2 <- vector()
  for ( i in 1:nf1 ){
    moi <- hitMat.filt1$tre[i]
    if ( hitMat.filt1$score[i] >= matrixThreshold[moi] ){
      keepers2 <- c(keepers2,i)
    }  
  }
  hitMat.filt2 <- hitMat.filt1[keepers2,]
  nf2 <- dim(hitMat.filt2)[1]
  hm <- hitMat.filt2
  return(hm)
}

## filter hm by matrixThreshold
getScoreFilteredHitMat <- function (hitMat, matrixThreshold ) {
  nf1 <- dim(hitMat)[1]
  keepers2 <- numeric()
  for ( i in 1:nf1 ){
    moi <- hitMat$tre[i]
    if ( hitMat$score[i] >= matrixThreshold[moi] ){
      keepers2 <- c(keepers2,i)
    }  
  }
  hitMat.filt2 <- hitMat[keepers2,]
  nf2 <- dim(hitMat.filt2)[1]
  hm <- hitMat.filt2
  return(hm)
}

## return 'TRE.OUT' filtered with a threshold matrix 
## If matrices are not found, that promoter is not returned
getScoreFilteredTE <- function ( TEin, matrixThreshold ){
  returnList <- list()
  for ( label in names(TEin) ){
    te <- TEin[[label]]
    seq.info <- te$seq.info
    hmStart <- te$tre
    hmFiltered <- getScoreFilteredHitMat(hmStart,matrixThreshold)
    if ( nrow(hmFiltered) > 0 ){
      teReturn <- list()
      teReturn$seq.info <- seq.info
      teReturn$tre <- hmFiltered
      returnList[[label]] <- teReturn
    }
  }
  return(returnList)
}

## return hitMat with a subset of matrices only
getMatrixFilteredHitMat <- function (hitMat, matrixSet ) {
  keepers <- which(hitMat$tre %in% matrixSet)
  if ( length(keepers) > 0 ){
    return(hitMat[keepers,])
  }

}

## return 'TRE.OUT' with a subset of matrices only
## If matrices are not found, that promoter is not returned
getMatrixFilteredTE <- function ( TEin, matrixSet ){
  returnList <- list()
  for ( label in names(TEin) ){
    te <- TEin[[label]]
    seq.info <- te$seq.info
    hmStart <- te$tre
    hmFiltered <- getMatrixFilteredHitMat(hmStart,matrixSet)
    if ( ! is.null(hmFiltered) ){
      if ( nrow(hmFiltered) > 0 ){
        teReturn <- list()
        teReturn$seq.info <- seq.info
        teReturn$tre <- hmFiltered
        returnList[[label]] <- teReturn
      }
    }
  }
  return(returnList)
}
    
    
    
  
  
getMatrixFilteredHitMat <- function (hitMat, matrixSet ) {
  keepers <- which(hitMat$tre %in% matrixSet)
  return(hitMat[keepers,])
}


## Does hit matrix have a pair within a maxsep
## also used minsep
hasClosePair <- function(mat1,mat2,hm,pwm.string="full",keep.order=TRUE,maxsep=100,minsep=0){

  hm.tre.full <- hm$tre ## we'll use this only to get exact matrixLengths
    
  if (pwm.string=="prefix"){
    hm.tre.here <- as.character(sapply(hm$tre,split1))
  } else if (pwm.string=="full") {
    hm.tre.here <- hm$tre
  } else if ( pwm.string=="family") {
    hm.tre.here <-  as.character(familyMap[hm$tre])
  }
  mat1.tss.dists <- hm$dist[which(hm.tre.here==mat1)]
  mat1.tss.lengths <- matrixLength[hm.tre.full[which(hm.tre.here==mat1)]]
  mat2.tss.dists <- hm$dist[which(hm.tre.here==mat2)]
  mat2.tss.lengths <- matrixLength[hm.tre.full[which(hm.tre.here==mat2)]]
  if ( length(mat1.tss.dists) >= 1  && length(mat2.tss.dists) >= 1 ){
    if ( pwm.string=="full"){
      mL1 <- matrixLength[mat1]
      mL2 <- matrixLength[mat2]
    } else if ( pwm.string =="prefix"  ){
      mL1 <- matrixLengthPrefix[mat1]
      mL2 <- matrixLengthPrefix[mat2]
    } else if ( pwm.string == "family" ){
      mL1 <- matrixLengthPrefix[mat1]  ## Not quite correct, but will have to do for now
      mL2 <- matrixLengthPrefix[mat2]  ## Not quite correct, but will have to do for now
    }
    min12sep <- minPairDistanceV(mat2.tss.dists,mat1.tss.dists,mat2.tss.lengths+minsep,mat1.tss.lengths+minsep,keep.order) ## Re-check the 1-2 here.
    ##min12sep <- minPairDistance(mat2.tss.dists,mat1.tss.dists,mL2+minsep,mL1+minsep,keep.order) ## Re-check the 1-2 here.
    if ( !is.null(min12sep) ){
      if ( min12sep <= maxsep ) {
      ##if ( (min12sep <= maxsep) & (min12sep >= (mL+minsep)) ) {
        return (TRUE)
      } else {
        return (FALSE)
      }
    } else { return(FALSE) }
   } else {
     return(FALSE)
   }
}

## 
## Essentially the same as hasClosePair, but returns the pair distance
closestPairDistance <- function(mat1,mat2,hm,pwm.string="full",keep.order=TRUE,minsep=0){
  
  hm.tre.full <- hm$tre ## we'll use this only to get exact matrixLengths
    
  if (pwm.string=="prefix"){
    hm.tre.here <- as.character(sapply(hm$tre,split1))
  } else if (pwm.string=="full") {
    hm.tre.here <- hm$tre
  } else if ( pwm.string=="family") {
    hm.tre.here <-  as.character(familyMap[hm$tre])
  }
  mat1.tss.dists <- hm$dist[which(hm.tre.here==mat1)]
  mat1.tss.lengths <- matrixLength[hm.tre.full[which(hm.tre.here==mat1)]]
  mat2.tss.dists <- hm$dist[which(hm.tre.here==mat2)]
  mat2.tss.lengths <- matrixLength[hm.tre.full[which(hm.tre.here==mat2)]]
  if ( length(mat1.tss.dists) >= 1  && length(mat2.tss.dists) >= 1 ){
    if ( pwm.string=="full"){
      mL1 <- matrixLength[mat1]
      mL2 <- matrixLength[mat2]
    } else if ( pwm.string =="prefix"  ){
      mL1 <- matrixLengthPrefix[mat1]
      mL2 <- matrixLengthPrefix[mat2]
    } else if ( pwm.string == "family" ){
      mL1 <- matrixLengthPrefix[mat1]  ## Not quite correct, but will have to do for now
      mL2 <- matrixLengthPrefix[mat2]  ## Not quite correct, but will have to do for now
    }
    min12sep <- minPairDistanceV(mat2.tss.dists,mat1.tss.dists,mat2.tss.lengths+minsep,mat1.tss.lengths+minsep,keep.order) ## Re-check the 1-2 here.
    ##min12sep <- minPairDistance(mat2.tss.dists,mat1.tss.dists,mL2+minsep,mL1+minsep,keep.order) ## Re-check the 1-2 here.
    if ( !is.null(min12sep) ){
      return(min12sep)
    } else {
      return(NaN)
    }
  } else {
    return(NaN)
  }
  
}

 
## tell me which target scans have a certain matrix 
whichScansHaveMatHit <- function(mat, treoutSubset ){
  returnVec <- character()
  for ( ensid in names(treoutSubset) ){
    mathits <- treoutSubset[[ensid]]$tre.hits$tre
    if ( mat %in% mathits ){
      returnVec <- c(returnVec,ensid)
    }
  }
  return(returnVec)
}


# filter scans with no hits  
removeNULLscans <- function(treoutSubset ){
  baddies <- which(unlist(lapply(lapply(treoutSubset,"[[","tre.hits"),is.null)))
  goodies <- setdiff(1:length(treoutSubset),baddies)
  treoutSubset[goodies]
}


##
## Sep17, Discovered that getSigPairs still uses em instead of m
## Have corrected this here below
## 

## assumed gobals: MatrixLength, indexSubset ,N, bb$m.pair.tally, m.pair.tally, and other pair tallies

## pwm.string can take two values
## full: use full form of matrix
## prefix: use prefix form of matrix
## keep.order can take two values
## TRUE: keep track of order
## FALSE: disregard order
getSigPairs <- function(hm,pval.threshold=0.01,maxdist=100,minsep=1,pwm.string="full",keep.order=TRUE,M,N){
 
  mpairs <- NULL
  nhits <- nrow(hm)
   
  for ( n in 1:(nhits-1) ){
    hm.soi <- hm[n,]; #the information for the particular TRE "site of interest"
    hm.up <- hm[(n+1):nhits,]; #looking backwards at all other TREs
    tdist <-  -hm.up[,'dist'] + hm.soi$dist #compute the distances between hits
    moi <-  hm.soi$tre
    soi.length <- matrixLength[moi]
     
    if ( pwm.string == "prefix" ){
      moi <- split1(moi)
    } else if ( pwm.string == "family" ){
      moi <- as.character(familyMap[moi])
    }
     
    crit1 <- (tdist <= maxdist) #filter out that are too far away from each other
    crit2 <- (tdist > (soi.length+minsep )) #filter out those that are too near to each other
 
    if ( TRUE %in% (crit1&crit2) ){ ## proceed only if there is a row that satisfies the criteria
      hm.up.filt <- hm.up[which(crit1&crit2),];
      near.sites <- hm.up.filt$tre
      if ( pwm.string == "prefix" ){
        near.sites <- as.character(sapply(near.sites,split1))
      } else if ( pwm.string == "family" ){
        near.sites <- as.character(familyMap[near.sites])
      }
 
      for ( pwm in near.sites ){ ## create pair collection
        mpairs <- rbind(mpairs,c(moi,pwm))
      }
       
    }
  }
   
 
 
  ## Unique site pairs
  if ( length(mpairs) ){
    mpairs <- unique(mpairs) ## unique rows
  }
 
   
  ## if order is disregarded
  if ( !keep.order ) {
    if ( length(mpairs) > 1 ){
      mpairs <- unique(t(apply(mpairs, 1, sort))) ## sort first within a row, then finds unique rows
    }
  }
   
  if ( length(mpairs) ){
    ## find appropriate pair tally
    if ( pwm.string=="full" & keep.order ){
      bb.pair.tally <- bb$m.pair.tally
      pair.tally <- m.pair.tally
    } else if ( pwm.string=="family" & keep.order ){
      bb.pair.tally <- bb$mf.pair.tally
      pair.tally <- mf.pair.tally
    } else if ( pwm.string=="full" & !keep.order ){
      bb.pair.tally <- bb$mdo.pair.tally
      pair.tally <- mdo.pair.tally
    } else if ( pwm.string=="family" & !keep.order ){
      bb.pair.tally <- bb$mdof.pair.tally
      pair.tally <- mdof.pair.tally
    }
    ## pvalue for each pair
    pvalstats <- fisherVector( matrixValsFromPairs(mpairs,bb.pair.tally),matrixValsFromPairs(mpairs,pair.tally),M,N)
    ## filter the pairs based on pvalue
    keep.indices <- which( (pvalstats$pvals<pval.threshold) & pvalstats$greaterThanExpected )
    if ( length(keep.indices) >= 1 ){
      nm.pairs <- length(keep.indices)
      mpairs <- mpairs[keep.indices,]
      pvals <- pvalstats$pvals[keep.indices]
      ms <- pvalstats$ms[keep.indices]
      expected <- pvalstats$expected[keep.indices]
      return(list(mpairs=mpairs,pvals=pvals,ms=ms,expected=expected))
    }
  }
  return(NULL)
   
}

getSigSingles <- function(hm,pval.threshold=0.01,maxdist=100,minsep=1,pwm.string="family",M,N){
 
  msingles <- NULL
  nhits <- nrow(hm)
 
  tre <- hm$tre
  found.mats <- unique(sort(tre))
  found.families <- unique(sort(as.character(familyMap[found.mats])))
 
  if ( pwm.string == "family" ){
    msingles = found.families
  } else if ( pwm.string == "full" ){
    msingles = found.mats
  } else {
    stop("Error: pwm.string not supported")
  }
   
  if ( length(msingles) ){
    ## find appropriate pair tally
    if ( pwm.string=="full"  ){
      bb.tally <- bb$m.tally
      tally <- m.tally
    } else {
      bb.tally <- bb$mf.tally
      tally <- mf.tally
    }
    ## pvalue for each pair
    pvalstats <- fisherVector(bb.tally[msingles],tally[msingles],M,N)
    ## filter the singles based on pvalue
    keep.indices <- which( (pvalstats$pvals<pval.threshold) & pvalstats$greaterThanExpected )
    if ( length(keep.indices) >= 1 ){
      nm.singles <- length(keep.indices)
      msingles <- msingles[keep.indices]
      pvals <- pvalstats$pvals[keep.indices]
      ms <- pvalstats$ms[keep.indices]
      expected <- pvalstats$expected[keep.indices]
      return(list(msingles=msingles,pvals=pvals,ms=ms,expected=expected))
    }
  }
  return(NULL)
   
}


## Mappings to TFs

getTFsForTE <- function( teSubset, tfsubset=transfac.tfs.expressed ){
  returnList <- list()
  for ( ensid  in names(teSubset) ) {
    psoi  <- as.character(repProbes.ncbiID[entrezIDofEnsemblID[ensid]])
    if ( length(psoi) > 1 ) {
      cat ("Trouble ahead. Multiple psois for this ensid,",ensid,"\n" )
    }
    hm1 <- teSubset[[ensid]]$tre
    tfset <- getTFsForHitMat(hm1, tfsubset )
    if ( length(tfset) >= 1 ){
      returnList[[psoi]] <- tfset
    }
  }
  returnList
}

getTFsForHitMat <- function (hitMat, tfsubset=transfac.tfs.expressed ) {
  allowed.tfs <- tfsubset ## (psois)
  mats1 <- hitMat$tre
  tfs <- unique(sort(as.character(unlist(pwm2tf.tte.ps[mats1]))))
  tfs.out <- as.character(cname.compare[intersect(tfs,allowed.tfs)])
  return(tfs.out)
}


## Miniminum distance between elements of a vector
##
## For ordered mode, think of
## "the distance to the closest a upstream of b"
##
## distance returned as a positive number
## fla, flb, Feature length associated with a and b, respectively
minPairDistance <- function(a,b,fla,flb,keep.order=TRUE){
  na <- length(a)
  nb <- length(b)
  min.dist <- NULL
  if ( na==1 | nb== 1){
    mm <- a - b
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(abs(mm))
    } else if ( keep.order & (TRUE %in% (mm <0) )) { # need to test that there actually is negative number
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  } else {
    a.repped <- t(matrix(rep(a,nb),nrow=na,ncol=nb))
    b.repped <- matrix(rep(b,na),nrow=nb,ncol=na)
    mm <- a.repped - b.repped
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(min(abs(mm)))
    } else if ( keep.order & (TRUE %in% (mm <0)) ){
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  }
  return(min.dist)
}



## Miniminum distance between elements of a vector
## Same as above, but fla and flb are vectors
##
## For ordered mode, think of
## "the distance to the closest a upstream of b"
##
## distance returned as a positive number
## fla, flb, Feature length associated with a and b, respectively
minPairDistanceV <- function(a,b,fla,flb,keep.order=TRUE){
  na <- length(a)
  nb <- length(b)
  min.dist <- NULL
  if ( na==1 | nb== 1){
    mm <- a - b
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(abs(mm))
    } else if ( keep.order & (TRUE %in% (mm <0) )) { # need to test that there actually is negative number
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  } else {
    a.repped <- t(matrix(rep(a,nb),nrow=na,ncol=nb))
    b.repped <- matrix(rep(b,na),nrow=nb,ncol=na)
    fla.repped <- t(matrix(rep(fla,nb),nrow=na,ncol=nb))
    flb.repped <- matrix(rep(flb,na),nrow=nb,ncol=na)

    mm <- a.repped - b.repped
    mm[which(mm>=0 & mm <= fla.repped )] <- Inf
    mm[which(mm<=0 & mm >= -flb.repped )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(min(abs(mm)))
    } else if ( keep.order & (TRUE %in% (mm <0)) ){
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  }
  return(min.dist)
}

