

## assumed globals: matrixLength

countPairs <- function ( geneIDs, pwms, maxsep, minsep, pdist ){

  nm <- length(pwms)
  m.tally <- integer(length=nm)
  names(m.tally) <- pwms

  fams <- unique(sort(as.character(familyMap[pwms])))
  nmf <- length(fams)
  mf.tally <- integer(length=nmf)
  names(mf.tally) <- fams
  
  ## Pairs, retaining order

  m.pair.tally <- matrix(0,nrow = nm,ncol = nm)
  colnames(m.pair.tally) <- pwms
  rownames(m.pair.tally) <- pwms
  
  mf.pair.tally <- matrix(0,nrow = nmf,ncol = nmf)
  colnames(mf.pair.tally) <- fams
  rownames(mf.pair.tally) <- fams
  
  ##
  ## Pairs, disregarding order
  ##

  mdo.pair.tally <- matrix(0,nrow = nm,ncol = nm)
  colnames(mdo.pair.tally) <- pwms
  rownames(mdo.pair.tally) <- pwms
  
  mdof.pair.tally <- matrix(0,nrow = nmf,ncol = nmf)
  colnames(mdof.pair.tally) <- fams
  rownames(mdof.pair.tally) <- fams

  for (entrezID in geneIDs ){

    ##cat("Looking at ",entrezID,"\n")

    ensmusgs <- ensemblIDsOfEntrezIDs[[entrezID]]
    
    m.tally.logical <- logical(length=nm)
    names(m.tally.logical) <- pwms
    mf.tally.logical <- logical(length=nmf)
    names(mf.tally.logical) <- fams
    
    m.pair.logical <- matrix(FALSE,nrow = nm,ncol = nm)
    rownames(m.pair.logical) <- pwms
    colnames(m.pair.logical) <- pwms
    
    mf.pair.logical <- matrix(FALSE,nrow = nmf,ncol = nmf)
    rownames(mf.pair.logical) <- fams
    colnames(mf.pair.logical) <- fams
  
    for ( ensmusg in ensmusgs ){
      
      gene.info <- TRE.OUT[[ensmusg]]$seq.info
      tre <- TRE.OUT[[ensmusg]]$tre.hits; #select the particular tf binding site info for a specific gen
      
      keep.indices <-  which(tre$score > matrixThreshold[tre$tre]) ## which hits are above threshold
      tre <- tre[keep.indices,] ## keep only those

      keep.indices.2 <- which(tre[,'tre'] %in% pwms) ## which hits are in the set of interest (pwms)
      tre <- tre[keep.indices.2,]
      
      NT <- nrow(tre); #the number of hits for this gene
      if(length(NT)!=0 ){ #check to make sure there actually are hits before processing further
        if(NT > 0  ){ #check to make sure there actually are hits before processing further
        
          tre[,'tre.original'] <- tre[,'tre']; #save the original names for later
          tre[,'tre'] <- apply(as.data.frame(tre[,'tre']),1,split1); #just use the prefix for the matrix
        
          found.mats <- unique(sort(tre[,'tre.original']))
          found.families <- unique(sort(as.character(familyMap[found.mats])))
        
          m.tally.logical <- m.tally.logical | ( pwms %in% found.mats )## logical will flip from FALSE to TRUE the first time, then remain on
          mf.tally.logical <- mf.tally.logical | ( fams %in% found.families )
        
          ##we will start closest to the TSS and work backwards
        
          for(tt in 1:(NT-1)){ #loop through all hits and find close neighbors
            tre.o <- tre[tt,]; #the information for the particular TRE
            start.length <- matrixLength[tre.o[,"tre.original"]];   
            tre.p <- tre[(tt+1):NT,]; #looking backwards at all other TREs
            tre.p$tdist <- -tre.p[,'dist'] + tre.o$dist; #compute the distances between hits
            tre.p$adist <- abs(round(0.5*(tre.p[,'dist'] + tre.o$dist))); #compute the average distance from TSS for hits
            tre.p.select <- which((tre.p$tdist<=maxsep & tre.p$tdist>=(minsep+start.length)) & tre.p$adist <= pdist);
            tre.p <- tre.p[tre.p.select,]; #filter out that are too far away from each other, too close to each other, or too far from TSS
            np <- nrow(tre.p); #the number of remaining hits within mdist from the first one
            if(np>0){

              moi <- tre.o[,"tre.original"]
              partners <-  unique(sort(tre.p$tre.original))
              m.pair.logical[moi,] <-  m.pair.logical[moi,] | (pwms %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
         
              moif <- as.character(familyMap[moi])
              partners <- unique(sort(as.character(familyMap[tre.p$tre.original])))    
              mf.pair.logical[moif,] <-  mf.pair.logical[moif,] | (fams %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
            };  #close test on np
          } ## close loop over NT

        } ## close test for NT > 1
      } ## close test for nonzero NT

    } ## close loop over ensembl IDs for the specific entrez ID

    m.tally <- m.tally + m.tally.logical # increment single tally, utilizing integer + logical = integer
    mf.tally <- mf.tally + mf.tally.logical
  
    m.pair.tally <- m.pair.tally + m.pair.logical ## increment pair tally, utilizing integer + logical = integer
    mdo.pair.tally <- mdo.pair.tally + (m.pair.logical | t(m.pair.logical) ) ## disregard order 
  
    mf.pair.tally <- mf.pair.tally + mf.pair.logical ## increment pair tally, utilizing integer + logical = integer
    mdof.pair.tally <- mdof.pair.tally + ( mf.pair.logical | t(mf.pair.logical) ) ## disregard order
  
  };  #close loop over entrez IDs
  
  return(list(m.tally=m.tally,mpair.tally=m.pair.tally,mf.tally=mf.tally,mdo.pair.tally=mdo.pair.tally,mf.pair.tally=mf.pair.tally,mdof.pair.tally=mdof.pair.tally))

} ## close function


## assumed globals: matrixLength

tallyPairs <- function ( indexSet, pwms, maxsep, minsep, pdist ){

  #### Initialize ############
  
  ##
  ## Retaining order of pair
  ##
  nem <- length(pwms)
  em.tally <- integer(length=nem)
  names(em.tally) <- pwms
  em.pair.tally <- matrix(0,nrow = nem,ncol = nem)
  colnames(em.pair.tally) <- pwms
  rownames(em.pair.tally) <- pwms
  
  pwmsPrefix <- unique(sort(as.character(sapply(pwms,split1))))

  nemp <- length(pwmsPrefix)
  emp.tally <- integer(length=nemp)
  names(emp.tally) <- pwmsPrefix
  emp.pair.tally <- matrix(0,nrow = nemp,ncol = nemp)
  colnames(emp.pair.tally) <- pwmsPrefix
  rownames(emp.pair.tally) <- pwmsPrefix

  nemf <- length(familyNames)
  emf.tally <- integer(length=nemf)
  names(emf.tally) <- familyNames
  emf.pair.tally <- matrix(0,nrow = nemf,ncol = nemf)
  colnames(emf.pair.tally) <- familyNames
  rownames(emf.pair.tally) <- familyNames
  
  ##
  ## Disregarding order
  ##
  emdo.pair.tally <- matrix(0,nrow = nem,ncol = nem)
  colnames(emdo.pair.tally) <- pwms
  rownames(emdo.pair.tally) <- pwms
  
  emdop.pair.tally <- matrix(0,nrow = nemp,ncol = nemp)
  colnames(emdop.pair.tally) <- pwmsPrefix
  rownames(emdop.pair.tally) <- pwmsPrefix

  emdof.pair.tally <- matrix(0,nrow = nemf,ncol = nemf)
  colnames(emdof.pair.tally) <- familyNames
  rownames(emdof.pair.tally) <- familyNames


  ## Beging looping through indices

  for ( ss in indexSet ){#loop over genes given by indices
  
    gene.info <- t(as.data.frame(TRE.OUT[[ss]]$seq.info)); # not used but may come in handy
    tre <- TRE.OUT[[ss]]$tre.hits; #select the particular tf binding site info for a specific gene

    keep.indices <- which(tre$tre %in% pwms) ## which hits correspond to 'expressed' matrices
    tre <- tre[keep.indices,] ## keep only those

    keep.indices.2 <-  which(tre$score > matrixThresholdA[tre$tre]) ## which hits are above threshold
    tre <- tre[keep.indices.2,] ## keep only those
 
    NT <- nrow(tre); #the number of hits for this gene
    if(length(NT)!=0 ){ #check to make sure there actually are hits before processing further
      if(NT > 0  ){ #check to make sure there actually are hits before processing further

        tre[,'tre.original'] <- tre[,'tre']; #save the original names for later
        tre[,'tre'] <- apply(as.data.frame(tre[,'tre']),1,split1); #just use the prefix for the matrix

        # single matrix ( or prefixed matrix ) tally 
        found.mats <- unique(sort(tre[,'tre.original']))
        found.prefix.mats <- unique(sort(tre[,'tre']))
        found.families <- unique(sort(as.character(familyMapVT1[found.mats])))

        em.tally <- em.tally + ( pwms %in% found.mats )
        emp.tally <- emp.tally + ( pwmsPrefix %in% found.prefix.mats )
        emf.tally <- emf.tally + ( familyNames %in% found.families ) 
        
        # initialize logical matrix describing whether pair was found 
        em.pair.logical <- matrix(FALSE,nrow = nem,ncol = nem)
        rownames(em.pair.logical) <- pwms
        colnames(em.pair.logical) <- pwms
        emp.pair.logical <- matrix(FALSE,nrow = nemp,ncol = nemp)
        rownames(emp.pair.logical) <- pwmsPrefix
        colnames(emp.pair.logical) <- pwmsPrefix
        emf.pair.logical <- matrix(FALSE,nrow = nemf,ncol = nemf)
        rownames(emf.pair.logical) <- familyNames
        colnames(emf.pair.logical) <- familyNames
        
        ##we will start closest to the TSS and work backwards
        
        for(tt in 1:(NT-1)){ #loop through all hits and find close neighbors
          tre.o <- tre[tt,]; #the information for the particular TRE
          start.length <- matrixLength[tre.o[,"tre.original"]];   
          tre.p <- tre[(tt+1):NT,]; #looking backwards at all other TREs
          tre.p$tdist <- -tre.p[,'dist'] + tre.o$dist; #compute the distances between hits
          tre.p$adist <- abs(round(0.5*(tre.p[,'dist'] + tre.o$dist))); #compute the average distance from TSS for hits
          tre.p.select <- which((tre.p$tdist<=maxsep & tre.p$tdist>=(minsep+start.length)) & tre.p$adist <= pdist);
          tre.p <- tre.p[tre.p.select,]; #filter out that are too far away from each other, too close to each other, or too far from TSS
          np <- nrow(tre.p); #the number of remaining hits within mdist from the first one
          if(np>0){

            moi <- tre.o[,"tre.original"]
            partners <-  unique(sort(tre.p$tre.original))
            em.pair.logical[moi,] <-  em.pair.logical[moi,] | (pwms %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
            
            moip <- tre.o[,"tre"]
            partners <-  unique(sort(tre.p$tre))
            emp.pair.logical[moip,] <-  emp.pair.logical[moip,] | (pwmsPrefix %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on

            moif <- as.character(familyMapVT1[moi])
            partners <- unique(sort(as.character(familyMapVT1[tre.p$tre.original])))    
            emf.pair.logical[moif,] <-  emf.pair.logical[moif,] | (familyNames %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
            
          }; #close test on np
        };

        em.pair.tally <- em.pair.tally + em.pair.logical ## increment pair tally, utilizing integer + logical = integer
        emp.pair.tally <- emp.pair.tally + emp.pair.logical ## increment pair tally, utilizing integer + logical = integer
        
        emdo.pair.tally <- emdo.pair.tally + (em.pair.logical | t(em.pair.logical) ) ## disregard order 
        emdop.pair.tally <- emdop.pair.tally + ( emp.pair.logical | t(emp.pair.logical) ) ## disregard order

        emf.pair.tally <- emf.pair.tally + emf.pair.logical ## increment pair tally, utilizing integer + logical = integer
        emdof.pair.tally <- emdof.pair.tally + ( emf.pair.logical | t(emf.pair.logical) ) ## disregard order
        
      }}; #close original if() on NT

     };  #close loop over genes


  return(list(em.tally=em.tally,em.pair.tally=em.pair.tally,emp.tally=emp.tally,emf.tally=emf.tally,emp.pair.tally=emp.pair.tally,emdo.pair.tally=emdo.pair.tally,emdop.pair.tally=emdop.pair.tally,emf.pair.tally=emf.pair.tally,emdof.pair.tally=emdof.pair.tally))
  
} ## close function


fisherMatrix <- function( mMatrix, MMatrix, n, N){

  matdim <- nrow(mMatrix)
  labels <- rownames(mMatrix)
  ## To be included: Checks for squareness, congruence of two mats
  
  pvals <- matrix(nrow=matdim,ncol=matdim)
  colnames(pvals) <- labels
  rownames(pvals) <- labels
  
  greaterThanExpected <- matrix(FALSE,nrow=matdim,ncol=matdim)
  colnames(greaterThanExpected) <- labels
  rownames(greaterThanExpected) <- labels

  ms <- matrix(0,nrow=matdim,ncol=matdim)
  colnames(ms) <- labels
  rownames(ms) <- labels

  expected <- matrix(0.0,nrow=matdim,ncol=matdim)
  colnames(expected) <- labels
  rownames(expected) <- labels
  
  for ( i in 1:matdim){
    for ( j in 1:matdim ) {
      
      m <- mMatrix[i,j]
      M <- MMatrix[i,j]
      test.mat <- rbind(c(m,M-m),c(n-m,N-n+m-M))
      pvalue <- fisher.test(test.mat)$p.value
      p.choice <- M/N
      p.group <- n/N
      count.expected <- p.choice*p.group*N
      if ( m > count.expected ){ greaterThanExpected[i,j] <- TRUE }
      pvals[i,j] <- pvalue
      ms[i,j] <- m
      expected[i,j] <- count.expected
    }
  }

  return(list(pvals=pvals,greaterThanExpected=greaterThanExpected,ms=ms,expected=expected))

}


### vector version 
fisherVector <- function( mVector, MVector, n, N){

  vecdim <- length(mVector)
  labels <- names(mVector)
  
  pvals <- numeric(length=vecdim)
  names(pvals) <- labels
  
  greaterThanExpected <- logical(length=vecdim) ## defaults to FALSE 
  names(greaterThanExpected) <- labels

  ms <- integer(length=vecdim)
  names(ms) <- labels

  expected <- numeric(length=vecdim)
  names(expected) <- labels
  
  for ( i in 1:vecdim){
    m <- mVector[i]
    M <- MVector[i]
    test.mat <- rbind(c(m,M-m),c(n-m,N-n+m-M))
    pvalue <- fisher.test(test.mat)$p.value
    p.choice <- M/N
    p.group <- n/N
    count.expected <- p.choice*p.group*N
    if ( m > count.expected ){ greaterThanExpected[i] <- TRUE }
    pvals[i] <- pvalue
    ms[i] <- m
    expected[i] <- count.expected
  }

  return(list(pvals=pvals,greaterThanExpected=greaterThanExpected,ms=ms,expected=expected))
}




##
## Get lowest pvalues
##

getLowestPvals <- function(pvalMat,N=NULL,cutoff=NULL){

  matdim <- nrow(pvalMat)
  labels <- rownames(pvalMat)
  
  sortobject <- sort(pvalMat,index.return=TRUE,decreasing=FALSE)
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
    stop.at <- max(which(sortedvals <= cutoff))
  }
  for ( ind in 1:stop.at ){
    ind1 <- ceiling(jj[ind]/matdim)
    ind2 <- jj[ind]-(ind1-1)*matdim
    partner1 <- labels[ind2]
    partner2 <- labels[ind1]
    p1 <- c(p1,partner1)
    p2 <- c(p2,partner2)
    ept <- c(ept, pvalMat[partner1,partner2])
  }
  pairs <- cbind(p1,p2)
  dimnames(pairs) <- NULL
  return (list(pairs=pairs,pvals=ept))
}



## Utility returning matrix element given matrix of pair row columns

matrixValsFromPairs <- function( pairMatrix, matrix ){

  if ( is.null(pairMatrix) | is.null(matrix) ) {
    stop("matrixValsFromPairs: Input Error")
  }

  vals <- vector()
  if ( length(pairMatrix) == 2 ){## single row of pair matrix
    vals <- matrix[pairMatrix[1],pairMatrix[2]]
  } else {
    for ( j in 1: nrow(pairMatrix) ){
      vals <- c(vals,matrix[pairMatrix[j,1],pairMatrix[j,2]])
    }
  }
  return(vals)
}
