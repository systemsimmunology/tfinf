

###
### Nov 2006: Use list construct to make 'objects'

###   sc
###   sc$metamsingles,sc$metapvals.singles,sc$metams.singles, and TF singles, to be named
###   pc
###   pc$metampairs, pc$metapvals, pc$metams and TF pairs, to be named

##source("./utilitiesInteractions.R")

### For a metapair list, find which promoters have the matrix mat
findSingleMatInPairHits <- function(mat,metampairs){
  nmeta <- length(metampairs)
  promovec <- character()
  for ( i in 1:nmeta ){
    if ( mat  %in% metampairs[[i]] ) { promovec <- c(promovec,names(metampairs[i])) }
  }
  return(promovec)
}

### For a metamsingles list, find which promoters have the matrix mat
findSingleMatHits <- function(mat,metamsingles){
  nmeta <- length(metamsingles)
  promovec <- character()
  for ( i in 1:nmeta ){
    if ( mat  %in% metamsingles[[i]] ) { promovec <- c(promovec,names(metamsingles[i])) }
  }
  return(promovec)
}

### For a metapair list, find which promoters have the matrix mat, at below a certain pvalue threshold
findSingleMatHitsPval <- function(mat,pval.thresh,metampairs,metapvals){
  nmeta <- length(metampairs)
  promovec <- character()
  for ( i in 1:nmeta ){
    cat( i,"\n")
    mpairs <- metampairs[[i]]
    pvalvec <- metapvals[[i]]
    fmp <- filterMpairs( pval.thresh, mpairs, pvalvec )
    if ( mat %in% fmp ){ promovec <- c(promovec,names(metampairs[i])) }
  }
  return(promovec)
}

## filter pairmatrix collection by pvalue
filterMpairCollection <- function ( pval.thresh, mpaircollection, pvalcollection ){
  returnList <- list()
  for ( pss in names(mpaircollection) ) {
    mpairmat <- mpaircollection[[pss]]
    pvalvec <- pvalcollection[[pss]]
    fm <- filterMpairs( pval.thresh, mpairmat, pvalvec )
    if ( !is.null(fm)) {
      returnList[[pss]] <- fm
    }
  }
  return(returnList)
}

## Filter pairmatrix by pvalue
filterMpairs  <- function( pval.thresh, mpairmat, pvalvec ){
  if ( length(which(pvalvec<=pval.thresh))>1 ){
    if ( !is.vector(mpairmat) ) {
      indices.keep <- which(pvalvec <= pval.thresh )
      return( mpairmat[indices.keep,])
    } else {
      if ( pvalvec <= pval.thresh ){ return(mpairmat)}
    }
  }
  return(NULL) ## if none of the above are met
}


## Filter single hit matrix by pvalue
filterMsingles  <- function( pval.thresh, msinglemat, pvalvec ){
  ## turns out we only need pvalvec, as it is labeled
  return(names(which(pvalvec <= pval.thresh )))
}
  

## Find partners to a single pwm in a pairmatrix
findMatPartners <- function( mat, mpairmat ){
  indices.matLeft <- which(mpairmat[,1]==mat)
  indices.matRight <- which(mpairmat[,2]==mat) 
  partners <- c(mpairmat[indices.matLeft,2],mpairmat[indices.matRight,1])
  return (partners)
}

## Give p-vales for the pairs involving matrix mat 
findMatPartnerPvals <- function( mat, mpairmat, pvalvec ){
  indices.matLeft <- which(mpairmat[,1]==mat)
  indices.matRight <- which(mpairmat[,2]==mat) 
  partners <- c(mpairmat[indices.matLeft,2],mpairmat[indices.matRight,1])
  pvals <- c(pvalvec[indices.matLeft],pvalvec[indices.matRight])
  names(pvals) <- partners
  return (pvals)
}

## filter PC collection by pvalue
filterPCbyPval <- function ( pc, pval.thresh ){

  mpaircollection <- pc$metampairs
  pvalcollection <- pc$metapvals
  emcollection <- pc$metams

  returnList <- list()
  returnList$metampairs <- list()
  returnList$metapvals <- list()
  returnList$metams <- list()

  for ( promoter in names(mpaircollection) ) {

    mpairmat <- mpaircollection[[promoter]]
    pvalvec <- pvalcollection[[promoter]]
    emvec <- emcollection[[promoter]]

    if ( length(which(pvalvec<=pval.thresh))>= 1 ){
      if ( !is.vector(mpairmat) ) {
        indices.keep <- which(pvalvec <= pval.thresh )
        if ( length(indices.keep) > 0 ){
          returnList$metampairs[[promoter]] <- mpairmat[indices.keep,]
          returnList$metapvals[[promoter]] <- pvalvec[indices.keep]
          returnList$metams[[promoter]] <- emvec[indices.keep]
        }
      } else { ## if single value
        if ( pvalvec <= pval.thresh ){
          returnList$metampairs[[promoter]] <- mpairmat
          returnList$metapvals[[promoter]] <- pvalvec
          returnList$metams[[promoter]] <- emvec
        }
      }
    }

  }

  if ( length(returnList$metapvals)==0 ){
    return ( NULL )
  } else {
    return(returnList)
  }

}

## filter SC collection by pvalue
filterSCbyPval <- function ( sc, pval.thresh ){

  msinglecollection <- sc$metamsingles
  pvalcollection <- sc$metapvals.singles
  emcollection <- sc$metams.singles

  returnList <- list()
  returnList$metamsingles <- list()
  returnList$metapvals.singles <- list()
  returnList$metams.singles <- list()

  for ( promoter in names(msinglecollection) ) {

    msingles <- msinglecollection[[promoter]]
    pvalvec <- pvalcollection[[promoter]]
    emvec <- emcollection[[promoter]]

    indices.keep <- which(pvalvec <= pval.thresh )
    if ( length(indices.keep) > 0 ){
      returnList$metamsingles[[promoter]] <- msingles[indices.keep]
      returnList$metapvals.singles[[promoter]] <- pvalvec[indices.keep]
      returnList$metams.singles[[promoter]] <- emvec[indices.keep]
    }
  }

  if ( length(returnList$metapvals.singles)==0 ){
    return ( NULL )
  } else {
    return(returnList)
  }

}


## filter PC collection by required matrices

## reqm is vector c(mat1,mat2) of two matrices

### Looks like this has yet to be completed

filterPCbyMatrices <- function ( pc , reqm ){

  mpaircollection <- pc$metampairs
  pvalcollection <- pc$metapvals
  emcollection <- pc$metams

  returnList <- list()
  returnList$metampairs <- list()
  returnList$metapvals <- list()
  returnList$metams <- list()

  for ( promoter in names(mpaircollection) ) {

    ##cat(promoter,"\n")   
    mpairmat <- mpaircollection[[promoter]]
    pvalvec <- pvalcollection[[promoter]]
    emvec <- emcollection[[promoter]]

   
    if ( !is.vector(mpairmat) ){
      for ( i in 1:nrow(mpairmat) ) {
        if ( identical(sort(reqm),sort(mpairmat[i,]))) {
          returnList$metampairs[[promoter]] <- mpairmat[i,]
          returnList$metapvals[[promoter]] <- pvalvec[i]
          returnList$metams[[promoter]] <- emvec[i]
        }
      }
    } else {
      if ( identical(sort(reqm),sort(mpairmat))) {
        returnList$metampairs[[promoter]] <- mpairmat
        returnList$metapvals[[promoter]] <- pvalvec
          returnList$metams[[promoter]] <- emvec
      }
    }
  }

  if ( length(returnList$metapvals)==0 ){
    return ( NULL )
  } else {
    return(returnList)
  }

}


## all pairs in two sets, maintaining order
expandPairs <- function(set1,set2){
  return.obj <- NULL
  for ( s1 in set1 ){
    for (s2 in set2 ){
      return.obj <- rbind(return.obj,c(s1,s2))
    }
  }
  return(return.obj)
}


## input PC, return TFs
## PC is indexed by ensids
## TFs indexed by psois

## possible criteria:
## 1. set of TFs
## specify tfs.subset if only  a subset of tfs is to be considered
## tf.subset.both=TRUE if both are required
## tf.subset.both=FALSE if only one is required ( Not yet implemented )
## 2. interactome

createTFsetFromPairs <- function ( pc, interactome=TRUE, tfsubset=transfac.tfs.expressed, tfsubset.both=TRUE, noConnections=NULL ){
  
  allowed.tfs <- setdiff(as.character(cname.compare[tfsubset]),noConnections)
  
  metampairs <- pc$metampairs
  metapvals <- pc$metapvals ## not used, but retained  in case we need it
  metams <- pc$metams ## not used, but retained  in case we need it

  ensids <- names(metampairs)

  ##
  ## Compute filtered gene sets by the revised cutoff
  ##

  returnList <- list()
  
  for ( ensid  in names(metampairs)  ) {

    mpairs <- metampairs[[ensid]]
    ppvals <- metapvals[[ensid]]
    ms <- metams[[ensid]]

    psoi  <- as.character(repProbes.ncbiID[entrezIDofEnsemblID[ensid]])
    if ( length(psoi) > 1 ) {
      cat ("Trouble ahead. Multiple psois for this ensid,",ensid,"\n" )
    }
    if ( length(psoi) > 1 ) {
      cat ("Trouble ahead. Multiple psois for this ensid:",ensid,"\n" )
    }
        
    ##
    ## Transcription factor pairs
    ##
    tfpaircollection <- character()

    if ( is.vector(mpairs ) ){ ## if one of one possible

      soi.tfs <- fam2tf.tte.gname[[mpairs[1]]]
      soi.tfs <- intersect( soi.tfs, allowed.tfs ) 
      near.tfs <- fam2tf.tte.gname[[mpairs[2]]]
      near.tfs <- intersect( near.tfs, allowed.tfs )
      if ( (length(soi.tfs)>0) & (length(near.tfs)>0) ){
        if ( interactome ) {
          pairs <- grabPairs(soi.tfs,near.tfs) 
        } else {
          pairs <- expandPairs(soi.tfs,near.tfs) 
        }
        if ( length(pairs) != 0 ) {
          same.tf.logical = pairs[,1]==pairs[,2]
          pairs <- pairs[ !same.tf.logical, ]
          ##if ( TRUE %in% (pairs[,1]==pairs[,2]) )stop("Error: self-pair")
          tfpaircollection <- rbind(tfpaircollection,pairs)
        }
      }
    } else { 

      for ( mpair.index in 1:nrow(mpairs) ){
        soi.tfs <- fam2tf.tte.gname[[mpairs[mpair.index,1]]]
        soi.tfs <- intersect( soi.tfs, allowed.tfs ) 
        near.tfs <- fam2tf.tte.gname[[mpairs[mpair.index,2]]]
        near.tfs <- intersect( near.tfs, allowed.tfs ) 
        if ( (length(soi.tfs)>0) & (length(near.tfs)>0) ){
          if ( interactome ) {
            pairs <- grabPairs(soi.tfs,near.tfs)
          } else {
            pairs <- expandPairs(soi.tfs,near.tfs)
          }
          if ( length(pairs) != 0 ) {
            same.tf.logical = pairs[,1]==pairs[,2]
            pairs <- pairs[ !same.tf.logical, ]
            ##if ( TRUE %in% (pairs[,1]==pairs[,2]) )stop("Error: self-pair")
            tfpaircollection <- rbind(tfpaircollection,pairs)
          }
        }
        
      }
    }

    rownames(tfpaircollection) <- NULL

    if ( length(tfpaircollection) > 1 ){
      tfpaircollection <- unique(t(apply(tfpaircollection, 1, sort))) ## sort first within a row, then finds unique rows
    }
  
    if ( length(tfpaircollection) > 0 ){

      ## Some psois correspond to multiple ensids
      ## We are going to treat all hypotheses for a given psoi on an equal footing, and stack them 
      if ( is.null(returnList[[psoi]]) ){
        returnList[[psoi]] <- tfpaircollection
      } else {
        returnList[[psoi]] <- rbind(returnList[[psoi]],tfpaircollection)
        returnList[[psoi]] <- unique(t(apply(returnList[[psoi]],1,sort)))
      }

    }
  }
    
  returnList
}

## input PC, return TFs
## PC is indexed by ensids
## TFs indexed by psois

## possible criteria:
## 1. set of TFs

createTFsetFromSingles <- function ( sc, tfsubset=transfac.tfs.expressed ){

  allowed.tfs <- as.character(cname.compare[tfsubset])
  
  metamsingles <- sc$metamsingles
  metapvals.singles <- sc$metapvals.singles ## not used, but retained  in case we need it
  metams.singles <- sc$metams.singles ## not used, but retained  in case we need it

  ensids <- names(metamsingles)

  ##
  ## Compute filtered gene sets by the revised cutoff
  ##

  returnList <- list()
  
  for ( ensid  in ensids  ) {

    msingles <- metamsingles[[ensid]]
    ppvals <- metapvals.singles[[ensid]]
    ms <- metams.singles[[ensid]]

    psoi  <- as.character(repProbes.ncbiID[entrezIDofEnsemblID[ensid]])
    if ( length(psoi) > 1 ) {
      cat ("Trouble ahead. Multiple psois for this ensid:",ensid,"\n" )
    }
    if ( is.na(psoi) ) {
      cat ("Trouble ahead. This ensid has no repProbe:",ensid,"\n" )
    }
        
    ##
    ## Transcription factor singles
    ##

    tfs <- unique(sort(as.character(unlist(fam2tf.tte.gname[msingles]))))
    tfs <- intersect( tfs, allowed.tfs )
  
    if ( length(tfs) > 0 ){

      ## Some psois correspond to multiple ensids
      ## We are going to treat all hypotheses for a given psoi on an equal footing, and concatenate them
      if ( is.null(returnList[[psoi]]) ){
        returnList[[psoi]] <- tfs
      } else {
        returnList[[psoi]] <- c(returnList[[psoi]],tfs)
        returnList[[psoi]] <- unique(sort(returnList[[psoi]]))
      }
    }
    
  }
    
  returnList
}
