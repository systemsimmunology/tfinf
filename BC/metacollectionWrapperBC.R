

##
## Use precomputed single.features and paired.features
## to populate objects in sigPairedSites.RData
##
## Jan 2009

##source("~/macrophage/TFinfluence/R/utilities.R")

## single. and paired. features indexed by psoi
## need to index by ensemblID

load("SingleFeatures.RData") ## single.features
load("PairedFeatures.RData") ## paired.features

load("featureMatrix.RData")
Mpair <- apply(featureMatrix,1,sum)
load("../sequence_data/featureMatrix.mf.RData")
Msingles <- apply(featureMatrix.mf,1,sum)

psois.have.single <- names(single.features)
psois.have.paired <- names(paired.features)
psois.all <- union(psois.have.single,psois.have.single)

ensids.paired <- unlist(ensemblIDsOfEntrezIDs[ncbiID[psois.have.paired]])

split2 <- function(splitme,splitchar='_'){ #parse binding site names by dropping "_" and everything behind it
 strsplit(splitme,split=splitchar,fixed=TRUE)[[1]];
}

##
## Pairs 
##
metampairs.bc <- list()
metapvals.bc <- list()
metams.bc <- list()
metans.bc <- list()
metaexpected.bc <- list()
for ( psoi in psois.have.paired ){

  ensid <- unlist(ensemblIDsOfEntrezIDs[ncbiID[psoi]])   ## at this stage, there should just be one of these
  matricks <-  paired.features[[psoi]]

  ## Pairs
  pairs.hyphenated <- colnames(matricks)
  pair.matrix <- t(sapply(pairs.hyphenated,split2,"-"))
  rownames(pair.matrix) <- NULL 
  metampairs.bc[[ensid]] <- pair.matrix  

  ## Pvals
  pvals <- matricks["pval.keep",]
  names(pvals) <- NULL
  metapvals.bc[[ensid]] <- pvals

  ## m (at ES maximum)
  ms <- matricks["ms.keep",]
  names(ms) <- NULL
  metams.bc[[ensid]] <- ms

  ## n (at ES maximum)
  ns <- matricks["ns.keep",]
  names(ns) <- NULL
  metans.bc[[ensid]] <- ns

  ## expected count (at ES maximum)
  bigMs <- Mpair[pairs.hyphenated]
  expecteds <- ns * bigMs/N
  names(expecteds) <- NULL
  metaexpected.bc[[ensid]] <- expecteds
  
}


##
## Singles
##

metamsingles.bc <- list()
metapvals.singles.bc <- list()
metams.singles.bc <- list()
metans.singles.bc <- list()
metaexpected.singles.bc <- list()
for ( psoi in psois.have.single ){

  ensid <- unlist(ensemblIDsOfEntrezIDs[ncbiID[psoi]])   ## at this stage, there should just be one of these
  matricks <-  single.features[[psoi]]

  ## Singles
  metamsingles.bc[[ensid]] <-  colnames(matricks)

  ## Pvals
  pvals <- matricks["pval.keep",]
  names(pvals) <- NULL
  metapvals.singles.bc[[ensid]] <- pvals

  ## m (at ES maximum)
  ms <- matricks["ms.keep",]
  names(ms) <- NULL
  metams.singles.bc[[ensid]] <- ms

  ## n (at ES maximum)
  ns <- matricks["ns.keep",]
  names(ns) <- NULL
  metans.singles.bc[[ensid]] <- ns

  ## expected count (at ES maximum)
  bigMs <- Msingles[colnames(matricks)]
  expecteds <- ns * bigMs/N
  names(expecteds) <- NULL
  metaexpected.singles.bc[[ensid]] <- expecteds
  
}

save(metampairs.bc,metapvals.bc,metams.bc,metaexpected.bc,metamsingles.bc,metapvals.singles.bc,metams.singles.bc,metaexpected.singles.bc,file="sigPairedSitesBC.RData")

## nbrs is the ***Set** of neighbours in psois. Named list. Names are also psoi.
#save(nbrs,maximumdist,metampairs,metapvals,metams,metaexpected,metamsingles,metapvals.singles,metams.singles,metaexpected.singles,file="sigPairedSites.RData")

