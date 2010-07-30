##
## Use precomputed single.features and paired.features
## to populate objects in sigPairedSites.RData
##

seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
load(paste(seq.dir,"SingleFeatures.RData",sep="/"))
load(paste(seq.dir,"PairedFeatures.RData",sep="/"))
load(paste(seq.dir,"featureMatrix.RData",sep="/"))
load(paste(seq.dir,"featureMatrix.mf.RData",sep="/"))
Mpair <- apply(featureMatrix,1,sum)
Msingles <- apply(featureMatrix.mf,1,sum)

psois.have.single <- names(single.features)
psois.have.paired <- names(paired.features)
psois.all <- union(psois.have.single,psois.have.single)

ensids.paired <- unlist(ensemblIDsOfEntrezIDs[ncbiID[psois.have.paired]])

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
  pair.matrix <- t(sapply(pairs.hyphenated,tokens,"-"))
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

ofile <- paste(seq.dir,"sigPairedSitesBC.RData",sep="/"))
save(metampairs.bc,metapvals.bc,metams.bc,metaexpected.bc,metamsingles.bc,metapvals.singles.bc,metams.singles.bc,metaexpected.singles.bc,file=ofile)
