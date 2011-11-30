
source("./metacollectionInitialize.R")

metampairs <- list()
metapvals <- list()
metams <- list()
metaexpected <- list()
metaM <- list()


metamsingles <- list()
metapvals.singles <- list()
metams.singles <- list()
metaexpected.singles <- list()
metaM.singles <- list()

nbrs <- list()
maximumdist <- numeric()

num.nbrs <- numeric()
counter <- 0 

##
##
## Important : Need to set N  
##


N <- length(ensemblIDsOfEntrezIDs) 

##
##
##  Parameters 
## 
##

pval.threshold <- 0.05
maxdist <- 100
maxsep <- 100
minsep <- 1
pdist <- 1000

min.number.correlates <- 3
max.number.correlates <- 30

min.correlation.dist <- 0.15

##ese = "16193" ## IL-6

##
## experiment with cutting down memory usagge??
## TRE.OUT <- TRE.OUT[cherries.ensids] Doesn't seem to help memory usage

## The Entrez IDs to predict
eids <- expressed.scanned.egid

## Probesets from which to find neighboring expression patterns

##psois.dn <- cherries
psois.dn <- expressed.scanned.ps
corMat <- cordistA[psois.dn,psois.dn]

## load expressedMats and expressedFamilyMats

for ( eid in eids  ) {

  counter <- counter + 1 
  write(c(eid,counter,counter/length(eids)),file='status.txt'); #report status

  ##cat("\nEntrezID:",eid,"\n")
  
  ensemblIDs <- ensemblIDsOfEntrezIDs[[eid]]

  ## get expression neighbors
  psoi <- as.character(repProbes.ncbiID[eid])
  
  psois <- getDistMatNeighbors(psoi,corMat,N=max.number.correlates,cutoff=min.correlation.dist)
  nbrs[[psoi]] <- psois
  if ( length(psois) >= 1 ){
    num.nbrs[psoi] <- length(psois)
    maximumdist[psoi] <- max(cordistA[psoi,psois])
  }
  
  if ( length(psois) >= min.number.correlates ){ ## for now, skip if insufficient number of expression correlates
    
    geneIDs <- as.character(ncbiID[psois])
    ## countPairs has TRE.OUT as a global
    bb <- countPairs(geneIDs,expressedMats,maxsep,minsep,pdist)

    ### Now, for each promoter, identify the paired sites, and determine if they are significant

    
    for ( ensemblID in ensemblIDs ){

      ##cat("Ensembl ID:",ensemblID,"\n")
      hm <- getFilteredHitMat(ensemblID,TRE.OUT,matrixThreshold=matrixThreshold)
      ret <- getSigPairs(hm,pval.threshold=pval.threshold,maxdist=maxdist,minsep=minsep,pwm.string="family",keep.order=FALSE,M=length(psois),N=N)
      ret2 <- getSigSingles(hm,pval.threshold=pval.threshold,maxdist=maxdist,minsep=minsep,pwm.string="family",M=length(psois),N=N)

      if ( !is.null(ret) ){
        metampairs[[ensemblID]] <- ret$mpairs
        metapvals[[ensemblID]] <- ret$pvals
        metams[[ensemblID]] <- ret$ms
        metaexpected[[ensemblID]] <- ret$expected
        metaM[[ensemblID]] <- length(psois)

      }
      if ( !is.null(ret2) ){
        metamsingles[[ensemblID]] <- ret2$msingles
        metapvals.singles[[ensemblID]] <- ret2$pvals
        metams.singles[[ensemblID]] <- ret2$ms
        metaexpected.singles[[ensemblID]] <- ret2$expected
        metaM.singles[[ensemblID]] <- length(psois)
      }
    }
  }
}

ofile <- file.path(Sys.getenv("TFINF"),"sequence_data/sigPairedSites.Nov2011.RData")
save(nbrs,maximumdist,metampairs,metapvals,metams,metaexpected,metamsingles,metapvals.singles,metams.singles,metaexpected.singles,file=ofile)

