## Required
## psois.all, nbrs, genemappings, featureMatrix
## n.pairs, N, Mpair, fampairs, Msteps
## ESmaxcubeMgrid

## If not prepared in an .RData file
source("./initializeBC.R")

load(paste(seq.dir,"featureMatrix.RData",sep="/"))
load(paste(seq.dir,"ESmaxcubeMgrid.RData",sep="/"))

maxn <-  600 ## The maximum N that we will consider
mst <- 20 # the stepsize on the M grid used in the randomization

Mpair <- apply(featureMatrix,1,sum)
allMs <- Mpair
fampairs <- names(Mpair)
n.pairs <- length(fampairs)

Msteps <- seq(mst,mst*round(N/mst),mst) ## steps used in randomization
numMsteps <- length(Msteps)

psois.all <- names(nbrs)

ptm <- proc.time()
paired.features <- list()
for ( psoi in psois.all   ){
##for ( psoi in repProbes.cname["Il12b"] ){
  nbs <- nbrs[[psoi]]
  nn <- length(nbs)

  if ( nn > maxn ){
    nbs <- nbs[1:maxn]
    nn <- maxn
  }
  
  if ( nn>5 ){

    nbs.eid <- as.character(ncbiID[nbs])
    fm.sub <- featureMatrix[,nbs.eid]
    mmat <- t(apply(fm.sub,1,cumsum))
    matn <- t(matrix(rep(1:nn,n.pairs),ncol=n.pairs))
    TPmat <- mmat
    FPmat <- matn - mmat
    FNmat <- Mpair - mmat
    TNmat <- N - TPmat - FPmat - FNmat
    FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
    TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
    ESmat <- TPRmat - FPRmat 
    ES <- ESmat
    Esmax <- apply(ES,1,max)
    indsatmax <- apply(ES,1,which.max)

    ## Need to index mmat, and matn by these
    ## Only non-loop way I could find is to index the matrices
    ## by single variables
    ## Due to the way R enumerates (along columns)
    ## need to do this for the transposed matrix
    transinds <- indsatmax + (0:(n.pairs-1))*nn
    msatmax <- t(mmat)[transinds]
    names(msatmax) <- fampairs
    nsatmax <- t(matn)[transinds]
    names(nsatmax) <- fampairs

    ## Determine which ESmax scores are significant
    pair.keep <- character()
    score.keep <- numeric()
    pval.keep <- numeric()
    ms.keep <- integer()
    ns.keep <- integer()
 
    for ( i in 1:n.pairs ){
      pval <- ESmaxPval(Esmax[i],allMs[i],nn,ESmaxcubeMgrid )
      if ( pval <= 0.05 ){
        pair.keep <- c(pair.keep,fampairs[i])
        score.keep <- c(score.keep,Esmax[i])
        names(pval) <- fampairs[i]
        pval.keep <- c(pval.keep,pval)
        ms.keep <- c(ms.keep,msatmax[i])
        ns.keep <- c(ns.keep,nsatmax[i])
      }
    }
    feature.save <- rbind(score.keep,pval.keep,ms.keep,ns.keep)
    paired.features[[psoi]] <- feature.save 
    
  } ## end if (nn>5)
  
} ## end loop over psois

save(paired.features,file=paste(seq.dir,"PairedFeatures.RData",sep="/"))

cat(proc.time() - ptm,"\n")

