## Required
## psois.all, nbrs, genemappings, featureMatrix
## n.pairs, N, Mpair, fampairs, Msteps
## ESmaxcubeMgrid

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
r.dir <- file.path(Sys.getenv("TFINF"),"R")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"annotation.objects.RData",sep="/"))
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(annot.dir,"allMouseMappings.August2009.RData",sep="/"))
source("./utilitiesBinaryClassification.R")

load(paste(seq.dir,"featureMatrix.mf.RData",sep="/"))
load(paste(seq.dir,"ESmaxcubeMgrid.RData",sep="/"))

Msingles <- apply(featureMatrix.mf,1,sum)
allMs <- Msingles
famsingles <- names(Msingles)
n.singles <- length(famsingles)

Msteps <- seq(20,20*round(N/20),20)
numMsteps <- length(Msteps)

psois.all <- names(nbrs)

ptm <- proc.time()
single.features <- list()

for ( psoi in psois.all   ){
##for ( psoi in repProbes.cname["Il12b"] ){
  nbs <- nbrs[[psoi]]
  nn <- length(nbs)

  if ( nn > 600 ){
    nbs <- nbs[1:600]
    nn <- 600
  }
  
  if ( nn>5 ){

    nbs.eid <- as.character(ncbiID[nbs])
    fm.sub <- featureMatrix.mf[,nbs.eid]
    mmat <- t(apply(fm.sub,1,cumsum))
    matn <- t(matrix(rep(1:nn,n.singles),ncol=n.singles))
    TPmat <- mmat
    FPmat <- matn - mmat
    FNmat <- Msingles - mmat
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
    transinds <- indsatmax + (0:(n.singles-1))*nn
    msatmax <- t(mmat)[transinds]
    names(msatmax) <- famsingles
    nsatmax <- t(matn)[transinds]
    names(nsatmax) <- famsingles
    
    ## Determine which ESmax scores are significant
    singles.keep <- character()
    score.keep <- numeric()
    pval.keep <- numeric()
    ms.keep <- integer()
    ns.keep <- integer()
    
    for ( i in 1:n.singles ){
      pval <- ESmaxPval(Esmax[i],allMs[i],nn,ESmaxcubeMgrid)
      if ( pval <= 0.05 ){
        singles.keep <- c(singles.keep,famsingles[i])
        score.keep <- c(score.keep,Esmax[i])
        names(pval) <- famsingles[i]
        pval.keep <- c(pval.keep,pval)
        ms.keep <- c(ms.keep,msatmax[i])
        ns.keep <- c(ns.keep,nsatmax[i])
      }
    }
    feature.save <- rbind(score.keep,pval.keep,ms.keep,ns.keep)
    single.features[[psoi]] <- feature.save 
    
  } ## end if (nn>5)
  
} ## end loop over psois

save(single.features,file=paste(seq.dir,"SingleFeatures.RData",sep="/"))

cat(proc.time() - ptm,"\n")

