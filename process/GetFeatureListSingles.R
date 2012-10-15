## Required
## psois.all, nbrs, genemappings, featureMatrix
## n.pairs, N, Mpair, fampairs, Msteps
## ESmaxcubeMgrid

## If not prepared in an .RData file
proc.dir <- file.path(Sys.getenv("TFINF"),"process")
out.dir <- file.path(Sys.getenv("TFINF"),"data/process_output")
load(paste(out.dir,"initializeBC.RData",sep="/"))

load(paste(out.dir,"featureMatrix.mf.RData",sep="/"))
load(paste(out.dir,"ESmaxcubeMgrid.RData",sep="/"))

maxn <-  600 ## The maximum N that we will consider
mst <- 20 # the stepsize on the M grid used in the randomization

Msingles <- apply(featureMatrix.mf,1,sum) ## For each matrix, the number of occurences (over all psois)
allMs <- Msingles
famsingles <- names(Msingles)
n.singles <- length(famsingles)

Msteps <- seq(mst,mst*round(N/mst),mst) ## steps used in randomization
numMsteps <- length(Msteps)

psois.all <- names(nbrs)

ptm <- proc.time()
single.features <- list()

for ( psoi in psois.all   ){
  nbs <- nbrs[[psoi]] ## set of neighbors of the psoi
  nn <- length(nbs) ## number of neigbhors

  if ( nn > maxn ){
    nbs <- nbs[1:maxn]
    nn <- maxn
  }

  ##
  ## Matrices have n.singles rows and nn columns
  ##
  if ( nn>5 ){
    nbs.eid <- as.character(ncbiID[nbs])
    fm.sub <- featureMatrix.mf[,nbs.eid] ## feature matrix for neighbor set
    mmat <- t(apply(fm.sub,1,cumsum))
    matn <- t(matrix(rep(1:nn,n.singles),ncol=n.singles))
    TPmat <- mmat
    FPmat <- matn - mmat
    FNmat <- Msingles - mmat
    TNmat <- N - TPmat - FPmat - FNmat
    FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
    TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
    ESmat <- TPRmat - FPRmat 
    ES <- ESmat ## ES curves for each matrix 
    Esmax <- apply(ES,1,max) ## maximum ES value of the curve
    indsatmax <- apply(ES,1,which.max) ## index at maximum

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

save(single.features,file=paste(out.dir,"SingleFeatures.RData",sep="/"))

cat(proc.time() - ptm,"\n")

