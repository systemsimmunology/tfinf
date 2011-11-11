
## If not prepared in an .RData file
source("./initializeBC.R")

## N should be obtained by preamble code (currently initializeBC.R )
## number of ensemblIDs for which we have scans, for which there are unique entrezIDs
## and which are in the expressed set

seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")

Msteps <- seq(20,20*floor(N/20),20)
numMsteps <- length(Msteps)

nearestMpoint <- function ( Mval ){
  index <- round(Mval/20)
  if ( index== 0 ){ index=1}
  return(index)
}

ndraws <- 1000
mcube <- array(0,dim=c(numMsteps,600,ndraws))
 
ptm <- proc.time()
for ( nm in 1:numMsteps ){
  bigM <- Msteps[nm]
  vectosample <- c(rep(TRUE,bigM),rep(FALSE,N-bigM))  
  for ( ns in 1:ndraws ){
    randvec <- sample(vectosample,600)
    mm <- cumsum(randvec)
    mcube[nm,,ns] <- mm
  }
}
proc.time() - ptm

ncube <- array(rep(sort(rep(1:600,numMsteps)),1000),dim=c(numMsteps,600,1000))
bigMcube <- array(rep(Msteps,600*1000),dim=c(numMsteps,600,1000))

## lead astray here, but aperm may come in handy later: array permutation
## usage: aperm( a, c(2,1,3)) permutes (transposes?) first two directions

TPcube <- mcube
FPcube <- ncube - mcube
FNcube <- bigMcube - mcube
rm(mcube,ncube,bigMcube)

TNcube <- N - TPcube - FPcube - FNcube
FPRcube <- FPcube/(FPcube+TNcube) ## False Positive Rate
TPRcube <- TPcube/(TPcube+FNcube) ## True Positive Rate = Recall = Sensitivity ! 
EScube <- TPRcube - FPRcube 
rm(TPcube,FPcube,FNcube,TNcube,TPRcube,FPRcube)
 
ESmaxcubeMgrid <- apply(EScube,c(3,1),cummax)
ESmaxcubeMgrid <- aperm(ESmaxcubeMgrid,c(3,1,2))
rm(EScube) 

save(ESmaxcubeMgrid,file=paste(seq.dir,"ESmaxcubeMgrid.RData",sep="/"))



