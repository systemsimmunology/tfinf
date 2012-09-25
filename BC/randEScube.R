##
## Construct a distribution of ESmax values for pairs of M and N
##
## Result is returned as a cube with the last dimensin over samplings
##
## If not prepared in an .RData file
source("./initializeBC.R")

maxn <-  600 ## The maximum N that we will consider
ndraws <- 1000 ## number of samples/trials
mst <- 20 # the stepsize on the M grid

## N should be obtained by preamble code (currently initializeBC.R )
## number of ensemblIDs for which we have scans, for which there are unique entrezIDs
## and which are in the expressed set

seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")

## Steps in M, the number of occurences of a features
Msteps <- seq(mst,mst*floor(N/mst),mst)
numMsteps <- length(Msteps)

## construct cube to fill in 
mcube <- array(0,dim=c(numMsteps,maxn,ndraws))
## Grid dimension 1: steps in M
## Grid dimension 2: each single N value up to maxN
## Grid dimnesion 3: one for each sampling
ptm <- proc.time()
for ( nm in 1:numMsteps ){ ## Loop over each possible M in the grid
  bigM <- Msteps[nm]
  vectosample <- c(rep(TRUE,bigM),rep(FALSE,N-bigM)) ## fixed vector of length M with M TRUEs rest FALSE
  for ( ns in 1:ndraws ){
    randvec <- sample(vectosample,maxn) ## randomized version of the above vector
    mm <- cumsum(randvec) ## occurences by cumulative sum
    mcube[nm,,ns] <- mm 
  }
}
proc.time() - ptm

## Cubes with same dimension of N values, and M values, respectively
ncube <- array(rep(sort(rep(1:maxn,numMsteps)),ndraws),dim=c(numMsteps,maxn,ndraws))
bigMcube <- array(rep(Msteps,maxn*ndraws),dim=c(numMsteps,maxn,ndraws))

##  Transform counts to ES object in GSEA
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
## usage: aperm( a, c(2,1,3)) permutes (transposes?) first two directions
## Not sure why aperm is needed, perhaps due to transformation in previous cummax step 
## End result is (still)
## Grid dimension 1: steps in M
## Grid dimension 2: each single N value up to maxN
## Grid dimnesion 3: one for each sampling
rm(EScube) 

save(ESmaxcubeMgrid,file=paste(seq.dir,"ESmaxcubeMgrid.RData",sep="/"))



