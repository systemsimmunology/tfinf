
load("mmstore.RData")

Msteps <- seq(20,4000,20)
numMsteps <- length(Msteps)

##N <- length(eids.univ)
N <- 4577

if ( FALSE ){ ## tester
  mms <- mmstore[101:105,1:20,]
  ncube <- array(rep(sort(rep(1:20,5)),1000),dim=c(5,20,1000))
  bigMcube <- array(rep(Msteps[101:105],20*1000),dim=c(5,20,1000))
  mcube <- mms
}

mcube <- mmstore
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

save(ESmaxcubeMgrid,file="ESmaxcubeMgrid.RData")



