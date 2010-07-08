N <- 4687 ## columns in featureMatrix.mf

Msteps <- seq(0,4680,20) ## needed for singles

nearestMpoint <- function ( Mval ){
  index <- round(Mval/20)
  if ( index== 0 ){ index=1}
  return(index)
}

ndraws <- 100
numMsteps <- length(Msteps)
mmstore <- array(0,dim=c(numMsteps,600,ndraws))
 
ptm <- proc.time()
for ( nm in 1:numMsteps ){
  bigM <- Msteps[nm]
  vectosample <- c(rep(TRUE,bigM),rep(FALSE,N-bigM))  
  for ( ns in 1:ndraws ){
    randvec <- sample(vectosample,600)
    mm <- cumsum(randvec)
    mmstore[nm,,ns] <- mm
  }
}
proc.time() - ptm
save(mmstore,file="mmstore.RData")

## Example of evaluation

allMs <- apply(featureMatrix,1,sum)

allMs["CREB-NFKB"]

> allMs["CREB-NFKB"]
CREB-NFKB 
     2536

## index of nearest m
> which(Msteps==2540)
[1] 127
mind <- nearestMpoint(2536)
Mval <- Msteps[mind]

## Il12b
nn
[1] 63


mmat <- t(mmstore[mind,1:nn,])
matn <- t(matrix(rep(1:nn,ndraws),ncol=ndraws))
TPmat <- mmat
FPmat <- matn - mmat
FNmat <- Mval - mmat
TNmat <- N - TPmat - FPmat - FNmat
FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
ESmat <- TPRmat - FPRmat 
ESmaxdistrib <- apply(ESmat,1,max)
Fn <- ecdf(ESmaxdstrib)
pval <- 1-Fn ( whateveritis )


## Same as above, but as a function
getpval <- function (esmaxval,em, nn ){
  mind <- nearestMpoint(em)
  Mval <- Msteps[mind]
  mmat <- t(mmstore[mind,1:nn,])
  matn <- t(matrix(rep(1:nn,ndraws),ncol=ndraws))
  TPmat <- mmat
  FPmat <- matn - mmat
  FNmat <- Mval - mmat
  TNmat <- N - TPmat - FPmat - FNmat
  FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
  TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
  ESmat <- TPRmat - FPRmat 
  ESmaxdistrib <- apply(ESmat,1,max)
  Fn <- ecdf(ESmaxdistrib)
  pval <- 1-Fn (esmaxval)
  return(pval)
}


## Il12b, "CREB-NFKB"
mind <- nearestMpoint(2536)
Mval <- Msteps[mind]
ESmaxdistrib <- ESmaxcubeMgrid[mind,63,]
Fn <- ecdf(ESmaxdistrib)
1 - Fn( value ) 

## Estimate p-value of ESmax value
## Uses "cube" of ESmaxvalues derived from random:
## ESmaxcubeMgrid
ESmaxPval <- function (esmaxval,em, nn ){
  mind <- nearestMpoint(em)
  Mval <- Msteps[mind]
  ESmaxdistrib <- ESmaxcubeMgrid[nn,63,]
  Fn <- ecdf(ESmaxdistrib)
  1 - Fn( value )
  pval <- 1-Fn (esmaxval)
  return(pval)
}


  
