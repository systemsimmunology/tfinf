

### Methods for estimating tau, using finite difference approximation for derivative

## 1. Grid of taus
## 2. OLS regression, including negative dotG as a term ( force constraints post fit )

## Could also consider package mgcv, which allows constrained regression
## http://www.stat.columbia.edu/~tzheng/tianblog/2006/07/fitting-constrained-least-square.html


## 2. is running faster than 1. 

##
## Method one, grid of taus
##
estimateBestTau1 <- function (targ,ivars,boost.vec,mus.mat.maxnormed){

  Thalf.grid <- 60*c(1./6.,0.5,1,2,4,10)
  rmsds.per.tau <- numeric(length(Thalf.grid))
  names(rmsds.per.tau) <- as.character(Thalf.grid)
  
  for ( Thalf in Thalf.grid ) {
    tau <- Thalf/log(2.)
    dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,mus.mat.maxnormed)
    ris <- regInfSelectD2F(dataToFit,s=1,zero.offset)
    coefs <- as.numeric(ris[-1])
    offset <- as.numeric(ris[1])
    outpred <- calculateOutpred(dataToFit$inputs,offset,coefs)
    cor <- calcCOR(outpred,dataToFit$outputs)
    rmsd <- calcRMSD(outpred,dataToFit$outputs)
    rmsds.per.tau[as.character(Thalf)] <- rmsd
  }

  best.tau  <- Thalf.grid[which.min(rmsds.per.tau)]/log(2.)
  return( as.numeric(best.tau ))
}


##
## Method 2, include tau as optimized variable
## Followed by 'throw into allowed tau region' if outside
##
estimateBestTau2 <- function (targ,ivars,boost.vec,mus.mat.maxnormed){

  dataToFit <- constructDataToFitTauOpt(t.index.min,t.index.max,targ,ivars,boost.vec,mus.mat.maxnormed)

  if ( length(ivars) > 1 ) {
    df <- as.data.frame(dataToFit)
  } else if ( length(ivars)==1) {
    df <- data.frame(cbind(drop(dataToFit$inputs),dataToFit$outputs))
  }
  
  colnames(df) <- c("NegGDot",ivars,"outputs")
  
  if ( !zero.offset ) {
    lmfit  <- lm(outputs ~ . ,df)
    lmo <- lmfit$coefficients[1]
    lmcs <- lmfit$coefficients[-1]
    names(lmcs) <- c("tau",ivars)
  } else {
    lmfit  <- lm(outputs ~ 0 + . ,df)
    lmo <- 0
    names(lmo) <- "(Intercept)"
    lmcs <- lmfit$coefficients
    names(lmcs) <- c("tau",ivars)
  }


  tau <- lmcs["tau"]
  tau10min <- 10./log(2.)
  tau10hr <- (10*60)/log(2.)
  if ( tau < tau10min ){
    tau <- tau10min
  } else if ( tau > tau10hr ){
    tau <- tau10hr
  }
  
  return(as.numeric(tau))

}


## out ~ G
## in ~ -dG/dt , T1, T2,..
## not to be used with offset, for now

constructDataToFitTauOpt <- function(t.index.min,t.index.max,targ,ivars,boost.vec, mus.mat.maxnormed){

  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  times <- array.times[1: dim(mus.mat.maxnormed)[2] ] ## up to: the maximum avaiable index
  inputs.boosted <- mus.mat.maxnormed[ivars,] ## initialize ( copy dimensions and labels )
  tcourse.out <- mus.mat.maxnormed[targ,]
  if ( is.vector(inputs.boosted) ){
    ## Single inputs 
    save.names <- names(inputs.boosted)
    inputs.boosted <- approx(times,mus.mat.maxnormed[ivars,],times-boost.vec[ivars],rule=2)$y
    names(inputs.boosted) <- save.names
    dataToFit <- makeOutputInputTauOpt(inputs.boosted,tcourse.out,t.index.min,t.index.max)
    ## an extra step to re-instate input names
    colnames(dataToFit$inputs) <- c("NegGDot",ivars)
    ## Possible long-term options
    ## 1. Pass named one-dimensional matrix in place of mus.mat.maxnormed[ivars,]
    ## 2. Pass ivars as argument to makeOuputInputTauOpt
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(times,mus.mat.maxnormed[ivar,],times-boost.vec[ivar],rule=2)$y
    }
    dataToFit <- makeOutputInputTauOpt(inputs.boosted[ivars,],tcourse.out,t.index.min,t.index.max)
  }
  dataToFit
}



makeOutputInputTauOpt <- function (inMat, outVec, t.index.min, t.index.max ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)

  tc.response <- outVec[t.index.min:(t.index.max-1)]
  NegGDot <- -1./(delta.t[t.index.min:(t.index.max-1)]) * (outVec[(t.index.min+1):t.index.max]-outVec[t.index.min:(t.index.max-1)])
  
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,t.index.min:(t.index.max-1)])
    inMatCol <- cbind(NegGDot, inMatCol )
  } else {

    inMatCol <- t(t(inMat[t.index.min:(t.index.max-1)]))

    ## for whatever reason, this creates a column
    
    ##inMatCol <- t(inMat[t.index.min:(t.index.max-1)])

    inMatCol <-  cbind(NegGDot, drop(inMatCol) )
  }

  
  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}
