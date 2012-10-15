
##
##
## For many of these functions t.index.max and t.index.min must be defined globally
##

#
# regInfSelect
# Uses l1ce lasso2 version of lasso
#

regInfSelect <- function(ovar,ivars, tau=30.,s=1.,boost.vec,zero.offset=FALSE,second.order=FALSE) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  tcourse.out <- lps.mus[ovar,]/max(lps.mus[ovar,])
  ##tcourse.out <- lps.mat.max1[ovar,]
  tcourse.in <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(tcourse.in) ){
    tcourse.in <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      tcourse.in[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }

  if ( !second.order ){
    dataToFit <- makeOutputInput.3.lars(tcourse.in,tcourse.out,tau)
  } else {
    dataToFit <- makeOutputInput.secondorder(tcourse.in,tcourse.out,tau)
  }
  df3 <- data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
  colnames(df3) <- c("output",ivars) 

  if ( !zero.offset ){
    l1ce.obj <- l1ce(output ~ .,df3,bound=s, sweep.out= ~1, absolute.t = FALSE)
    coefs <- l1ce.obj$coefficients
  } else { ## zero offset
    l1ce.obj <- l1ce(output ~ 0 + . ,df3,bound=s, absolute.t = FALSE, sweep.out = NULL )
    coefs <- c(0,l1ce.obj$coefficients)
  }
  
  names(coefs) <- c("(Intercept)",ivars)
  coefs
}

regInfSelectD2F <- function(dataToFit,s=1,zero.offset=FALSE) {

  ivars <- colnames(dataToFit$inputs)
  dataFrame <- as.data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
  colnames(dataFrame)[1] <- c("output")

  if ( !zero.offset ){
    l1ce.obj <- l1ce(output ~ .,dataFrame,bound=s, sweep.out= ~1, absolute.t = FALSE)
    coefs <- l1ce.obj$coefficients
  } else { ## zero offset
    l1ce.obj <- l1ce(output ~ 0 + . ,dataFrame,bound=s, absolute.t = FALSE, sweep.out = NULL )
    coefs <- c(0,l1ce.obj$coefficients)
  }
  
  names(coefs) <- c("(Intercept)",ivars)
  coefs
}


## minimum s by gcv
## warning, No check for 'simple convex curve'
## Take one step left of minimum. We do not know the variability in gcv currently
## otherwise, would use same protocol as with Halo inferelator
## type = c("OPT", "Tibshirani")
getBestS <- function(dataToFit,zero.offset=FALSE,type = "OPT" ){
  
  dataFrame <- as.data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
  colnames(dataFrame)[1] <- c("output")
  ## ivars should be preserved as column names
  approx.zero = 0.0001
  sseq=c(approx.zero,seq(0.1,1,0.1)) ## sequence of shrinkage parameters. Exact 0 gives an error in gcv.l1ce
  
  if ( !zero.offset ){
    gcv.matrix <- gcv(l1ce(output ~ .,dataFrame,bound=sseq, sweep.out= ~1, absolute.t = FALSE),type=type)
  } else {
    gcv.matrix <- gcv(l1ce(output ~ 0 + .,dataFrame,bound=sseq, sweep.out=NULL, absolute.t = FALSE),type=type)
  }
  min.index <- which.min(gcv.matrix[,"gcv"])
  s.index <- min.index -1
  if ( s.index == 0 ){ s.index <- 1 } ## Negative s is meaningless
  rel.bound.min.gcv <- gcv.matrix[s.index,"rel.bound"]
  gcv.value <- gcv.matrix[s.index,"gcv"]
  if ( rel.bound.min.gcv == approx.zero ){ rel.bound.min.gcv <- 0} ## exact 0 works fine in l1ce 
  returnList <- list()
  returnList$s <- rel.bound.min.gcv
  returnList$gcv <- gcv.value
  return(returnList)
}

## Refit dataToFit to OLS.
## The input value s is the lasso value that defines the model complexity 
## s is not used for actual model fit
modReFit <- function(dataToFit,s=1,zero.offset=FALSE,targ,t.index.min,t.index.max,tau){

  ## model, to be returned
  returnMod <- list()
  
  ## orignal data frame (unselected)
  dataFrame <- as.data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
  colnames(dataFrame)[1] <- c("output")
  ivars.orig <- colnames(dataToFit$inputs) ## ivars in the original data frame, unselected
  
  ## First, determine which inputs survived
  if ( !zero.offset ){
    coefs.post.lasso <- coefficients(l1ce(output ~.,dataFrame,bound=s,sweep.out=~1,absolute.t = FALSE))[-1]
  } else {
    coefs.post.lasso <- coefficients(l1ce(output ~0+.,dataFrame,bound=s,sweep.out=NULL,absolute.t = FALSE))
  }
  
  ivars.keepers <- as.character(ivars.orig[which(coefs.post.lasso != 0 )])
  ## After finding which variables are kept, s is no longer used
  
  ## 2nd , refit parameters
  ## Note the bound=1 below, meaning OLS
  if ( length(ivars.keepers) > 0 ){
  
    dataToFit.inputs.new <- drop(dataToFit$inputs[,ivars.keepers]) ## may need cases for length(ivars.keepers)
    dataFrame.new <- as.data.frame(cbind(dataToFit$outputs,dataToFit.inputs.new )) ## may need to investigate headings etc. 
    colnames(dataFrame.new) <- c("output",ivars.keepers)
    if ( !zero.offset ){
      ris <- coefficients(l1ce(output ~.,dataFrame.new,bound=1,sweep.out=~1,absolute.t = FALSE))
      betas <- as.numeric(ris[-1])
      offset <- as.numeric(ris[1]) 
    } else {
      ris <- coefficients(l1ce(output ~0+.,dataFrame.new,bound=1,sweep.out=NULL,absolute.t = FALSE))
      betas <- as.numeric(ris)
      offset <- 0. 
    }
    
    ## Quality measures
    outpred <- calculateOutpred(dataToFit.inputs.new,offset,betas) # should work for 1 or two keepers
    cor.val <- calcCOR(outpred,dataToFit$outputs)
    rmsd.val <- calcRMSD(outpred,dataToFit$outputs)

    ## targ is not known!!, nor is t.index.min,t.index.max,tau
    
    start.val <- lps.mat.max1[targ,t.index.min] ## The RHS unfortunately contains a global variable
    scale.fac <- max.intensity[targ] ## The RHS unfortunately contains a global variable
    wt.predicted.fullscale <- calculateTimeEvolPrediction.6(t.index.min,t.index.max,tau,outpred,start.val,scale.fac)



    cor.unrolled <- calcCOR(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max]) ## The RHS unfortunately contains a global variable
    rmsd.unrolled <- calcRMSD(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max]) ## The RHS unfortunately contains a global variable
    
    
    returnMod$tfs <- ivars.keepers
    returnMod$tfs.cname <- as.character(cname.compare[ivars.keepers])
    returnMod$betas <- betas
    returnMod$offset <- offset
    returnMod$rmsd <- rmsd.val
    returnMod$cor <- cor.val
    returnMod$rmsd.unrolled <- rmsd.unrolled
    returnMod$cor.unrolled <- cor.unrolled 
  } else {
    returnMod <- NULL
  }
  
  returnMod
}

#
# regInfSelect
# Uses lars version of lasso
#

regInfSelect.old <- function(ovar,ivars, tau=30.,s=1.,boost.vec) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  tcourse.out <- lps.mus[ovar,]/max(lps.mus[ovar,])
  ##tcourse.out <- lps.mat.max1[ovar,]
  tcourse.in <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(tcourse.in) ){
    tcourse.in <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      tcourse.in[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }

  if ( !is.null(dim(tcourse.in)) ){
    dataToFit <- makeOutputInput.3.lars(tcourse.in,tcourse.out,tau)
    object <- lars(dataToFit$inputs,dataToFit$outputs,type="lasso")
    coefs <- coef.lars(object,s=s,mode="fraction")
    subvec <- apply(dataToFit$inputs,2,mean)
    offset.shift <- -sum(coefs*subvec) 
    offset <- mean(dataToFit$outputs)+offset.shift
    coefs <- c(offset,coefs)
    names(coefs) <- c("(Intercept)",ivars)
  } else { ## seems l1ce does not fail here!
    dataToFit <- makeOutputInput.3.lars(tcourse.in,tcourse.out,tau) 
    df3 <- data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
    l1ce.obj <- l1ce(X1~X2,df3,bound=s, absolute.t = FALSE)
    coefs <- l1ce.obj$coefficients
    
    ## zero offset
    ##l1ce.obj <- l1ce(X1 ~ 0 + . ,df,bound=s, sweep.out = NULL, absolute.t = FALSE, sweep.out = NULL)
    ##coefs <- c(0,l1ce.obj$coefficients)

    names(coefs) <- c("(Intercept)",ivars)
  }
  coefs
}

##
## Squashed version of regInfSelect
## 
regInfSelect.4 <- function(ovar,ivars, tau=30.,s=1.,boost.vec) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  tcourse.out <- lps.mat.max1[ovar,]
  tcourse.in <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(tcourse.in) ){
    tcourse.in <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      tcourse.in[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }
  outputs <- t.course.out
  ofixed <- mean(outputs[1:2])
  outsubo  <- outputs - ofixed
  
  if ( !is.null(dim(tcourse.in)) ){
    dataToFit <- makeOutputInput.4.lars(tcourse.in,tcourse.out,tau) ## tcourse.in is the unsquashed timecourse
    object <- lars(dataToFit$inputs,outsubo,type="lasso") ## not totally sure whether intercept is handled properly
    coefs <- coef.lars(object,s=s,mode="fraction")
  } else { ## seems l1ce does not fail here!
    dataToFit <- makeOutputInput.4.lars(tcourse.in,tcourse.out,tau) 
    df3 <- data.frame(cbind(outsubo,drop(dataToFit$inputs)))
    l1ce.obj <- l1ce(X1~0+X2,df3,bound=s, absolute.t = FALSE)
    coefs <- l1ce.obj$coefficients
    names(coefs) <- ivars
  }
  coefs
}




regInfSelect.old <- function(ovar,ivars, tau=30.,s=1.,target.pull=0) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  tcourse.out <- lps.mat.max1[ovar,t.index.min:t.index.max]
  tcourse.out <- approx(array.times,lps.mat.max1[targ,],array.times+target.pull,rule=2)$y
  tcourse.in <- lps.mat.max1[ivars,t.index.min:t.index.max]  ## input values maximized at 1 to have comparable weights

  if ( !is.null(dim(tcourse.in)) ){
    dataToFit <- makeOutputInput.3.lars(tcourse.in,tcourse.out,tau)
    object <- lars(dataToFit$inputs,dataToFit$outputs,type="lasso")
    coefs <- coef.lars(object,s=s,mode="fraction")
  } else { ## seems l1ce does not fail here!
    dataToFit <- makeOutputInput.3.lars(tcourse.in,tcourse.out,tau) 
    df3 <- data.frame(cbind(dataToFit$outputs,drop(dataToFit$inputs)))
    l1ce.obj <- l1ce(X1~X2,df3,bound=s, absolute.t = FALSE)
    coefs <- l1ce.obj$coefficients
    names(coefs) <- c("(Intercept)",ivars)
  }
  coefs
}


##
## Calculate RMSD for transformed variables
##
findRMSD <- function(coefs,targ,ivars,tau=30.,s=1.,boostvec,zero.offset=FALSE) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  inputs.boosted <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }
  tcourse.out <- lps.mus[targ,]/max(lps.mus[targ,])
  dataToFit <- makeOutputInput.3.lars(inputs.boosted,tcourse.out,tau)
  
  ris <- regInfSelect(targ,ivars,tau,s,boost.vec,zero.offset) ## This is crazy. Shouldn't have to re-predict.
  coefs <- as.numeric(ris[-1])
  pred.pids <- names(ris[-1])
  offset <- as.numeric(ris[1])
  
  if (!(dim(dataToFit$inputs)[1]==1)){
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- offset + drop(dataToFit$inputs) * coefs ## this expression is specific to l1ce fit
  }
  return (sqrt(sum((outpred-dataToFit$outputs)^2))/length(outpred))
}

##
## Calculate RMSD for transformed variables
##
findRMSD.2 <- function(outpred,dataToFit){
  return (sqrt(sum((outpred-dataToFit$outputs)^2))/length(outpred))
}

## 
##
calcRMSD <- function(predicted,measured){
  return (sqrt(sum((predicted-measured)^2)/length(predicted)))
  ## Jan 2007: Have used the one below for a while but the 'N' should be under the root
  ## return (sqrt(sum((predicted-measured)^2))/length(predicted))
}

##
## Calculate correlation coefficient of transformed variables
##
findCOR <- function(coefs,targ,ivars,tau=30.,s=1.,boostvec,zero.offset=FALSE) {
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  inputs.boosted <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }
  tcourse.out <- lps.mus[targ,]/max(lps.mus[targ,])
  dataToFit <- makeOutputInput.3.lars(inputs.boosted,tcourse.out,tau)

  ris <- regInfSelect(targ,ivars,tau,s,boost.vec,zero.offset) ## This is crazy. Shouldn't have to re-predict.
  coefs <- as.numeric(ris[-1])
  pred.pids <- names(ris[-1])
  offset <- as.numeric(ris[1])
  
  if (!(dim(dataToFit$inputs)[1]==1)){
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- offset + drop(dataToFit$inputs) * coefs ## this expression is specific to l1ce fit
  }
  
  return (cor(outpred,dataToFit$outputs))

}

## 
##
calcCOR <- function(predicted,measured){
  return (cor(predicted,measured))
}

##
## Calculate correlation coefficient of transformed variables
##
findCOR.2 <- function(outpred,dataToFit){
  return (cor(outpred,dataToFit$outputs))
}


findRMSD.old <- function(coefs,targ,ivars,tau=30.,target.pull=0) { # older version with target pull
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  targ.exp.pulled <- approx(array.times,lps.mat.max1[targ,],array.times+target.pull,rule=2)$y
  dataToFit <- makeOutputInput.3.lars(lps.mat.max1[ivars,t.index.min:t.index.max],targ.exp.pulled[t.index.min:t.index.max],tau)

  if (!(dim(dataToFit$inputs)[1]==1)){
    subvec <- apply(dataToFit$inputs,2,mean)
    offset.shift <- -sum(coefs*subvec) 
    offset <- mean(dataToFit$outputs)+offset.shift
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- coefs[1] + dataToFit$inputs * coefs[2] ## this expression is specific to l1ce fit
  }
  return (sqrt(sum((outpred-dataToFit$outputs)^2))/length(outpred))
}

#
# Construct response data frame
# inMat should be a row-wise slice
#

makeOutputInput.lars <- function (inMat, outVec, tau, squash.mode = "none" ) {

  delta.t <- c(20.,20.,20.,20.,40.)
  tc.response <- tau/delta.t * (outVec[2:6]-outVec[1:5])+outVec[1:5]
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,1:5])
  } else {
    inMatCol <- t(inMat[1:5])
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}


##
##  This version for the long lps time-course
##

makeOutputInput.2.lars <- function (inMat, outVec, tau, squash.mode = "none" ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,240.,960.)
  
  tc.response <- tau/(delta.t[1:(t.index.max-1)]) * (outVec[2:t.index.max]-outVec[1:(t.index.max-1)])+outVec[1:(t.index.max-1)]
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,1:(t.index.max-1)])
  } else {
    inMatCol <- t(inMat[1:(t.index.max-1)])
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}

##
##  This version for the long lps time-course, including hr 6, hr 18
##

makeOutputInput.3.lars <- function (inMat, outVec, tau, squash.mode = "none" ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  
  tc.response <- tau/(delta.t[t.index.min:(t.index.max-1)]) * (outVec[(t.index.min+1):t.index.max]-outVec[t.index.min:(t.index.max-1)])+outVec[t.index.min:(t.index.max-1)]
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,t.index.min:(t.index.max-1)])
  } else {
    inMatCol <- t(inMat[t.index.min:(t.index.max-1)])
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}

##
## 2nd order approximation to derivative 
## 
## t.index.min t.index assumed known

makeOutputInput.secondorder <- function (inMat, outVec, tau ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)

  hminus <- delta.t[t.index.min:(t.index.max-2)]
  hplus <- delta.t[(t.index.min+1):(t.index.max-1)]
  Gnp1 <- outVec[(t.index.min+2):t.index.max]
  Gn <- outVec[(t.index.min+1):(t.index.max-1)]
  Gnm1 <-  outVec[t.index.min:(t.index.max-2)]
  
  Gdot <- ( hminus^2 * Gnp1 + (hplus^2-hminus^2)* Gn - hplus^2 * Gnm1 ) /
    hplus/hminus/(hplus+hminus)

  names(Gdot) <- conds[t.index.min:(t.index.max-2)]

  ## Use finite difference for last point
  ##GdotFinal <- (outVec[t.index.max]-outVec[t.index.min])/delta.t[(t.index.max-1)]
  ##names(GdotFinal) <- conds[t.index.max-1]
  ##Gdot <- c(Gdot,GdotFinal)

  tc.response <- tau*Gdot+ outVec[(t.index.min+1):(t.index.max-1)]

  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,(t.index.min+1):(t.index.max-1)])
  } else {
    inMatCol <- t(inMat[(t.index.min+1):(t.index.max-1)])
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}


g<-function(x){1/(1+exp(-x))}
ginv<-function(y){log(y/(1-y))}
##
## Squashed version of version .3.
##
makeOutputInput.4.lars <- function (inMat, outVec, tau, squash.mode = "squash" ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  
  tc.response <- tau/(delta.t[t.index.min:(t.index.max-1)]) * (outVec[(t.index.min+1):t.index.max]-outVec[t.index.min:(t.index.max-1)])+outVec[t.index.min:(t.index.max-1)]

  ## Before applying ginv, we need this vector to be in the range (0,1)
  v1 <- tc.response
  v2 <- v1-min(v1)
  v3 <- v2/max(v2)
  eps <- 0.01
  v4 <- (v3-0.5)*(1-2*eps)+0.5
  v5 <- ginv(v4) 
  outputs.squished <- v5
  
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,t.index.min:(t.index.max-1)])
  } else {
    inMatCol <- t(inMat[t.index.min:(t.index.max-1)])
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- outputs.squished
  
  inOut
  
}

## Same as .3, but with t.index.min and t.index.max as input variables
makeOutputInput.5.lars <- function (inMat, outVec, tau, t.index.min, t.index.max ) {

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  
  tc.response <- tau/(delta.t[t.index.min:(t.index.max-1)]) * (outVec[(t.index.min+1):t.index.max]-outVec[t.index.min:(t.index.max-1)])+outVec[t.index.min:(t.index.max-1)]
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,t.index.min:(t.index.max-1)])
  } else {
    inMatCol <- t(t(inMat[t.index.min:(t.index.max-1)]))
    ## for whatever reason, this creates a column
  }

  inOut <- list()
  inOut$inputs <- inMatCol
  inOut$outputs <- tc.response
  
  inOut
  
}

#
# regInfSelect.0
#
# Lasso2 version had weird problem, probably due to the fourmoola construct

regInfSelect.0 <- function(outvar,invars) {

  outvar.gname <- gid2gname(outvar)
  invars.gname <- gid2gname(invars)

  tc.out <- cube["LPS",outvar,]
  tc.in <- cube["LPS",invars,]
  
  dfToFit <- makeOutputInput.lasso(tc.in,tc.out,tau)

  ## Reinstate output variable name
  ## Convert to common names, as probeset column names are seen to cause problem for lm etc.
  colnames(dfToFit) <- c(outvar.gname, invars.gname)
  
  fourmoola  <- as.formula(paste(outvar.gname,"~","."))

  lasso.result <- l1ce(fourmoola, dfToFit,bound = 0.5, absolute.t = FALSE)

  lasso.result$coefficients 
  
}








#
# Construct response data frame
# inMat should be a row-wise slice
#

makeOutputInput.lasso <- function (inMat, outVec, tau, squash.mode = "none" ) {

  delta.t <- c(20.,20.,20.,20.,40.)
  tc.response <- tau/delta.t * (outVec[2:6]-outVec[1:5])+outVec[1:5]
  if ( ! is.null(dim(inMat)) ) { ## test for the collapsing matrix 
    inMatCol <-  t(inMat[,1:5])
  } else {
    inMatCol <- t(inMat[1:5])
  }
  outVec <- tc.response
  m2 <- t(rbind(tc.response,inMatCol))
  ##colnames(m2) <- c("Y",gid2gname(rownames(inMat))) # we've lost the output variable name
  df2 <- as.data.frame(m2)

  df2
  
}


#
# return list of incoming TFs
# as a list of character arrays
#
# input is a paired set of tfs and targets 

getIncomingWires <- function(tfs, targs){

  targs.unique <- unique(sort(targs))

  outlist <- list()

  for ( targ in targs.unique ){

    
    tinds <- which( targs == targ )
    outlist[[targ]] <- sort(unique(tfs[tinds]))

  }

  outlist

}

#
# return list of TF targets
# as a list of character arrays
#
# input is a paired set of tfs and targets 

getOutgoingWires <- function(tfs, targs){

  tfs.unique <- unique(sort(tfs))
  outlist <- list()

  for ( tf in tfs.unique ){
    
    tinds <- which( tfs == tf )
    outlist[[tf]] <- sort(unique(targs[tinds]))
  }
  outlist

}



##---------------------------------------------------------

plotPrediction <- function(target.pid,pred.pids,pred.out){

  profile.true <- lps.mat[target.pid,1:6] ################ Note indexing!!!!!!!!!!!!!

  profile.pred <- 0.5 * max.intensity[target.pid] * lps.mat.max1[pred.pids,]
  
  ymax <- max(c(pred.out, profile.true))
  xvals <- c(0.,20.,40.,60.,80.,120.)

  main=gid2gname(target.pid)

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  plot(xvals,profile.true, ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)

  axis(1)
  axis(2)

  lines(c(20.,40.,60.,80.,120.) ,pred.out, col='blue',lty=2,lwd=3)
  ##lines(c(20.,40.,60.,80.,120.) ,pred.out.drop1, col='green',lty=2,lwd=3)


  ##  colorvec = c('red','magenta','pink','cyan','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(pred.pids))
  lty <- rep(lty, length.out=length(pred.pids))

  
  for ( i in 1:length(pred.pids) ){
    
    lines( xvals, profile.pred[i,], col=col[i],lty=lty[i],lwd=3)

  }

  ## Color and linetype

  ## Color and line-type indexes

  legend <- gid2gname(pred.pids)

  
  lg <- c(1.5,0.9*max.intensity[target.pid] )
  ## lg <- c(n.col-2.5,1.2 )
  legend(lg[1],lg[2], legend=legend, col=col, lty=lty)

  

  par(op)

}
##---------------------------------------------------------
## Longer time version of plotPrediction
##
## pred.pids = pids of predicted influencors
## pred.out = predicted *profile* of output gene

plotPrediction.2 <- function(target.pid,pred.pids,pred.out){

  profile.true <- lps.mus[target.pid,t.index.min:t.index.max]

  ## rescale profile of predicted TFs, for plotting purposes
  profile.pred <- 0.5 * max.intensity[target.pid] * lps.mat.max1[pred.pids,t.index.min:t.index.max]
  
  ymax <- max(c(pred.out, profile.true))
  xvals <- c(0.,20.,40.,60.,80.,120.,240.,480.,1440.)[t.index.min:t.index.max]

  main=paste(gid2gname(target.pid),":",target.pid)

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  plot(xvals,profile.true[t.index.min:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)

  axis(1)
  axis(2)

  lines(c(20.,40.,60.,80.,120.,240.,480.,1440.) ,pred.out, col='blue',lty=2,lwd=3)
  ##lines(c(20.,40.,60.,80.,120.) ,pred.out.drop1, col='green',lty=2,lwd=3)

  colorvec = c('red','magenta','pink','cyan','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(pred.pids))
  lty <- rep(lty, length.out=length(pred.pids))

  
  for ( i in 1:length(pred.pids) ){
    lines( xvals, profile.pred[i,t.index.min:t.index.max], col=colorvec[i],lty=1,lwd=3)
#    lines( xvals, profile.pred[i,], col=col[i],lty=lty[i],lwd=3)

  }

  ## Color and linetype

  ## Color and line-type indexes

  legend <- gid2gname(pred.pids)

  
  lg <- c(1.5,0.9*max.intensity[target.pid] )
  ## lg <- c(n.col-2.5,1.2 )
  ##legend(lg[1],lg[2], legend=legend, col=col, lty=lty, lwd=3)
  legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=3)
  

  par(op)

}

##---------------------------------------------------------
## Longer time version of plotPrediction
##
## pred.pids = pids of predicted influencors
## pred.out = predicted *profile* of output gene

## target.error. Whether or not to plot standard errors on the measured expressed
## pred.out.error. Can contain predictions from randomization trials

plotPrediction.3 <- function(target.pid,pred.coefs,pred.out,target.error=FALSE,pred.out.error=NULL){

  pred.pids <- names(pred.coefs)
  coefs <- as.numeric(pred.coefs)
  
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  
  profile.true <- lps.mus[target.pid,1:t.index.max] ## we'll plot the initial values, even if they were not use for the fit

  if ( target.error ){
    error.true <- lps.muStdErrs[target.pid,1:t.index.max] ## the measurement error of the 'true' curve
  }
  
  ## rescale profile of predicted TFs (scale determined by target intensity ) for plotting purposes
  profile.pred <- 0.5 * max.intensity[target.pid] * lps.mat.max1[pred.pids,1:t.index.max]
  
  ymax <- max(c(profile.true,pred.out,as.vector(pred.out.error)))
  if ( target.error ){
    ymax <- max(c(profile.true+error.true,pred.out,as.vector(pred.out.error)))
  }
  ## seems like NULL does not affect max


  xvals <- 1:t.index.max
  
  main=paste(gid2gname(target.pid),":",target.pid)

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  
  if ( !target.error ){
    plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',main=main,lwd=3,axes=FALSE)
  } else {
    plotCI(x=xvals,y=profile.true[1:t.index.max],uiw=error.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',main=main,lwd=3,axes=FALSE)
  }
    
  axislabels <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  axis(1,1:t.index.max,labels=axislabels[1:t.index.max])
  axis(2)

  ## Plot randomization trials if needed
  if ( !is.null(pred.out.error) ){
    n.rands <- nrow(pred.out.error)
    for ( n in 1:n.rands ){
      lines(xvals[(t.index.min+1):t.index.max],pred.out.error[n,], col='#93D6FF',lty=1,lwd=1)
    }
    ## Yuk. We must redraw the plots we drew over 
    if ( !target.error ){
      lines(xvals,profile.true[1:t.index.max], col='blue',lwd=3)
    } else {
      plotCI(x=xvals,y=profile.true[1:t.index.max],uiw=error.true[1:t.index.max], col='blue',type='l',add=TRUE)
    } 
  }
 
  
  lines(xvals[(t.index.min+1):t.index.max],pred.out, col='blue',lty=2,lwd=3)

  colorvec = c('red','magenta','pink','cyan','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(pred.pids))
  lty <- rep(lty, length.out=length(pred.pids))

  if ( length(pred.pids) > 1 ){
    for ( i in 1:length(pred.pids) ){
      lines( xvals, profile.pred[i,1:t.index.max], col=colorvec[i],lty=1,lwd=3)
    }
  }
  if ( length(pred.pids) == 1 ){
    lines( xvals, profile.pred[1:t.index.max], col=colorvec[1],lty=1,lwd=3)
  }

  s1 <- gid2gname(pred.pids)
  s2 <- as.character(round(coefs,digits=1))
  s2 <- paste("(",s2,")",sep="")
  legend <- paste(s1,s2)
  
  lg <- c(1.5,0.9*max.intensity[target.pid] )
  legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=3)
  

  par(op)

}


##---------------------------------------------------------

calculatePrediction <- function( target.pid, tau ){

  pred.pids <- names(collect[[target.pid]])
  coefs <- as.numeric(collect[[target.pid]])

  
  tc.out <- lps.mat.max1[target.pid,]
  tc.in <- lps.mat.max1[pred.pids,]

  dataToFit <- makeOutputInput.lars(tc.in,tc.out,tau)
  subvec <- apply(dataToFit$inputs,2,mean)
  offset.shift <- -sum(coefs*subvec)
  offset <- mean(dataToFit$outputs)+offset.shift

  outpred <- offset + drop(dataToFit$inputs %*% coefs)

  outpred.unrolled <- delta.t / tau * outpred + ( 1 - delta.t/tau) * dataToFit$outputs[1:5]

  pred.out <- max.intensity[target.pid] * outpred.unrolled
  names(pred.out) <- c("min20","min40","min60","min80","min120")

  pred.out
}



##---------------------------------------------------------
## LPS long time-course

calculatePrediction.2 <- function( regInfSelectResult, target.pid, tau ){

  pred.pids <- names(regInfSelectResult)
  coefs <- as.numeric(regInfSelectResult)

  tc.out <- lps.mat.max1[target.pid,t.index.min:t.index.max]
  tc.in <- lps.mat.max1[pred.pids,t.index.min:t.index.max]

  dataToFit <- makeOutputInput.2.lars(tc.in,tc.out,tau)
  subvec <- apply(dataToFit$inputs,2,mean)
  offset.shift <- -sum(coefs*subvec)
  offset <- mean(dataToFit$outputs)+offset.shift

  outpred <- offset + drop(dataToFit$inputs %*% coefs)
  delta.t <- c(20.,20.,20.,20.,40.,120.,240.,960.)
  outpred.unrolled <- delta.t[1:(t.index.max-1)] / tau * outpred + ( 1 - delta.t[1:(t.index.max-1)]/tau) * dataToFit$outputs[1:(t.index.max-1)]

  pred.out <- max.intensity[target.pid] * outpred.unrolled
  names(pred.out) <- c("min20","min40","min60","min80","min120","hr4","hr8","hr24")[1:(t.index.max-1)]

  pred.out
}

##---------------------------------------------------------
## LPS long time-course, including hr6, hr8

calculatePrediction.3 <- function( regInfSelectResult, target.pid, tau ){

  pred.pids <- names(regInfSelectResult)
  coefs <- as.numeric(regInfSelectResult)

  tc.out <- lps.mat.max1[target.pid,t.index.min:t.index.max]
  tc.in <- lps.mat.max1[pred.pids,t.index.min:t.index.max]

  dataToFit <- makeOutputInput.3.lars(tc.in,tc.out,tau)
  subvec <- apply(dataToFit$inputs,2,mean)
  offset.shift <- -sum(coefs*subvec)
  offset <- mean(dataToFit$outputs)+offset.shift

  outpred <- offset + drop(dataToFit$inputs %*% coefs)
  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  outpred.unrolled <- delta.t[t.index.min:(t.index.max-1)] / tau * outpred + ( 1 - delta.t[t.index.min:(t.index.max-1)]/tau) * dataToFit$outputs[t.index.min:(t.index.max-1)]

  pred.out <- max.intensity[target.pid] * outpred.unrolled
  names(pred.out) <- c("min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")[t.index.min:(t.index.max-1)]

  pred.out
}

##---------------------------------------------------------
## LPS full time evolution
##---------------------------------------------------------

calculateTimeEvolPrediction.3 <- function(targ,ivars,tau,s=1.,boost.vec, zero.offset = FALSE){

  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)

  inputs.boosted <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  tcourse.out <- lps.mus[targ,]/max(lps.mus[targ,])
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
    dataToFit <- makeOutputInput.3.lars(inputs.boosted,tcourse.out,tau)
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
    dataToFit <- makeOutputInput.3.lars(inputs.boosted[ivars,],tcourse.out,tau)
  }


  ris <- regInfSelect(targ,ivars,tau,s,boost.vec,zero.offset) ## This is crazy. Shouldn't have to re-predict.
  coefs <- as.numeric(ris[-1])
  pred.pids <- names(ris[-1])
  offset <- as.numeric(ris[1])
  if (!(dim(dataToFit$inputs)[1]==1)){
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- offset + dataToFit$inputs * coefs
  }

  ## (redundant?) preamble over 
  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  term1 <- delta.t[t.index.min:(t.index.max-1)] / tau * outpred ## vector for the first term
  term2.multiplier <-( 1 - delta.t[t.index.min:(t.index.max-1)]/tau)
  conds <- c("min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  wt.predicted <- numeric(length=(t.index.max-t.index.min))
  names(wt.predicted) <- conds[t.index.min:(t.index.max-1)]


  wt.predicted[1] <- term1[1] + term2.multiplier[1] * tcourse.out[t.index.min] ## not sure about the last here
  for ( i in 2:(t.index.max-t.index.min) ){
    wt.predicted[i] <- term1[i] + term2.multiplier[i] * wt.predicted[i-1]
  }
  
  wt.predicted.fullscale <- max(lps.mus[targ,]) * wt.predicted

  return(wt.predicted.fullscale)
}

##
## consructDataToFit <- 
## 
constructDataToFit <- function(t.index.min,t.index.max,tau,targ,ivars,boost.vec, mus.mat.maxnormed){
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  times <- array.times[1: dim(mus.mat.maxnormed)[2] ] ## up to: the maximum avaiable index
  inputs.boosted <- mus.mat.maxnormed[ivars,] ## initialize ( copy dimensions and labels )
  tcourse.out <- mus.mat.maxnormed[targ,]
  if ( is.vector(inputs.boosted) ){
    ## Single inputs
    save.names <- names(inputs.boosted)
    inputs.boosted <- approx(times,mus.mat.maxnormed[ivars,],times-boost.vec[ivars],rule=2)$y
    names(inputs.boosted) <- save.names
    dataToFit <- makeOutputInput.5.lars(inputs.boosted,tcourse.out,tau,t.index.min,t.index.max)
    ## an extra step to re-instate input names
    colnames(dataToFit$inputs) <- ivars
    ## Possible long-term options
    ## 1. Pass named one-dimensional matrix in place of mus.mat.maxnormed[ivars,]
    ## 2. Pass ivars as argument to makeOuputInput.5.lars
  } else {
    ## Multiple inputs 
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(times,mus.mat.maxnormed[ivar,],times-boost.vec[ivar],rule=2)$y
    }
    dataToFit <- makeOutputInput.5.lars(inputs.boosted[ivars,],tcourse.out,tau,t.index.min,t.index.max)
  }
  dataToFit
}

##
## calculateOutpred
##
## inputs is typically dataToFit$inputs
calculateOutpred <- function(inputs,offset,coefs){
  if ( length(coefs)>1 ){
    outpred <- offset + drop(inputs %*% coefs)
  } else {
    outpred <- offset + inputs * coefs
  }
  drop(outpred)
}
   

##---------------------------------------------------------
## LPS full time evolution
## This version has offset and betas as input
##---------------------------------------------------------


calculateTimeEvolPrediction.5 <- function(t.index.min,t.index.max,coefs,offset,targ,ivars,tau,s=1.,boost.vec, mus.mat.maxnormed){

  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  times <- array.times[1: dim(mus.mat.maxnormed)[2] ] ## up to: the maximum avaiable index
  
  inputs.boosted <- mus.mat.maxnormed[ivars,] ## initialize ( copy dimensions and labels )

  tcourse.out <- mus.mat.maxnormed[targ,]
  
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(times,mus.mat.maxnormed[ivars,],times-boost.vec[ivars],rule=2)$y
    dataToFit <- makeOutputInput.5.lars(inputs.boosted,tcourse.out,tau,t.index.min,t.index.max)
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(times,mus.mat.maxnormed[ivar,],times-boost.vec[ivar],rule=2)$y
    }
    dataToFit <- makeOutputInput.5.lars(inputs.boosted[ivars,],tcourse.out,tau,t.index.min,t.index.max)
  }

  if ( length(coefs)>1 ){
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- offset + dataToFit$inputs * coefs
  }

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  term1 <- delta.t[t.index.min:(t.index.max-1)] / tau * outpred ## vector for the first term
  term2.multiplier <-( 1 - delta.t[t.index.min:(t.index.max-1)]/tau)
  conds <- c("min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  wt.predicted <- numeric(length=(t.index.max-t.index.min))
  names(wt.predicted) <- conds[t.index.min:(t.index.max-1)]


  wt.predicted[1] <- term1[1] + term2.multiplier[1] * tcourse.out[t.index.min] ## not sure about the last here
  for ( i in 2:(t.index.max-t.index.min) ){
    wt.predicted[i] <- term1[i] + term2.multiplier[i] * wt.predicted[i-1]
  }

  ## The universal 'scale' for the profiles is the maximum lps intensity, max.intensity
  wt.predicted.fullscale <- max.intensity[targ] * wt.predicted

  return(wt.predicted.fullscale)

}

## start.val is the initial value, typically tcourse.out[t.index.min] 
calculateTimeEvolPrediction.6 <- function(t.index.min,t.index.max,tau,outpred,start.val,scale.fac){

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  term1 <- delta.t[t.index.min:(t.index.max-1)] / tau * outpred ## vector for the first term
  term2.multiplier <-( 1 - delta.t[t.index.min:(t.index.max-1)]/tau)
  conds <- c("min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  wt.predicted <- numeric(length=(t.index.max-t.index.min))
  names(wt.predicted) <- conds[t.index.min:(t.index.max-1)]

  wt.predicted[1] <- term1[1] + term2.multiplier[1] * start.val 
  for ( i in 2:(t.index.max-t.index.min) ){
    wt.predicted[i] <- term1[i] + term2.multiplier[i] * wt.predicted[i-1]
  }

  ## The universal 'scale' for the profiles is the maximum lps intensity, max.intensity
  wt.predicted.fullscale <- scale.fac * wt.predicted
  
  return(wt.predicted.fullscale)

}



## For a given TF probeset, find which targets it is a candidate for
findGenesWithTFBS <- function (tfpsoi) {
  targ.collection <- vector()
  for ( targ in targs.expressed ){
     if ( tfpsoi %in% cands[[targ]] ){
      targ.collection <- c(targ.collection,targ)
    }
  }
  return (targ.collection)
}











calculateTimeEvolPrediction.3.old <- function(targ,ivars,tau,target.pull=0 ){ ## target pull version
  
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  targ.exp.pulled <- approx(array.times,lps.mat.max1[targ,],array.times+target.pull,rule=2)$y
  dataToFit <- makeOutputInput.3.lars(lps.mat.max1[ivars,t.index.min:t.index.max],targ.exp.pulled[t.index.min:t.index.max],tau)

  regInfSelectResult <- regInfSelect(targ,ivars,tau,s=1.,target.pull) ## This is crazy. Shouldn't have to re-predict.
  pred.pids <- names(regInfSelectResult)
  coefs <- as.numeric(regInfSelectResult)

  if (!(dim(dataToFit$inputs)[1]==1)){
    subvec <- apply(dataToFit$inputs,2,mean)
    offset.shift <- -sum(coefs*subvec)
    offset <- mean(dataToFit$outputs)+offset.shift
    outpred <- offset + drop(dataToFit$inputs %*% coefs)
  } else {
    outpred <- coefs[1] + dataToFit$inputs * coefs[2]
  }

  delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
  term1 <- delta.t[t.index.min:(t.index.max-1)] / tau * outpred ## vector for the first term
  term2.multiplier <-( 1 - delta.t[t.index.min:(t.index.max-1)]/tau)
  conds <- c("min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  wt.predicted <- numeric(length=(t.index.max-1))
  names(wt.predicted) <- conds[t.index.min:(t.index.max-1)]
  wt.predicted[1] <- term1[1] + term2.multiplier[1] * targ.exp.pulled[1] ## not sure about the last here
  for ( i in (t.index.min+1):(t.index.max-1) ){
    wt.predicted[i] <- term1[i] + term2.multiplier[i] * wt.predicted[i-1]
  }
  wt.predicted.fullscale <- max.intensity[targ] * wt.predicted

  ## boost back forward
  wt.predicted.boosted <- approx(array.times[-1],wt.predicted.fullscale,array.times[-1]-target.pull,rule=2)$y

  wt.predicted.boosted
}

