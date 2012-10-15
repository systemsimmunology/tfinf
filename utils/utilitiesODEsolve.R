
## Requires
## library(odesolve)

## Can precompute this for all TFs , but at the moment we are not smart enough to do this !!!
## The time range indicates which times are used to guide the interpolation - the resulting
## functions 'work' over all times (defined over all times )
generateTFfunctions <- function ( ivars, n.cands=2,t.index.min=1,t.index.max=11,mus.mat.maxnormed=lps.mat.max1 ) {

  inputs.boosted <- mus.mat.maxnormed[ivars,t.index.min:t.index.max] ## initialize ( copy dimensions and labels )
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(array.times[t.index.min:t.index.max],mus.mat.maxnormed[ivars,t.index.min:t.index.max],array.times[t.index.min:t.index.max]-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(array.times[t.index.min:t.index.max],mus.mat.maxnormed[ivar,t.index.min:t.index.max],array.times[t.index.min:t.index.max]-boost.vec[ivar],rule=2)$y
    }
  }
  
  ## Piecewise linear *functions* for tfs

  if ( n.cands == 1) {
    tf <- approxfun(array.times[t.index.min:t.index.max],inputs.boosted,rule=2)
    return(list(tf=tf))
  } else if ( n.cands  == 2){
    tf1 <- approxfun(array.times[t.index.min:t.index.max],inputs.boosted[1,],rule=2)
    tf2 <- approxfun(array.times[t.index.min:t.index.max],inputs.boosted[2,],rule=2)
    return(list(tf1=tf1,tf2=tf2))
  }
 
}

## Gdot
generateGproduction <- function ( tfFuncs, n.cands=2 ){

  if ( n.cands == 1){
    
    tf <- tfFuncs$tf
    Gproduction <- function ( t, G, parms ){
      c1 <- parms[1]
      c2 <- parms[2]
      Gdot <- -c2*G + c1*tf(t)
      return(list(Gdot))
    }
    
  } else if ( n.cands == 2 ){
    
    tf1 <- tfFuncs$tf1
    tf2 <- tfFuncs$tf2
    Gproduction <- function ( t, G, parms ){
      c1 <- parms[1]
      c2 <- parms[2]
      c3 <- parms[3]
      Gdot <- -c3*G + c1*tf1(t)+c2*tf2(t)
      return(list(Gdot))
    }
    
  } ## end n.cands == 2

  return (Gproduction)
    
}

constructFitFunction <- function(psoi){

  Gvals <- lps.mat.max1.exp[psoi,]

  ## Gdot
  Gproduction <- function ( t, G, parms ){
    c1 <- parms[1]
    c2 <- parms[2]
    c3 <- parms[3]
    Gdot <- -c3*G + c1*tf1(t)+c2*tf2(t)
    return(list(Gdot))
  }

  ## del(Gdot)/del(parms), meant to speed up ODE solution
  gradGproduction <- function ( t, G, parms ){
    c1 <- parms[1]
    c2 <- parms[2]
    c3 <- parms[3]
    Gdot <- -c3*G + c1*tf1(t)+c2*tf2(t)
    return(c(tf1(t),tf2(t),-G))
  }

  G0 <- Gvals[1]
  
  Gfit <- function ( parms ){
    out <- lsoda(G0,array.times,Gproduction,jacfunc=gradGproduction,parms,hmin=0.001,rtol=1e-4)
    sse <- sum((out[,2]-Gvals)^2)
    return(sse)
  }

  return(Gfit)

}

## Gdot
Gproduction <- function ( t, G, parms ){
  c1 <- parms[1]
  c2 <- parms[2]
  c3 <- parms[3]
  Gdot <- -c3*G + c1*tf1(t)+c2*tf2(t)
  return(list(Gdot))
}

gradGproduction <- function ( t, G, parms ){
  c1 <- parms[1]
  c2 <- parms[2]
  c3 <- parms[3]
  Gdot <- -c3*G + c1*tf1(t)+c2*tf2(t)
  return(c(tf1(t),tf2(t),-G))
}


## The sse that is minimized
Gfit <- function ( parms ){
  out <- lsoda(G0,array.times,Gproduction,jacfunc=gradGproduction,parms,hmin=0.001,rtol=1e-4)
  sse <- sum((out[,2]-Gvals)^2)
  return(sse)
}

## The sse that is minimized
Gfit.eode <- function ( parms ){
  Gpred <- Gvec( parms )
  sse <- sum((Gpred-Gvals)^2)
  return(sse)
}
 
## compute
steppes <- function(ndiv=50,t.index.min=1,t.index.max=11){
  at <- array.times
  dt <- delta.t
  ov <- 0.
  for ( i in 1:10 ){
    ov <- c(ov,seq(at[i],at[i+1],dt[i]/ndiv)[-1])
  }
  return(ov)
}  

plotPrediction.ode <- function(target.pid,coefs,taufit,pred.out,xvals.dense,pred.out.dense,t.index.min=1,t.index.max=11) {

  xvals <- 0:10 ## reference scale for the plot

  ivars <- names(coefs)
  coefs <- as.numeric(coefs)
  
  profile.true <- lps.mus[target.pid,1:t.index.max] ## we'll plot the initial values, even if they were not use for the fit

  profile.ivars <- 0.5 * max.intensity[target.pid] * lps.mat.max1[ivars,1:t.index.max]
  
  ymax <- max(c(profile.true,pred.out))

  thalf <- round(taufit*log(2))
  main=paste(gid2gname(target.pid),", Probeset: ",target.pid,", T-half=",thalf," min",sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',main=main,lwd=3,axes=FALSE)

  axislabels <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  axis(1,xvals[1:t.index.max],labels=axislabels[1:t.index.max])
  axis(2)

  points(xvals[t.index.min:t.index.max],profile.true[t.index.min:t.index.max],col='Blue', pch=19 )

  ## Predicted profile 
  points(xvals[t.index.min:t.index.max],pred.out, col='Blue', pch=19 )
  ## 2nd slot may need modification
  lines(xvals.dense,pred.out.dense, col='Blue', lty=2 )

  ## TFS

  colorvec = c('red','green','pink','cyan','magenta')
  ##colorvec = c('cyan','magenta','pink','red','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(ivars))
  lty <- rep(lty, length.out=length(ivars))

  if ( length(ivars) > 1 ){
    for ( i in 1:length(ivars) ){
      lines( xvals, profile.ivars[i,1:t.index.max], col=colorvec[i],lty=1,lwd=3)
    }
  }
  if ( length(ivars) == 1 ){
    lines( xvals, profile.ivars[1:t.index.max], col=colorvec[1],lty=1,lwd=3)
  }
  
  s1 <- gid2gname(ivars)
  s2 <- as.character(round(coefs,digits=1))
  s2 <- paste("(",s2,")",sep="")
  legend <- paste(s1,s2)

  lg <- c(1.5,0.9*max.intensity[target.pid] )
  legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=3)

  par(op)
}

plotPrediction.ode.hack <- function(target.pid,coefs,taufit,pred.out,xvals.dense,pred.out.dense,t.index.min=1,t.index.max=11){

  xvals <- 0:10 ## reference scale for the plot

  ivars <- names(coefs)
  coefs <- as.numeric(coefs)
  
  profile.true <- lps.mus[target.pid,1:t.index.max] ## we'll plot the initial values, even if they were not use for the fit

  profile.ivars <- 0.5 * max.intensity[target.pid] * lps.mat.max1[ivars,1:t.index.max]
  
  ymax <- max(c(profile.true,pred.out))

  thalf <- round(taufit*log(2))
  main=paste(gid2gname(target.pid),":",target.pid,":T=",thalf,sep="")

  ##op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  op <- par(cex=0.525)
  
  plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',main=main,lwd=1,lty=1,axes=FALSE)

  axislabels <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  axis(1,xvals[1:t.index.max],labels=axislabels[1:t.index.max])
  axis(2)

  points(xvals[t.index.min:t.index.max],profile.true[t.index.min:t.index.max],col='Blue', pch=17 )

  ## Predicted profile 
  points(xvals[t.index.min:t.index.max],pred.out, col='Blue', pch=19 )
  ## 2nd slot may need modification
  lines(xvals.dense,pred.out.dense, col='Blue', lty=2 )
  
  ## TFS

  colorvec = c('brown','magenta','red','pink','green','black')
  ##colorvec = c('red','black','pink','cyan','green','magenta')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(ivars))
  lty <- rep(lty, length.out=length(ivars))

  if ( length(ivars) > 1 ){
    for ( i in 1:length(ivars) ){
      lines( xvals, profile.ivars[i,1:t.index.max], col=colorvec[i],lty=1,lwd=1)
    }
  }
  if ( length(ivars) == 1 ){
    lines( xvals, profile.ivars[1:t.index.max], col=colorvec[1],lty=1,lwd=1)
  }
  
  s1 <- gid2gname(ivars)
  s2 <- as.character(round(coefs,digits=1))
  s2 <- paste("(",s2,")",sep="")
  legend <- paste(s1,s2)

  lg <- c(1.5,0.9*max.intensity[target.pid] )
  legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=1)

  par(op)
}


plotPrediction.ode.minimal <- function(target.pid,coefs,taufit,pred.out,xvals.dense,pred.out.dense,t.index.min=1,t.index.max=11){

  xvals <- 0:10 ## reference scale for the plot

  ivars <- names(coefs)
  coefs <- as.numeric(coefs)
  
  profile.true <- lps.mus[target.pid,1:t.index.max] ## we'll plot the initial values, even if they were not use for the fit

  profile.ivars <- 0.5 * max.intensity[target.pid] * lps.mat.max1[ivars,1:t.index.max]
  
  ymax <- max(c(profile.true,pred.out))

  thalf <- round(taufit*log(2))
  #main=paste(gid2gname(target.pid))
  main=""
  #main=paste(gid2gname(target.pid),":",target.pid,":T=",thalf,sep="")

  op <- par(font=2,lwd=5,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab="",ylab="",main=main,lwd=5,axes=FALSE)

  axislabels <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  ##axis(1,xvals[1:t.index.max],labels=rep("",t.index.max))
  ##axis(2)

  points(xvals[t.index.min:t.index.max],profile.true[t.index.min:t.index.max],col='Blue', pch=19 )

  ## Predicted profile 
  points(xvals[t.index.min:t.index.max],pred.out, col='Blue', pch=19 )
  ## 2nd slot may need modification
  lines(xvals.dense,pred.out.dense, col='Blue', lty=2 )

  ## TFS

  colorvec = c('red','green','pink','cyan','magenta')
  ##colorvec = c('cyan','magenta','pink','red','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(ivars))
  lty <- rep(lty, length.out=length(ivars))

  if ( length(ivars) > 1 ){
    for ( i in 1:length(ivars) ){
      lines( xvals, profile.ivars[i,1:t.index.max], col=colorvec[i],lty=1,lwd=5)
    }
  }
  if ( length(ivars) == 1 ){
    lines( xvals, profile.ivars[1:t.index.max], col=colorvec[1],lty=1,lwd=5)
  }
  
  # s1 <- gid2gname(ivars)
  #s2 <- as.character(round(coefs,digits=1))
  #s2 <- paste("(",s2,")",sep="")
  #legend <- paste(s1,s2)

  #lg <- c(1.5,0.9*max.intensity[target.pid] )
  #legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=3)

  par(op)
}



##
## Currently only works with specific TFs
## 
computeODEMods <- function( mods, n.cands=2){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]

    ## TEMPORARY: JUST ONE SET OF BETAS
    b1zero <- mod$betas[1,1]
    b2zero <- mod$betas[1,2]
    tauzero <- 800/log(2.)

    Gfit <- constructFitFunction(psoi)

    ## Do the fit 
    ##parzero <- c(b1zero,b2zero,1.)/tauzero
    ## Xpo1 is senstive to starting point ( -2 <-> 1.9 )
    parzero <-  c(-2,4,1)/1000
    
    control <- list()
    control[["trace"]] <- 2
    control[["reltol"]] <- 1e-6
    fit <- optim(parzero,Gfit,control=control)

    taufit  <- 1./fit$par[3]
    bfit <- fit$par[1:2]*taufit
    fval <- fit$value
    rmsd.ode <- sqrt(fval)/sqrt(11) ## make this more general

    returnList[[psoi]] <- NULL

    returnList[[psoi]]$tfs <- mod$tfs
    returnList[[psoi]]$tfs.cname <- mod$tfs.cname
    returnList[[psoi]]$betas <- mod$betas   
    returnList[[psoi]]$offsets <- mod$offsets
    returnList[[psoi]]$rmsds <- mod$rmsds
    returnList[[psoi]]$cors <- mod$cors
    returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled
    returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled
    returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter
    returnList[[psoi]]$gcv <- mod$gcv

    returnList[[psoi]]$betas.ode <- as.numeric(bfit)
    returnList[[psoi]]$tau.ode <- taufit
    returnList[[psoi]]$rmsd.ode <- rmsd.ode
  }

  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}


##
## 
## 
computeODEModsAnalytic <- function( mods, n.cands=2){

  returnList <- list()
  psois <- names(mods)

  ## Can be used for constraints
  tau24 <- (24*60)/log(2.)
  tau10hr <- (10*60)/log(2.)
  tau10min <- 10./log(2.)
  
  counter <- 0 
  for ( psoi in psois ){

    ##cat("psoi:",psoi,"\n")
    counter <- counter + 1
    write(c(psoi,counter,counter/length(psois)),file='status.txt'); #report status

    mod <- mods[[psoi]]

    returnList[[psoi]] <- NULL
        
    ## Copy original list contents
    returnList[[psoi]]$tfs <- mod$tfs
    returnList[[psoi]]$tfs.cname <- mod$tfs.cname
    returnList[[psoi]]$betas <- mod$betas   
    returnList[[psoi]]$offsets <- mod$offsets
    returnList[[psoi]]$rmsds <- mod$rmsds
    returnList[[psoi]]$cors <- mod$cors
    returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled
    returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled
    returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter
    returnList[[psoi]]$gcv <- mod$gcv

    if ( !is.null( returnList[[psoi]]$binding.partner )){
      returnList[[psoi]]$binding.partner <- mod$binding.partner
    }
    
    ## Compute ODE-based solutions        
    G1 <- lps.mat.max1.exp[psoi,1]

    ## Single input
    if ( n.cands == 1 ) {
      n.sings <- length(mod$tfs)    

      for ( n in 1:n.sings ){

        tfs <- mod$tfs[n]
##        cat(cc[tfs],"\n")
        
        b1zero <- mod$betas[n]
        tauzero <- 800/log(2.)
        parzero <- c(b1zero,1.)/tauzero
        
        tf.contribution <- getMuNuVecs( tfs , 1, 11,boost.vec,n.cands=1)
        
        GVec <- GvecConstructorOneBeta(G1, tf.terms=tf.contribution )
        sse <- sseConstructor(psoi,GVec)
        gradSSE <- gradSSEConstructor(psoi, tf.contribution,GVec,n.cands=1)
        
        control <- list()
        control[["trace"]] <- 0
        control[["reltol"]] <- 1e-8
        fit <- optim(parzero,fn=sse,gr=gradSSE,control=control,method="BFGS")
        

        if ( (fit$par[2] < 1./tau10hr) | (fit$par[2] > 1./tau10min ) ){
          fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-10000,1/tau10hr),upper=c(10000,1/tau10min))
##          fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-Inf,1/tau10hr),upper=c(Inf,1/tau10min))
        }
        
      
        taufit  <- 1./fit$par[2]
        bfit <- fit$par[1]*taufit
        fval <- fit$value
        rmsd.ode <- sqrt(fval)/sqrt(t.index.max-t.index.min+1) 

        
        returnList[[psoi]]$betas.ode <- c(returnList[[psoi]]$betas.ode,as.numeric(bfit))
        returnList[[psoi]]$tau.ode <- c(returnList[[psoi]]$tau.ode,taufit)
        returnList[[psoi]]$rmsd.ode <- c(returnList[[psoi]]$rmsd.ode,rmsd.ode)
        returnList[[psoi]]$convergence.ode <- c(returnList[[psoi]]$convergence.ode,fit$convergence)
        if ( is.null(fit$message) ){
          fit.message <- "" ## required for vector order below
        } else {
          fit.message <- fit$message
        }
        returnList[[psoi]]$message.ode <- c(returnList[[psoi]]$message.ode,fit.message)

      }
    } ## end tasks for one input
    
    if ( n.cands == 2 ) {
      n.pairs <- nrow(mod$tfs)    
      for ( n in 1:n.pairs ){

        tfs <- mod$tfs[n,]
        
        b1zero <- mod$betas[n,1]
        b2zero <- mod$betas[n,2]
        tauzero <- 800/log(2.)
        parzero <- c(b1zero,b2zero,1.)/tauzero
        
        tf.contribution <- getMuNuVecs( tfs , 1, 11,boost.vec)
        
        GVec <- GvecConstructor(G1, tf.terms=tf.contribution )
        sse <- sseConstructor(psoi,GVec)
        gradSSE <- gradSSEConstructor(psoi, tf.contribution,GVec)
        
        control <- list()
        control[["trace"]] <- 0
        control[["reltol"]] <- 1e-8
        fit <- optim(parzero,fn=sse,gr=gradSSE,control=control,method="BFGS")
        

        if ( (fit$par[3] < 1./tau10hr) | (fit$par[3] > 1./tau10min ) ){
          fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-Inf,-Inf,1/tau10hr),upper=c(Inf,Inf,1/tau10min))
        }
        
      
        taufit  <- 1./fit$par[3]
        bfit <- fit$par[1:2]*taufit
        fval <- fit$value
        rmsd.ode <- sqrt(fval)/sqrt(t.index.max-t.index.min+1) ## make this more general
        
        returnList[[psoi]]$betas.ode <- rbind(returnList[[psoi]]$betas.ode,as.numeric(bfit))
        returnList[[psoi]]$tau.ode <- c(returnList[[psoi]]$tau.ode,taufit)
        returnList[[psoi]]$rmsd.ode <- c(returnList[[psoi]]$rmsd.ode,rmsd.ode)
        returnList[[psoi]]$convergence.ode <- c(returnList[[psoi]]$convergence.ode,fit$convergence)
        if ( is.null(fit$message) ){
          fit.message <- "" ## required for vector order below
        } else {
          fit.message <- fit$message
        }
        returnList[[psoi]]$message.ode <- c(returnList[[psoi]]$message.ode,fit.message)
      }
    } ## end tasks for two inputs
  } ## end loop over psois
      
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}
#######################################################
## 
## Evaluate optima on a grid of delay values  
## 
#######################################################

## mod is a 'starting model' containing first guesses for beta

delayGridOptima <- function ( tf1, tf2 , psoi, mod, verbose=FALSE ){

  nucloc.delays <- c(0,30)
  expression.delays <- c(0,30,60,120,240)

  tf1.delays <- expression.delays
  tf2.delays <- expression.delays
  if ( tf1 %in% rc[nucloc.input.set] ){ tf1.delays <- nucloc.delays }
  if ( tf2 %in% rc[nucloc.input.set] ){ tf2.delays <- nucloc.delays }
  tf1.delays.string <- paste(tf1.delays,"min",sep="")
  tf2.delays.string <- paste(tf2.delays,"min",sep="")

  sseMat <- matrix(NA,nrow=length(tf1.delays),ncol=length(tf2.delays))
  rownames(sseMat) <- tf1.delays.string
  colnames(sseMat) <- tf2.delays.string
  
  row.ind <- as.numeric(which( (mod$tfs[,1] == tf1) & (mod$tfs[,2] == tf2) ))
  b1zero <- mod$betas[row.ind,1]
  b2zero <- mod$betas[row.ind,2]
  tauzero <- 800/log(2.)
  
  parzero <- c(b1zero,b2zero,1.)/tauzero
  
  control <- list()
  control[["trace"]] <- 0
  control[["reltol"]] <- 1e-8

  ivars <- mod$tfs[row.ind,]

  #Save for posterity
  boost.vec.orig.tf1 <- boost.vec[tf1]
  boost.vec.orig.tf2 <- boost.vec[tf2]

  tau24 <- (24*60)/log(2.)
  tau10hr <- (10*60)/log(2.)
  tau10min <- 10./log(2.)
       
  for ( i in 1:length(tf1.delays) ){
    for ( j in 1:length(tf2.delays) ){

      if ( verbose ){
        cat(tf1.delays.string[i],"\n")
        cat(tf2.delays.string[j],"\n")
      }
      
      boost.vec[tf1] <- tf1.delays[i]
      boost.vec[tf2] <- tf2.delays[j]
      
      G1 <- lps.mat.max1.exp[psoi,1]
      tf.contribution <- getMuNuVecs( ivars, 1, 11, boost.vec)
      GVec <- GvecConstructor(G1, tf.contribution )
      
      sse <- sseConstructor(psoi,GVec)
      gradSSE <- gradSSEConstructor(psoi, tf.contribution,GVec)
    
      parzero <- c(b1zero,b2zero,1.)/tauzero

      control <- list()
      if ( verbose ){
        control[["trace"]] <- 2
      } else {
        control[["trace"]] <- 0
      }
      control[["reltol"]] <- 1e-8
      control[["maxit"]] <- 100
      fit <- optim(parzero,fn=sse,gr=gradSSE,control=control,method="BFGS")
    
      if ( (fit$par[3] < 1./tau10hr) | (fit$par[3] > 1./tau10min ) ){
        fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-Inf,-Inf,1/tau10hr),upper=c(Inf,Inf,1/tau10min))
      }
      
      taufit  <- 1./fit$par[3]
      bfit <- fit$par[1:2]*taufit
      names(bfit) <- ivars
      fval <- fit$value
      rmsd.ode <- sqrt(fval)/sqrt(11) ## make this more general
      
      ##cat("tau:",taufit, "b1:",bfit[1], "b2:",bfit[2],"fval:",fval,"\n")
    
      sseMat[i,j] <- rmsd.ode
    } ## end j loop
  }## end i loop

  ## For now be conservative and keep the original values
  boost.vec[tf1] <- boost.vec.orig.tf1
  boost.vec[tf2] <- boost.vec.orig.tf2

  return( sseMat )
 
}

###################################################################
# See above.
## Cube of values, for multiple genes
###################################################################


delayGridCube <- function ( tf1, tf2, modset, verbose=FALSE ){
  
  psois <- names(modset)
  
  nucloc.delays <- c(0,30)
  expression.delays <- c(0,30,60,120,240)

  tf1.delays <- expression.delays
  tf2.delays <- expression.delays
  if ( tf1 %in% rc[nucloc.input.set] ){ tf1.delays <- nucloc.delays }
  if ( tf2 %in% rc[nucloc.input.set] ){ tf2.delays <- nucloc.delays }
  tf1.delays.string <- paste(tf1.delays,"min",sep="")
  tf2.delays.string <- paste(tf2.delays,"min",sep="")
  
  nd1 <- length(tf1.delays.string)
  nd2 <- length(tf2.delays.string)
  
  sseOptCube <- rep(NA,length(psois)*nd1*nd2)
  
  dim(sseOptCube) <- c(length(psois),nd1,nd2)
  
  dimnames(sseOptCube)[[3]] <- tf2.delays.string
  dimnames(sseOptCube)[[2]] <- tf1.delays.string
  dimnames(sseOptCube)[[1]] <- psois
  
  for ( psoi in psois ){
    ##cat(psoi,"\n")
    mod <- modset[[psoi]]
    sseOptCube[psoi,,] <- delayGridOptima(tf1,tf2,psoi,mod,verbose=FALSE)  
  }
  
  sseOptCube
  
}

#######################################################
## Plot 'contact sheets' with predictions
## filename must include full path
gridPredPlotWrapperODE <- function(mods,file="gridODE.ps",nx=1,ny=1, n.cands=2) {

  ### KEEP AN EYE ON THIS VARIABLE ##########
  zero.offset <- TRUE
  
  psois <- names(mods)
  
  ## could do an input parameter check
  
  ## psois follow order presented to this function 
  
  ## PNG - seem to work only one page at a time, and gives crummy resolution
  ## png(file,height=800,width=1100,quality=100)
  
  ## EPS - may work! EPS not rendered in Word ( mac )
  ## postscript(file,onefile = FALSE,paper='letter')

  ## PS
  postscript(file,paper='letter')

  ##Can't get pdf to do 'landscape'
  ##pdf(file,paper='letter',height=8.5,width=11)
 
  par(mfrow=c(nx,ny),mar=(c(3, 2, 2, 1)+0.1),mgp=c(1.75, 0.75, 0))
  par(mfrow=c(nx,ny),mar=(c(2, 2, 1, 1)+0.1),mgp=c(1.75, 0.75, 0))
  
  counter <- 0
  for ( psoi in psois ) {
    ##    for ( psoi in psois[1:min(length(psois),nx*ny)] ) {
    ##cat(psoi,"\n")
    
    if ( n.cands == 1 ){

      
      n.sings <- length(mods[[psoi]]$betas)
      
      for ( i in 1:n.sings ){
        
        ivars <- drop(mods[[psoi]]$tfs[i])
        
        betas <- drop(mods[[psoi]]$betas[i])
        offset <- drop(mods[[psoi]]$offsets[i])
        
        names(betas) <- ivars
        
        betas.ode <- drop(mods[[psoi]]$betas.ode[i])
        names(betas.ode) <- ivars
        tau.ode <- drop(mods[[psoi]]$tau.ode[i])
        
        parsfit <- c(betas.ode,1)/tau.ode
        
        scale.fac <- max.intensity[psoi]
        G0 <- lps.mat.max1.exp[psoi,t.index.min]

        ##cat("ivars",ivars,"\n")
        Gproduction <- generateGproduction(generateTFfunctions(ivars,n.cands=1),n.cands=1)    
        
        wt.pred.array.times <- scale.fac*lsoda(G0,array.times,Gproduction,parsfit)[,2]
        ndiv=50.
        times.dense <- steppes(ndiv=ndiv)
        xvals.dense <- seq(0,10,1/ndiv)

        wt.pred.dense <- scale.fac*lsoda(G0,times.dense,Gproduction,parsfit)[,2]
        
        pred.out <- max.intensity[psoi]*wt.pred.array.times ## needs units
        pred.out.dense <- max.intensity[psoi]*wt.pred.dense ## needs units ?
        
        plotPrediction.ode.hack(psoi,betas.ode,tau.ode,wt.pred.array.times, xvals.dense,wt.pred.dense)

      }
    } ## end n.cands = 1
    
    
    
    
    if ( n.cands == 2 ){
    
      n.pairs <- nrow(mods[[psoi]]$betas)
      
      for ( i in 1:n.pairs ){
        
        ivars <- drop(mods[[psoi]]$tfs[i,])
        
        ##cat(i, ivars,"\n")
        
        betas <- drop(mods[[psoi]]$betas[i,])
        offset <- drop(mods[[psoi]]$offsets[i])
        
        names(betas) <- ivars
        
        betas.ode <- drop(mods[[psoi]]$betas.ode[i,])
        names(betas.ode) <- ivars
        tau.ode <- drop(mods[[psoi]]$tau.ode[i])
        
        parsfit <- c(betas.ode,1)/tau.ode
        
        scale.fac <- max.intensity[psoi]
        G0 <- lps.mat.max1.exp[psoi,t.index.min]
        
        Gproduction <- generateGproduction(generateTFfunctions(ivars))    
        
        wt.pred.array.times <- scale.fac*lsoda(G0,array.times,Gproduction,parsfit)[,2]
        ndiv=50.
        times.dense <- steppes(ndiv=ndiv)
        xvals.dense <- seq(0,10,1/ndiv)
        wt.pred.dense <- scale.fac*lsoda(G0,times.dense,Gproduction,parsfit)[,2]
        
        pred.out <- max.intensity[psoi]*wt.pred.array.times ## needs units
        pred.out.dense <- max.intensity[psoi]*wt.pred.dense ## needs units ?
        
        plotPrediction.ode.hack(psoi,betas.ode,tau.ode,wt.pred.array.times, xvals.dense,wt.pred.dense)
      }
    } ## end n.cands = 2

    
  } ## end psois loop
  
  dev.off()
  
  
}

## August 2007
## mod can be a single element of a modvec, or an older style mod object. Supply targ for that reason.
## Assumed that there is an ODE fit

## Note the two different inputs for the profiles
## mus.mat.maxnormed.exp is expression only
## mus.mat.maxnormed can include adjustments for nuclear localization

calculateTimeEvolPrediction.ode <- function(targ,mod,t.index.min,t.index.max, mus.mat.maxnormed.exp,mus.mat.maxnormed ){

  ivars <- mod$tfs
  betas.ode <- mod$betas.ode
  tau.ode <- mod$tau.ode
  parsfit <- c(betas.ode,1)/tau.ode       
  scale.fac <- max.intensity[targ] ## This is a single scale.fac for this gene ( but is obtained from LPS)
  G0 <- mus.mat.maxnormed.exp[targ,t.index.min] ## this used to be lps.mat.max1.exp : keep that in mind

  n.tfs <- length(mod$tfs)
  
  if ( n.tfs==1 ){
    Gproduction <- generateGproduction(generateTFfunctions(ivars,n.cands=1,t.index.min=t.index.min,t.index.max=t.index.max,mus.mat.maxnormed=mus.mat.maxnormed),n.cands=1)    
  }  else if ( n.tfs==2 ) {        
    Gproduction <- generateGproduction(generateTFfunctions(ivars,n.cands=2,t.index.min=t.index.min,t.index.max=t.index.max,mus.mat.maxnormed=mus.mat.maxnormed))    
  } else {
    stop("Oops, not set up for other than 1 or 2 regulator ODE solutions")
  }
  
  wt.pred.array.times <- scale.fac*lsoda(G0,array.times,Gproduction,parsfit)[t.index.min:t.index.max,2]
  ##pred.out <- max.intensity[targ]*wt.pred.array.times ## needs units
  
  ##Dense, in time predictions, leave out for now
  ##times.dense <- steppes(ndiv=ndiv)
  ##xvals.dense <- seq(0,10,1/ndiv)
  ##wt.pred.dense <- scale.fac*lsoda(G0,times.dense,Gproduction,parsfit)[,2]
  ## ndiv=50.
  ##pred.out.dense <- max.intensity[targ]*wt.pred.dense ## needs units ?

  return(wt.pred.array.times)
  
}
