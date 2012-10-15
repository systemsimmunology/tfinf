
## Global

array.times <- c(0, 20, 40, 60, 80,120,240,360,480,1080,1440)

t <- array.times

getMuNuVecs <- function( ivars,t.index.min,t.index.max,boost.vec=boost.vec, n.cands=2 ){

  inputs.boosted <- lps.mat.max1[ivars,] ## initialize ( copy dimensions and labels )
  if ( is.vector(inputs.boosted) ){
    inputs.boosted <- approx(array.times,lps.mat.max1[ivars,],array.times-boost.vec[ivars],rule=2)$y
  } else {
    for ( ivar in ivars ){
      inputs.boosted[ivar,] <- approx(array.times,lps.mat.max1[ivar,],array.times-boost.vec[ivar],rule=2)$y
    }
  }
  
  ti <- t.index.min
  tf <- t.index.max
  
  ## This should be generalized to a diferent number of inputs

  if ( length(ivars) == 1 ){
    
    X <- inputs.boosted
    mux <- X[ti:(tf-1)]-t[ti:(tf-1)]*(X[(ti+1):tf]-X[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])
    nux <- (X[(ti+1):tf]-X[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])
    
    mus <- as.matrix(t(mux))
    nus <- as.matrix(t(nux))
    
    return(list(mus=mus,nus=nus))
    ##return(list(mux=mux,nux=nux,muy=muy,nuy=nuy))
  }
    
  if ( length(ivars) == 2 ){
    
    X <- inputs.boosted[1,]
    mux <- X[ti:(tf-1)]-t[ti:(tf-1)]*(X[(ti+1):tf]-X[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])
    nux <- (X[(ti+1):tf]-X[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])
    
    Y <- inputs.boosted[2,]
    muy <- Y[ti:(tf-1)]-t[ti:(tf-1)]*(Y[(ti+1):tf]-Y[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])
    nuy <- (Y[(ti+1):tf]-Y[ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)])

    ib <- inputs.boosted
    ## Can we do both at once?
    ## This looks really confusing, in part because t() is transpose and t[] is time
    muu <- ib[,ti:(tf-1)]-t(t[ti:(tf-1)]*t(ib[,(ti+1):tf]-ib[,ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)]))
    nuu <- t(t(ib[,(ti+1):tf]-ib[,ti:(tf-1)])/(t[(ti+1):tf]-t[ti:(tf-1)]))
    
    return(list(mus=muu,nus=nuu))
    ##return(list(mux=mux,nux=nux,muy=muy,nuy=nuy))
  }
    
}

### params is the vector
### c1 = b1/tau
### c2 = b2/tau
### c3 = 1/tau.

### for the ODE,  the intermediates use in the solution are
### d = c1 ( "delta" )
### e = c2 ( "epsilon" )
### l = -c3  ( "lambda"  )



## Function returns vector value
##G1 <- lps.mat.max1[psoi,1]
## G1 , or psoi, should probably be an input variable

##Gvals <- lps.mat.max1[psoi,]
##G1 <- Gvals[1]


####### GLOBALS BELOW
####### t


## Constructs the Gvec function from the value G1
## The Gvec function is 11-vector valued, defined on the par space

GvecConstructor <- function ( G1, tf.terms ){

  mux <- tf.terms$mus[1,]
  muy <- tf.terms$mus[2,]
  nux <- tf.terms$nus[1,]
  nuy <- tf.terms$nus[2,]
  
  ## Given G1 the values at all other array times
  Gvec <- function ( params ) {
    
    d <- as.numeric(params[1])
    e <- as.numeric(params[2])
    l <- as.numeric(-params[3])
    
    Gout <- numeric(length=10)

    Gout[1] <- G1
  
    mu <- d * mux + e * muy
    nu <- d * nux + e * nuy
    
    ## Vector of constant result from differential eq.
    ## c <- numeric(length=10)

    for ( n in 1:10 ){
      
      expo <- exp(l*(t[n+1]-t[n]))
      termg <- Gout[n]*expo
      termt <- nu[n]/l*(t[n]*expo-t[n+1])
      termk <- (nu[n]+l*mu[n])/l/l * ( expo - 1 )
      Gout[n+1] <- termg + termt + termk
      
      ## older, and equivalent
      ##c[n] <- exp(-l*t[n])*(Gout[n]+(nu[n]+l*mu[n])/l/l + nu[n]*t[n]/l)
      ##Gout[n+1] <- -(nu[n]+l*mu[n])/l/l - nu[n]*t[n+1]/l + c[n]*exp(l*t[n+1])

    }
  
    Gout

  }

  return(Gvec)
  
}


## GVecConstructor for single input

GvecConstructorOneBeta <- function ( G1, tf.terms ){

  mux <- tf.terms$mus[1,]
  nux <- tf.terms$nus[1,]
  
  ## Given G1 the values at all other array times
  Gvec <- function ( params ) {
    
    d <- as.numeric(params[1])
    l <- as.numeric(-params[2])
    
    Gout <- numeric(length=10) ## needs to  be generlized

    Gout[1] <- G1
  
    mu <- d * mux
    nu <- d * nux

    for ( n in 1:10 ){
      
      expo <- exp(l*(t[n+1]-t[n]))
      termg <- Gout[n]*expo
      termt <- nu[n]/l*(t[n]*expo-t[n+1])
      termk <- (nu[n]+l*mu[n])/l/l * ( expo - 1 )
      Gout[n+1] <- termg + termt + termk
      
    }
  
    Gout

  }

  return(Gvec)
  
}

##
## Derivative w.r.t first fitting parameter ( beta1/tau )
##

del1GvecConstructor <- function ( tf.terms ) {

  mux <- tf.terms$mus[1,]
  muy <- tf.terms$mus[2,]
  nux <- tf.terms$nus[1,]
  nuy <- tf.terms$nus[2,]

  del1GVec <- function ( params ) {
  
    d <- as.numeric(params[1])
    e <- as.numeric(params[2])
    l <- as.numeric(-params[3])
    
    del1Gout <- numeric(length=11) ## not sure about this
    del1Gout[1] <- 0 ## not sure about this either 
  
    ## Partial derivatives
    dmudd <- mux
    dnudd <- nux
  
    for ( n in 2:11 ){

      expo <- exp(l*(t[n]-t[n-1]))
      dtermgdd <-  del1Gout[n-1] * expo ## ??
      dtermtdd <- dnudd[n-1]/l*(t[n-1]*expo-t[n])
      dtermkdd <- (dnudd[n-1]+l*dmudd[n-1])/l/l * ( expo - 1 )
      
      del1Gout[n] <- dtermgdd + dtermtdd + dtermkdd 
 
    }
  
    del1Gout
  }

  return( del1GVec )

}

##
## Derivative w.r.t first fitting parameter ( beta1/tau )
##

del1GvecConstructorOneBeta <- function ( tf.terms ) {

  mux <- tf.terms$mus[1,]
  nux <- tf.terms$nus[1,]

  del1GVec <- function ( params ) {
  
    d <- as.numeric(params[1])
    l <- as.numeric(-params[2])
    
    del1Gout <- numeric(length=11) ## not sure about this
    del1Gout[1] <- 0 ## not sure about this either 
  
    ## Partial derivatives
    dmudd <- mux
    dnudd <- nux
  
    for ( n in 2:11 ){

      expo <- exp(l*(t[n]-t[n-1]))
      dtermgdd <-  del1Gout[n-1] * expo ## ??
      dtermtdd <- dnudd[n-1]/l*(t[n-1]*expo-t[n])
      dtermkdd <- (dnudd[n-1]+l*dmudd[n-1])/l/l * ( expo - 1 )
      
      del1Gout[n] <- dtermgdd + dtermtdd + dtermkdd 
 
    }
  
    del1Gout
  }

  return( del1GVec )

}


##
## Derivative w.r.t second fitting parameter ( beta2/tau )
##
del2GvecConstructor <- function ( tf.terms ) {
  
  mux <- tf.terms$mus[1,]
  muy <- tf.terms$mus[2,]
  nux <- tf.terms$nus[1,]
  nuy <- tf.terms$nus[2,]
  
  del2Gvec <- function ( params ) {

    d <- as.numeric(params[1])
    e <- as.numeric(params[2])
    l <- as.numeric(-params[3])
    
    del2Gout <- numeric(length=11) ## not sure about this
    del2Gout[1] <- 0 ## not sure about this either 
  
    ## Partial derivatives
    dmude <- muy
    dnude <- nuy
    
    for ( n in 2:11 ){
      expo <- exp(l*(t[n]-t[n-1]))
      dtermgde <-  del2Gout[n-1] * expo ## ??
      dtermtde <- dnude[n-1]/l*(t[n-1]*expo-t[n])
      dtermkde <- (dnude[n-1]+l*dmude[n-1])/l/l * ( expo - 1 )
      del2Gout[n] <- dtermgde + dtermtde + dtermkde 
      
    }
    
    del2Gout

  }

  return ( del2Gvec )
  
  
}


##
## Derivative w.r.t third fitting parameter ( 1 /tau )
##
## Note that the l paramter differs in sign (only) from 1/tau
## Therefore overall sign flip in last step

del3GvecConstructor <- function ( tf.terms, Gvalues ){

  mux <- tf.terms$mus[1,]
  muy <- tf.terms$mus[2,]
  nux <- tf.terms$nus[1,]
  nuy <- tf.terms$nus[2,]
                                     
  del3Gvec <- function ( params ) {

    d <- as.numeric(params[1])
    e <- as.numeric(params[2])
    l <- as.numeric(-params[3])

    del3Gout <- numeric(length=11) ## not sure about this
    del3Gout[1] <- 0 ## not sure about this either 

    mu <- d * mux + e * muy
    nu <- d * nux + e * nuy
    
    for ( n in 2:11 ){
      expo <- exp(l*(t[n]-t[n-1]))
      dexpodl <- (t[n]-t[n-1])*expo
      
      ## derivative of termg <- Gout[n-1]*expo
      dtermgdl <-  del3Gout[n-1] * expo + Gvalues[n-1] * dexpodl

      ## derivative of termt <- nu[n-1]/l*(t[n-1]*expo-t[n])
      dtermtdl <- -  nu[n-1]*(t[n-1]*expo-t[n])/l/l +  nu[n-1]/l*t[n-1]*dexpodl

      ## derivative of termk <- (nu[n-1]+l*mu[n-1])/l/l * ( expo - 1 )
      dtermkdl <- (-2*nu[n-1]/l/l/l-mu[n-1]/l/l)*( expo - 1 ) + (nu[n-1]+l*mu[n-1])/l/l* dexpodl 

      del3Gout[n] <- dtermgdl + dtermtdl + dtermkdl 
      
    }

    del3Gout <- -del3Gout ## see comment above 
    
    del3Gout

  }
  return ( del3Gvec )
}

## For fitting single input, this fit is del2, not del3 as above
del2GvecConstructorOneBeta <- function ( tf.terms, Gvalues ){

  mux <- tf.terms$mus[1,]
  nux <- tf.terms$nus[1,]
                                     
  del2Gvec <- function ( params ) {

    d <- as.numeric(params[1])
    l <- as.numeric(-params[2])

    del2Gout <- numeric(length=11) ## not sure about this
    del2Gout[1] <- 0 ## not sure about this either 

    mu <- d * mux
    nu <- d * nux
    
    for ( n in 2:11 ){
      expo <- exp(l*(t[n]-t[n-1]))
      dexpodl <- (t[n]-t[n-1])*expo
      
      ## derivative of termg <- Gout[n-1]*expo
      dtermgdl <-  del2Gout[n-1] * expo + Gvalues[n-1] * dexpodl

      ## derivative of termt <- nu[n-1]/l*(t[n-1]*expo-t[n])
      dtermtdl <- -  nu[n-1]*(t[n-1]*expo-t[n])/l/l +  nu[n-1]/l*t[n-1]*dexpodl

      ## derivative of termk <- (nu[n-1]+l*mu[n-1])/l/l * ( expo - 1 )
      dtermkdl <- (-2*nu[n-1]/l/l/l-mu[n-1]/l/l)*( expo - 1 ) + (nu[n-1]+l*mu[n-1])/l/l* dexpodl 

      del2Gout[n] <- dtermgdl + dtermtdl + dtermkdl 
      
    }

    del2Gout <- -del2Gout ## see comment above 
    
    del2Gout

  }
  return ( del2Gvec )
}





####################################################################################


## The analytic gradient of the sse that is minimized

## consructs sse function
sseConstructor <- function ( psoi, GvecFunc ) {
  
  Gvals <- lps.mat.max1.exp[psoi,]
  
  sse <- function ( parms ){
    sse <- sum((GvecFunc(parms)-Gvals)^2)
    return( sse )
  }

  return ( sse )
}


gradSSEConstructor <- function ( psoi, tf.terms, GvecFunc, n.cands=2 ) {
  
  Gvals <- lps.mat.max1.exp[psoi,]
  
  gradSSE <- function ( parms ){
    
    if ( n.cands == 1 ){

      del1Gvec <- del1GvecConstructorOneBeta(tf.terms)
      Gpred <- GvecFunc(parms)
      del2Gvec <- del2GvecConstructorOneBeta(tf.terms,Gpred)  
      
      dG1 <- 2*sum( (Gpred-Gvals)*del1Gvec(parms) )
      dG2 <- 2*sum( (Gpred-Gvals)*del2Gvec(parms) )
      return( c(dG1,dG2) )
      
    } else if ( n.cands == 2 ) {
      
      del1Gvec <- del1GvecConstructor(tf.terms)
      del2Gvec <- del2GvecConstructor(tf.terms)
      Gpred <- GvecFunc(parms)
      del3Gvec <- del3GvecConstructor(tf.terms,Gpred)  
      
      dG1 <- 2*sum( (Gpred-Gvals)*del1Gvec(parms) )
      dG2 <- 2*sum( (Gpred-Gvals)*del2Gvec(parms) )
      dG3 <- 2*sum( (Gpred-Gvals)*del3Gvec(parms) )
    
      return( c(dG1,dG2,dG3) )
    }
  }
  
  return( gradSSE )
}



## Ts are half-lives **** in hours ***********
## nb1, nb2, nT number of steps
sseOnGrid <- function(sseFunc, b1start,b1end,b1step,b2start,b2end,b2step,Tstart,Tend,Tstep){
  
  nb1 <- (b1end-b1start)/b1step + 1 ## The number of b1 points to be evaluated
  nb2 <- (b2end-b2start)/b2step + 1 ## The number of b2 points to be evaluated
  nT <- (Tend-Tstart)/Tstep + 1 ## The number of b3 points to be evaluated
  
  nAll <- nb1*nb2*nT
  cat("About to evaluate grid of ",nAll,"=",nb1,"x",nb2,"x",nT," sse values.\n", sep="")
  
  b1s <- seq(b1start,b1end,b1step)
  b2s <- seq(b2start,b2end,b2step)

  Ts <- seq(Tstart,Tend,Tstep) * 60  ############ Note this 
  
  
  outCube <- rep(NA,nb1*nb2*nT )
  dim(outCube) <- c(nb1,nb2,nT)
  dimnames(outCube)[[3]] <- Ts
  dimnames(outCube)[[2]] <- b2s
  dimnames(outCube)[[1]] <- b1s

  
  for ( T in Ts ){
    tau <- T/log(2.)
    c3 <- 1/tau
    for ( b1 in b1s ){
      c1 <- b1/tau  
      for ( b2 in b2s ){
        c2 <- b2/tau
        pars <- c(c1,c2,c3)
        outCube[as.character(b1),as.character(b2),as.character(T)] <- sse(pars)
      }
    }
  }
  
  outCube
  
}

