## Randomly assign regulators to genes, to get null model disribution for rmsd score
## Simplest version
## For each number of regs n.reg
## choose random gene. choose random regs. do the fit. compile the distribution

source(file.path(Sys.getenv("UTILS_DIR"),"selectionUtilities.R"))
source(file.path(Sys.getenv("UTILS_DIR"),"utilitiesODEsolve.R"))
source(file.path(Sys.getenv("UTILS_DIR"),"analyticODE.R"))

load(file.path(Sys.getenv("EXP_DIR"),"all.ps.list.objects.RData"))
load(file.path(Sys.getenv("PO_DIR"),"tteMaps.RData"))
load(file.path(Sys.getenv("AUX_DIR"),"boost.vec.RData"))
load(file.path(Sys.getenv("PO_DIR"),"scaled.mus.objects.RData"))
load(file.path(Sys.getenv("SEQ_DIR"),"allMouseMappings.RData"))
load(file.path(Sys.getenv("PO_DIR"),"expressed.scanned.ensembl.RData")) 

tte <- transfac.tfs.expressed

expressed.scanned.egid <- unique(as.character(entrezIDofEnsemblID[expressed.scanned.ensembl]))

target.list <- paste(expressed.scanned.egid,"_at",sep="")
## includes the targets in pdna.curated, according to a check

t.index.min <- 1
t.index.max <- 11
tau <- 600/log(2.)
zero.offset <- TRUE
tau10min <- 10./log(2.)
tau10hr <- (10*60)/log(2.)

if ( FALSE ){
## dev code for using models to get target list
load(paste(ddata.dir,"models.ode.RData",sep="/"))  ## targets used to be derived from pdnaModels.RData. Now use mods.
mods.ode.sings <- c("mods.curated.ode","mods.hs.et.05.ode","mods.hs.et.01.ode","mods.hs.et.001.ode","mods.enrs.001.ode","mods.singles.fromCorTFPairs.01.ode","mods.enrp.sings.01.ode")
targs.sings <- character()
for ( mod in mods.ode.sings ){
  targs.sings <- unique(c(targs.sings,names(sto(mod))))
}
targs.dubs <- names(mods.enrp.dubs.01.ode)
}



library(odesolve)
library(lars)
library(lasso2)

n.trials <- 100

#tf.pool <- setdiff(c(bigDowns,bigUps, midDowns,midUps),rc["E2f5"])
#tf.pool <- setdiff(tte,rc["E2f5"])

tf.pool <- tte

n.targs <- length(target.list)

collection.1reg <- matrix(nrow=n.targs,ncol=n.trials)
rownames(collection.1reg) <- target.list

n.regs <- 1

zero.offset <- TRUE

for ( j in 1:n.targs ){ 
  
  targ <- target.list[j]
  
  for ( i in 1:n.trials ){

    if ( ( i %% 20 ) == 0 ){
      write(c(j,i),file='status.1.reg.txt'); #report status   
    }

    if ( ( i %% 500 ) == 0 ){
      cat("Rand. no.", i ,"\n")
    }
      
    ivars <- sample(tf.pool, n.regs )
    
    ##targ <- sample(targ.pool, 1)
    
    ## We're going to calculate the finite difference approximate solution, just to mimic
    ## initial guess for 'full ODE'
    dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,lps.mat.max1)
    ris <- regInfSelectD2F(dataToFit,s=1,zero.offset)
    coefs <- as.numeric(ris[-1])
    offset <- as.numeric(ris[1])
    
    ## end finite difference fit, beging full ODE fit
    
    G1 <- lps.mat.max1[targ,1]
    
    b1zero <- coefs
    tauzero <- 600/log(2.)
    parzero <- c(b1zero,1.)/tauzero
    
    tf.contribution <- getMuNuVecs( ivars , 1, 11,boost.vec,n.cands=n.regs)
    
    if ( n.regs == 1 ){
      GVec <- GvecConstructorOneBeta(G1, tf.terms=tf.contribution )
    } else if ( n.regs == 2 ){
      GVec <- GvecConstructor(G1, tf.terms=tf.contribution )
    }
    sse <- sseConstructor(targ,GVec)
    gradSSE <- gradSSEConstructor(targ, tf.contribution,GVec,n.cands=n.regs)
    
    control <- list()
    control[["trace"]] <- 0
    control[["reltol"]] <- 1e-8
    fit <- optim(parzero,fn=sse,gr=gradSSE,control=control,method="BFGS")
    
    if ( n.regs == 1 ) {
      if ( (fit$par[2] < 1./tau10hr) | (fit$par[2] > 1./tau10min ) ){
        fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-10000,1/tau10hr),upper=c(10000,1/tau10min))
        taufit  <- 1./fit$par[2]
        bfit <- fit$par[1]*taufit
        names(bfit) <- ivars
      }
    } else if ( n.regs == 2 ){
      if ( (fit$par[3] < 1./tau10hr) | (fit$par[3] > 1./tau10min ) ){
        fit <- optim(parzero,fn=sse,gr=gradSSE,method="L-BFGS-B",lower=c(-10000,-10000,1/tau10hr),upper=c(10000,10000,1/tau10min))
        ## These are not needed, but are stored here. Note indices if used.  
        taufit  <- 1./fit$par[3]
        bfit <- fit$par[1:2]*taufit
        names(bfit) <- ivars
      }
    }
    
    
    fval <- fit$value
    rmsd.ode <- sqrt(fval)/sqrt(t.index.max-t.index.min+1) 
    
##    collecthion[n.regs,i] <- rmsd.ode
    collection.1reg[targ,i] <- rmsd.ode
    
    if ( FALSE ){
      Gproduction <- generateGproduction(generateTFfunctions(ivars))    
      scale.fac <- max.intensity[targ]
      wt.pred.array.times <- scale.fac * lsoda(G1,array.times,Gproduction,fit$par)[,2]
      ndiv=50.
      times.dense <- steppes(ndiv=ndiv)
      xvals.dense <- seq(0,10,1/ndiv)
      wt.pred.dense <- scale.fac * lsoda(G1,times.dense,Gproduction,fit$par)[,2]
      pred.out <- wt.pred.array.times ## needs units
      pred.out.dense <- wt.pred.dense ## needs units ?
      x11()
      plotPrediction.ode(targ,bfit,taufit,wt.pred.array.times, xvals.dense,wt.pred.dense)
    }
    
  } ## end loop over random models for fixed n.regs
}

save(collection.1reg,file=file.path(Sys.getenv("PO_DIR"),"collection.1reg.RData"))


