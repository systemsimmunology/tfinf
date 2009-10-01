##
## June 2007
## Pipeline from pdna objects to model fitted objects
##

##
## Runc calcModels to get OLS solutions
##  

library(lars)
library(lasso2)
library(odesolve)

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
interact.dir <- file.path(Sys.getenv("TFINF"),"interaction_data")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
ddata.dir <- file.path(Sys.getenv("TFINF"),"derived_data")

load(paste(interact.dir,"pdnaModels.30September2009.RData",sep="/"))
load(paste(exp.dir,"scaled.mus.objects.RData",sep="/"))
##load(paste(ddata.dir,"boost.vec.refit.RData",sep="/"))

load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))
load(paste(ddata.dir,"boost.vec.RData",sep="/"))


load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(exp.dir,"all.mus.objects.RData",sep="/"))
load(paste(exp.dir,"scaled.mus.objects.RData",sep="/"))

source("./utilitiesTAC.R")
source("./utilitiesFiniteDiff3D.R")
source("./selectionUtilities.R")

##max.intensity <- apply(lps.mus,1,max) 
##lps.mat.max1 <- lps.mus/max.intensity ## this one will be modified for nucloc
##lps.mat.max1.exp <- lps.mus/max.intensity ## this one will retain expression vars 
##source("./createBoostVec.R")

cat ("Starting calcModels runs \n")

t.index.min <- 1
t.index.max <- 11
tau <- 600./log(2.)

mods.enrp.05 <- calcModels(pdna.enrp.05,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=2,zero.offset=TRUE,tau.estimate=TRUE)

mods.enrp.01 <- calcModels(pdna.enrp.01,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=2,zero.offset=TRUE,tau.estimate=TRUE)
mods.enrs.001 <- calcModels(pdna.enrs.001,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)


## January 2008
## Not sure how this arose, but Irf8 -> Irf8 returns empty list
mods.enrs.001 <- mods.enrs.001[which(unlist(lapply(mods.enrs.001,length))>0)]


mods.hs.et.5 <- calcModels(pdna.hs.et.5,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)

mods.hs.et.1 <- calcModels(pdna.hs.et.1,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)

mods.hs.et.05 <- calcModels(pdna.hs.et.05,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)
mods.hs.et.01 <- calcModels(pdna.hs.et.01,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)
mods.hs.et.001 <- calcModels(pdna.hs.et.001,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)

mods.curated <- calcModels(pdna.curated,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)

## Test for reproduciblity
##res <- calcModels(pdna.hs.et.001[1:2],t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)
##rez <- calcModels(pdna.enrp.01[1:2],t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=2,zero.offset=TRUE,tau.estimate=TRUE)


##
## For two-input cases, treat correlated inputs
##
cat ("Treating correlated inputs, for paired sites \n")

tte.cormat <- cor(t(lps.mat.max1[tte,]))
tte.cormat[lower.tri(tte.cormat)] <- 0
diag(tte.cormat) <- 0
### tried using NAs, without much luck
## highly correlated inputs
corTFPairs <- matrixExtrema(tte.cormat,cutoff=0.8,decreasing=TRUE)$pairs
anticorTFPairs <- matrixExtrema(-tte.cormat,cutoff=0.8,decreasing=TRUE)$pairs

result.01 <- pairToSingleMods(mods.enrp.01,rbind(corTFPairs,anticorTFPairs))
mods.enrp.sansCorTFPairs.01 <- result.01$mods.pairs
pdna.singles.fromCorTFPairs.01 <- result.01$targsAndCands.singles
mods.singles.fromCorTFPairs.01 <- calcModels(pdna.singles.fromCorTFPairs.01,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=1,zero.offset=TRUE,tau.estimate=TRUE)

##
## Model shrinkage
##

cat("Model shrinkage/selection begins\n") 
## ( This one takes a while, so use sparingly ! )
result.ss.01 <- shrinkageSelectModels(mods.enrp.sansCorTFPairs.01, t.index.min, t.index.max, tau=0, boost.vec,lps.mat.max1, n.cands=2, zero.offset=TRUE )
mods.enrp.dubs.01 <- result.ss.01$mods.pairs
mods.enrp.sings.01 <- result.ss.01$mods.singles

mods.preode.strings <- c("mods.curated","mods.hs.et.05","mods.hs.et.01","mods.hs.et.001","mods.enrs.001","mods.singles.fromCorTFPairs.01","mods.enrp.sings.01","mods.enrp.dubs.01")
save(list=mods.preode.strings,file=paste(ddata.dir,"models.preode.30Sep2009.RData",sep="/"))

##
## Non-linear ODE fitting
## 

cat("Beginning non-linear ODE fitting\n")

mods.enrp.dubs.01.ode <- computeODEModsAnalytic(mods.enrp.dubs.01,n.cands=2)

mods.enrp.sings.01.ode <- computeODEModsAnalytic(mods.enrp.sings.01,n.cands=1)
mods.singles.fromCorTFPairs.01.ode <- computeODEModsAnalytic( mods.singles.fromCorTFPairs.01, n.cands=1)  
mods.enrs.001.ode <- computeODEModsAnalytic(mods.enrs.001,n.cands=1)

mods.hs.et.05.ode  <- computeODEModsAnalytic(mods.hs.et.05, n.cands=1)
mods.hs.et.01.ode  <- computeODEModsAnalytic(mods.hs.et.01, n.cands=1)
mods.hs.et.001.ode <- computeODEModsAnalytic(mods.hs.et.001,n.cands=1)

mods.curated.ode <- computeODEModsAnalytic(mods.curated,n.cands=1)

#mods.bind.ode <- computeODEModsAnalytic(mods.bind,n.cands=1)

mods.ode.strings <- c("mods.curated.ode","mods.hs.et.05.ode","mods.hs.et.01.ode","mods.hs.et.001.ode","mods.enrs.001.ode","mods.singles.fromCorTFPairs.01.ode","mods.enrp.sings.01.ode","mods.enrp.dubs.01.ode")
save(list=mods.ode.strings,file=paste(ddata.dir,"models.ode.30Sep2009.RData",sep="/"))

###
### RMSD filter
###

source("./randomTFs.R") ## Overall TFs, for all genes (not gene-specific)

cat("Fitering on RMSD\n")

mods.enrp.dubs.01.ode.rmsf <- filterModsByRMSD( mods.enrp.dubs.01.ode, rmsd.threshold=rands.2regs.0.05, rmsd.type="fullode", scale=FALSE,n.cands=2)
mods.curated.ode.rmsf <- filterModsByRMSD( mods.curated.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.hs.et.05.ode.rmsf <- filterModsByRMSD( mods.hs.et.05.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.hs.et.01.ode.rmsf <- filterModsByRMSD( mods.hs.et.01.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.hs.et.001.ode.rmsf <- filterModsByRMSD( mods.hs.et.001.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.enrs.001.ode.rmsf <- filterModsByRMSD( mods.enrs.001.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.singles.fromCorTFPairs.01.ode.rmsf <- filterModsByRMSD( mods.singles.fromCorTFPairs.01.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)
mods.enrp.sings.01.ode.rmsf <- filterModsByRMSD( mods.enrp.sings.01.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)

##mods.bind.ode.rmsf <- filterModsByRMSD( mods.bind.ode, rmsd.threshold=rands.1reg.0.05, rmsd.type="fullode", scale=FALSE,n.cands=1)

mods.rmsf.strings <- c("mods.curated.ode.rmsf","mods.hs.et.05.ode.rmsf","mods.hs.et.01.ode.rmsf","mods.hs.et.001.ode.rmsf","mods.enrs.001.ode.rmsf","mods.singles.fromCorTFPairs.01.ode.rmsf","mods.enrp.sings.01.ode.rmsf","mods.enrp.dubs.01.ode.rmsf")

save(list=mods.rmsf.strings,file=paste(ddata.dir,"models.rmsf.30Sep2009.RData",sep="/"))

## This may be confusing. Fisrt 01 applies to .. , second to randomization scores
mods.enrp.dubs.01.1.ode.rmsf <- filterModsByRMSD( mods.enrp.dubs.01.ode, rmsd.threshold=rands.2regs.0.10, rmsd.type="fullode", scale=FALSE,n.cands=2)

mods.enrp.dubs.01.2.ode.rmsf <- filterModsByRMSD( mods.enrp.dubs.01.ode, rmsd.threshold=rands.2regs.0.20, rmsd.type="fullode", scale=FALSE,n.cands=2)

mods.enrp.dubs.01.5.ode.rmsf <- filterModsByRMSD( mods.enrp.dubs.01.ode, rmsd.threshold=rands.2regs.0.50, rmsd.type="fullode", scale=FALSE,n.cands=2)
