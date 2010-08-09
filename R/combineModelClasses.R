
source("./utilitiesMods.R")

##annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
ddata.dir <- file.path(Sys.getenv("TFINF"),"derived_data")

##load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(ddata.dir,"models.rmsf.2Aug2010.RData",sep="/"))

## Contains c("mods.curated.ode.rmsf","mods.hs.et.05.ode.rmsf","mods.hs.et.01.ode.rmsf","mods.hs.et.001.ode.rmsf","mods.enrs.001.ode.rmsf","mods.singles.fromCorTFPairs.01.ode.rmsf","mods.enrp.sings.01.ode.rmsf","mods.enrp.dubs.01.ode.rmsf")

singles.tags <- c("mods.curated","mods.hs.et.05","mods.enrs.001","mods.singles.fromCorTFPairs.01","mods.enrp.sings.01")
doubles.tags <- "mods.enrp.dubs.01"
model.tags <- c(singles.tags,doubles.tags)
modnames <- paste(model.tags,".ode.rmsf",sep="")

modvec.1 <- createModVecFromMultMod ( modnames[1], sto(modnames[1]), n.cands=1)
modvec.2 <- createModVecFromMultMod ( modnames[2], sto(modnames[2]), n.cands=1)
modvec.3 <- createModVecFromMultMod ( modnames[3], sto(modnames[3]), n.cands=1)
modvec.4 <- createModVecFromMultMod ( modnames[4], sto(modnames[4]), n.cands=1)
modvec.5 <- createModVecFromMultMod ( modnames[5], sto(modnames[5]), n.cands=1)
modvec.6 <- createModVecFromMultMod ( modnames[6], sto(modnames[6]), n.cands=2)

modvec.now <- mergeModVecs(modvec.1,modvec.2)
modvec.now <- mergeModVecs(modvec.now,modvec.3)
modvec.now <- mergeModVecs(modvec.now,modvec.4)
modvec.now <- mergeModVecs(modvec.now,modvec.5)
modvec.now <- mergeModVecs(modvec.now,modvec.6)
modvec.e  <- modvec.now
rm(modvec.now)

save(modvec.e,file=paste(ddata.dir,"modvec.e.RData",sep="/"))

osx.path.public <- "/Volumes/isb-2/Vesteinn/MacTransNet/RegulationHypotheses/"
linux.path.public <- "/net/public/Vesteinn/MacTransNet/RegulationHypotheses/"
windows.path.public <- "\\\\Isb-2\\PUBLIC\\Vesteinn\\MacTransNet\\RegulationHypotheses\\"

osx.path.local <- "/Users/thorsson/macrophage/TFinfluence/R/Networks/"
linux.path.local <- "/users/thorsson/macrophage/TFinfluence/R/Networks/"
windows.path.local <- "\\\\Isb-1\\thorsson\\macrophage\\TFinfluence\\R\\Networks\\"

prefix <- "All"

dir.osx.local <- paste(osx.path.local,prefix,"/",sep="")
dir.linux.local <- paste(linux.path.local,prefix,"/",sep="")
dir.win.local <- paste(windows.path.local,prefix,"\\",sep="")

dir.osx.public <- paste(osx.path.public,prefix,"/",sep="")
dir.linux.public <- paste(linux.path.public,prefix,"/",sep="")
dir.win.public <- paste(windows.path.public,prefix,"\\",sep="")

createCytoscapeFilesMV(subnet.cytokines,prefix,outputdir=dir.osx.local,linuxdir=dir.linux.public,windir=dir.win.public,osxdir=dir.osx.public,ppints=FALSE ) 
 
predPlotWrapperModVec(  subnet.cytokines.best ,savedir=dir.osx.local )

createBioTapestryFilesMV(subnet.cytokines,prefix="Cytokines",outputdir=dir.osx.local,regionvec=regionvec)

createBioTapestryFilesMV(subnet.cytokines.best,prefix="Cytokines",outputdir=dir.osx.local,regionvec=regionvec) 

createBioTapestryFilesMV(subnet.cytokines,prefix="Cytokines",outputdir=dir.osx.local,regionvec=regionvec)

