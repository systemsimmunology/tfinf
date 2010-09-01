library(odesolve)
load("/Users/thorsson/tfinf/derived_data/modvec.e.RData")
load("/Users/thorsson/tfinf/expression_data/scaled.mus.objects.RData")
load("/Users/thorsson/tfinf/expression_data/all.mus.objects.RData")
load("/Users/thorsson/tfinf/annotations/representativeProbes.RData")
load("/Users/thorsson/tfinf/annotations/annotation.objects.RData")
load("/Users/thorsson/tfinf/derived_data/boost.vec.RData")
source("/Users/thorsson/tfinf/R/utilitiesODEsolve.R")
source("/Users/thorsson/tfinf/R/utilitiesMods.R")
source("/Users/thorsson/tfinf/R/utilitiesAnnotation.R")
delta.t <- c(20.,20.,20.,20.,40.,120.,120.,120.,600.,360.)
array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
cc <- cname.compare

## June 2010 - Cytokines
if (FALSE){
ca.ps <- as.character(repProbes.ncbiID[as.character(read.table("/Users/thorsson/data/GeneOntology/CytokineActivity.tsv")$V1)])
modvec <- sliceModVecByTargs(ca.ps,modvec.e)
predPlotWrapperModVec(modvec,savedir="JUNE2010")
createModPropFileMV(modvec,"ModelAttributes.tsv")
createEdgeFileMV(modvec,"EdgeAttributes.tsv")
createModfileLookupFileMV(modvec,"ProfileImages.txt","/proj/ilyalab/Vesteinn/forHector/ProfilePlots","png")
}
## July 2010 - Full models
modvec <- modvec.e
createModPropFileMV(modvec,"ModelAttributes.tsv")
createEdgeFileMV(modvec,"EdgeAttributes.tsv")
createModfileLookupFileMV(modvec,"ProfileImages.txt","ProfilePlots","png")
predPlotWrapperModVec(modvec,savedir="ProfilePlots")
