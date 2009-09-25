
load("~/macrophage/AffyArray/newAffy/annotation.objects.RData")
source("~/macrophage/RShared/tallyUtilities.R")

##commented out April/May 2008. Hope these are no longer needed
##load("variousObjects.RData")
##Load TF network distance
##load("NetworkDistance/TFNetDist.RData") ### this loads tf.dist.cn 
##noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information

##
## Load TRE.OUT
##

load("/proj/ilyalab/Vesteinn/data/TransScan/RcodeFromDan/Parsed.Scan.eset.allMouseSeptember2009.RData")

##load("../../../data/TransScan/RcodeFromDan/Parsed.Scan.eset.allMouseSeptember2009.RData")

##load("../../../data/TransScan/RcodeFromDan/Parsed.Scan.eset.allMouseMarch2008.RData")


TRE.OUT <- TRE.OUT.eset
## Sep 17, 2009: Length 4742
expressed.scanned.egid <- unique(as.character(unlist(lapply(lapply(TRE.OUT,"[[","seq.info"),"[[","entrez.id")))) ## entrez IDs
load("~/macrophage/AffyArray/newAffy/representativeProbes.RData")
expressed.scanned.ps <- sort(as.character(repProbes.ncbiID[expressed.scanned.egid]))


## May 12, 2008: Length 4869=4880-11
## The 11 that did not pass through SliceTREOUT.R
## seem all to have dashes in the gene symbol
## and may hence have missed out on parsing 
##Not found: 52270 
##Not found: 24429 
##Not found: 73421 
##Not found: 73411 
##Not found: 61232 
##Not found: 16206 
##Not found: 55413 
##Not found: 53835 
##Not found: 72806 
##Not found: 71281 
##Not found: 8686 

##
## Load m.tally,m.pair.tally,mf.tally,mdo.pair.tally,mf.pair.tally,mdof.pair.tally
##


load("/proj/ilyalab/Vesteinn/data/TransScan/RcodeFromDan/hitTally.allMouse.September2009.RData")

##load("../../../data/TransScan/RcodeFromDan/hitTally.allMouse.thresholdA.RData")

##load("/users/thorsson/data/TransScan/RcodeFromDan/hitTally.allMouse.q9999.RData")

## make sure tallyUtilities.R have been loaded are loaded 

## Read all mouse matrices and lengths
gg <-read.table("../../RShared/matrixLengths",as.is=TRUE); matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg);
mouseMatrices <- names(matrixLength)
 
## load familyMap, a mapping from vector of matrix names to family strings
load("../../RShared/familyMapVT3.RData")
familyMap <- familyMapVT3; rm(familyMapVT3);
familyNames <- unique(sort(as.character(familyMap)))

## Load file of threshold for each matrix
load("../../RShared/matrixThresholdA.RData")
matrixThreshold <- matrixThresholdA; 

maxsep <- 100; #the maximum allowable distance between hits to be retained
minsep <- 1; #don't retain pairs if they are less than this distance apart
pdist <- 1000; #don't retain hits if their average distance is more than this from the TSS

## see initialize.R for expressed.scanned.ensembl, etc

## Comment in or out as needed
##targExpressedScans <- TRE.OUT[expressed.scanned.ensembl]
##TRE.OUT <- targExpressedScans
##rm(TRE.OUT)

##targScans <- targExpressedScans

# load lps.mus
load("~/macrophage/AffyArray/newAffy/all.mus.objects.RData")
setdiff(ls()[grep("mus",ls())],"lps.mus") ## remove others from memory

# load lps.ratios
load("~/macrophage/AffyArray/newAffy/all.ratios.objects.RData")
setdiff(ls()[grep("ratios",ls())],"lps.ratios") ## remove others from memory

### Memory problem here? three copies?
## 
esp <- expressed.scanned.ps 
 ## 
cordistA <- 1-cor(t(lps.mus[esp,]))

lps.ratios.w.zero <- cbind(rep(0.,length(esp)),lps.ratios[esp,])
colnames(lps.ratios.w.zero)[1] <- "min0"
cordistR <- 1-cor(t(lps.ratios.w.zero)) #
## 
##targsToPredict <- esp  ## this may be excessively many targets to predict

##targsToPredict <- names(TFset)



