
annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
  
load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(annot.dir,"allMouseMappings.August2010.RData",sep="/"))

source("./tallyUtilities.R")
source("./utilitiesHitMat.R")


##commented out April/May 2008. Hope these are no longer needed
##load("variousObjects.RData")
##Load TF network distance
##load("NetworkDistance/TFNetDist.RData") ### this loads tf.dist.cn 
##noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information

##
## Load TRE.OUT
##

load(paste(seq.dir,"Parsed.Scan.eset.allMouseAug2010.RData",sep="/"))

TRE.OUT <- TRE.OUT.eset
expressed.scanned.egid <- unique(as.character(unlist(lapply(lapply(TRE.OUT,"[[","seq.info"),"[[","entrez.id")))) ## entrez IDs

load(paste(annot.dir,"representativeProbes.RData",sep="/"))
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

load(paste(seq.dir,"hitTally.allMouse.Aug2010.RData" ,sep="/"))

## Read all mouse matrices and lengths
gg <-read.table(paste(annot.dir,"matrixLengths" ,sep="/"),as.is=TRUE);
matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg);
mouseMatrices <- names(matrixLength)
 
## load familyMap, a mapping from vector of matrix names to family strings
load(paste(annot.dir,"familyMapVT3.RData" ,sep="/"))
familyMap <- familyMapVT3; rm(familyMapVT3);
familyNames <- unique(sort(as.character(familyMap)))

## Load file of threshold for each matrix
load(paste(annot.dir,"matrixThresholdA.RData" ,sep="/"))
matrixThreshold <- matrixThresholdA; 

maxsep <- 100; #the maximum allowable distance between hits to be retained
minsep <- 1; #don't retain pairs if they are less than this distance apart
pdist <- 1000; #don't retain hits if their average distance is more than this from the TSS

## Comment in or out as needed
##targExpressedScans <- TRE.OUT[expressed.scanned.ensembl]
##TRE.OUT <- targExpressedScans
##rm(TRE.OUT)

##targScans <- targExpressedScans

# load lps.mu
load(paste(exp.dir,"all.mus.objects.RData" ,sep="/"))
# load lps.ratios
load(paste(exp.dir,"all.ratios.objects.RData" ,sep="/"))

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



