
rm(list=ls());

#########################
split1 <- function(splitme,splitchar='_'){ #parse binding site names by dropping "_" and everything behind it
 strsplit(splitme,split=splitchar,fixed=TRUE,extended=FALSE)[[1]][1];
}
#########################

#load in the data
load('Parsed.Scan.car.allMouseMarch2008.RData'); #load in the parsed TF binding site data
TRE.OUT <- TRE.OUT.car.diffexp
ensemblIDs.car <- names(TRE.OUT)

load("../Parsed.Scan.eset.allMouseMarch2008.RData")
TRE.OUT <- TRE.OUT.ensids
ensids <- names(TRE.OUT)

## Read all mouse matrices and lengths
gg <-read.table("~/macrophage/RShared/matrixLengths",as.is=TRUE); matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg); 
mouseMatrices <- names(matrixLength)

## Load file of threshold for each matrix 
load("~/macrophage/RShared/matrixThresholdA.RData")
matrixThreshold <- matrixThresholdA; rm(matrixThresholdA); 

## load familyMap, a mapping from vector of matrix names to family strings
load("~/macrophage/RShared/familyMapVT2.RData")
addendum <- "ISRE"
names(addendum) <- "ISRE_01"
familyMapVT2 <- c(familyMapVT2,addendum)
familyMap <- familyMapVT2; rm(familyMapVT2); 
familyNames <- unique(sort(as.character(familyMap)))

# load ensemblIDsOfEntrezIDs
load('/Users/thorsson/data/TransScan/allMouseUpstreamMarch2008/allMouseMappings.RData')

## Load required functions
source("/Users/thorsson/macrophage/RShared/tallyUtilities.R")

entrezIDs.car <- as.character(entrezIDofEnsemblID[ensemblIDs.car])

entrezIDs <- unique(unlist(entrezIDofEnsemblID[ensids]))
## otherwise known as expressed.scanned.egid ( no hyphens) 

load("~/macrophage/TFinfluence/R/tteMaps.RData")


save.image()
