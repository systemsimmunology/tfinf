
rm(list=ls());

seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")

load(paste(seq.dir,"Parsed.Scan.Results.allMouseMay2010.dat",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))

## Read all mouse matrices and lengths
gg <-read.table(paste(annot.dir,"matrixLengths",sep="/"),as.is=TRUE); matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg); 
mouseMatrices <- names(matrixLength)

## Load file of threshold for each matrix 
load(paste(annot.dir,"matrixThresholdA.RData",sep="/"))
matrixThreshold <- matrixThresholdA; rm(matrixThresholdA); 

## load familyMap, a mapping from vector of matrix names to family strings
load(paste(annot.dir,"familyMapVT3.RData",sep="/"))
familyMap <- familyMapVT3; rm(familyMapVT3); 
familyNames <- unique(sort(as.character(familyMap)))

# load entrezIDs,ensemblIDsOfEntrezIDs
load(paste(annot.dir,"allMouseMappings.August2009.RData",sep="/"))
entrezIDs <- names(ensemblIDsOfEntrezIDs)
