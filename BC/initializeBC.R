
## Create .RData file for runnnig multiple scripts in BC directory

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
r.dir <- file.path(Sys.getenv("TFINF"),"R")

load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(annot.dir,"allMouseMappings.August2010.RData",sep="/"))

source(paste(r.dir,"tallyUtilities.R",sep="/"))
source("./utilitiesBinaryClassification.R")

load(paste(seq.dir,"Parsed.Scan.eset.RData",sep="/"))
TRE.OUT <- TRE.OUT.eset
ensids <- names(TRE.OUT)
entrezIDs <- unique(unlist(entrezIDofEnsemblID[ensids]))
## otherwise known as expressed.scanned.egid ( no hyphens) 

redundancy <- unlist(lapply(ensemblIDsOfEntrezIDs[entrezIDs],length))
good.eids <- names(which(redundancy==1)) ## Gene IDs with non-redundant ensembl IDs
good.ensids <- as.character(unlist(ensemblIDsOfEntrezIDs[good.eids]))

ensids <- good.ensids 
eids <- good.eids
entrezIDs <- eids
N <- length(entrezIDs) ## Universe size

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

## Expression neighbors
eids.all <- entrezIDofEnsemblID[good.ensids]

load(paste(annot.dir,"representativeProbes.RData",sep="/"))
psois.all  <- as.character(repProbes.ncbiID[eids.all])
## The simpler version below also is identical (July2010)
##psois.all <- paste(eids.all,"_at",sep="")  ## improve this

load(paste(exp.dir,"all.mus.objects.RData",sep="/"))
t.index.min <- 1
t.index.max <- 11
cordistH <- 1-cor(t(lps.mus[psois.all,t.index.min:t.index.max]))

maxn <-  600 
mincor <- 0.1
system.time( nbrs <- getNbrs(psois=psois.all,cloud=rownames(cordistH), distmat=cordistH, threshold=mincor))

## Use this if running command line version
save.image()
