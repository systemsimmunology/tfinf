
## Preparation for BC runs

exp.dir <- file.path(Sys.getenv("TFINF"),"data/expression")
seq.dir <- file.path(Sys.getenv("TFINF"),"data/motifscan_input")
aux.dir <- file.path(Sys.getenv("TFINF"),"aux")
utils.dir <- file.path(Sys.getenv("TFINF"),"utils")
out.dir <- file.path(Sys.getenv("TFINF"),"data/process_output")

load(paste(exp.dir,"representativeProbes.RData",sep="/"))
load(paste(exp.dir,"annotation.objects.RData",sep="/"))
load(paste(out.dir,"tteMaps.RData",sep="/"))
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(seq.dir,"allMouseMappings.RData",sep="/"))
source(paste(utils.dir,"utilitiesBinaryClassification.R",sep="/"))

load(paste(out.dir,"Parsed.Scan.eset.RData",sep="/"))
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
gg <-read.table(paste(aux.dir,"matrixLengths" ,sep="/"),as.is=TRUE);
matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg);
mouseMatrices <- names(matrixLength)

 
## load familyMap, a mapping from vector of matrix names to family strings
load(paste(aux.dir,"familyMapVT3.RData" ,sep="/"))
familyMap <- familyMapVT3; rm(familyMapVT3);
familyNames <- unique(sort(as.character(familyMap)))

## Load file of threshold for each matrix
load(paste(aux.dir,"matrixThresholdA.RData" ,sep="/"))
matrixThreshold <- matrixThresholdA; 

## Expression neighbors
eids.all <- entrezIDofEnsemblID[good.ensids]

load(paste(exp.dir,"representativeProbes.RData",sep="/"))
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

save.image(file=paste(out.dir,"initializeBC.RData",sep="/"))
