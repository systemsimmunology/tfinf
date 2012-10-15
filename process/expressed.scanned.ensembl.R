##
## load annotations, mappings, probe set lists
##
exp.dir <- file.path(Sys.getenv("TFINF"),"data/expression")
seq.dir <- file.path(Sys.getenv("TFINF"),"data/motifscan_input")
out.dir <- file.path(Sys.getenv("TFINF"),"data/process_output")

load(paste(exp.dir,"representativeProbes.RData",sep="/"))
load(paste(exp.dir,"annotation.objects.RData",sep="/"))
load(paste(exp.dir,"all.ps.list.objects.RData",sep="/"))
## Read Ensembl-Entrez Mappings 
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(seq.dir,"allMouseMappings.RData",sep="/"))

## Some of the scanned genes are not on the array. Remove those.
eids.on.array <- names(repProbes.ncbiID)
eids.scanned <- names(ensemblIDsOfEntrezIDs) ## entrezIDs
eids.on.both <- intersect(eids.on.array,eids.scanned)
scanned.ps <- as.character(repProbes.ncbiID[eids.on.both])
## these last ones are unique: probe sets that are not non-repProbes are not included

## June 2008. Hyphenated gene names were not incorporated correctly in scanning objects
## Sep 2009: May mean that they are not in TRE.OUT, e.g. entrez:Krtap5-5,114666=ensembl:"73785"=
## Sep 2009: TRE.OUT generator (read.parsed.files.R) now accomodates hyphenated names 
##hasCharacter <- function( instring, querychar ){querychar %in% strsplit(instring,split="")[[1]]} ## does R really not have this?
##dashed.genenames.ps <- names(which(sapply(cname.compare[scanned.ps],hasCharacter,"-")))
##scanned.ps <- setdiff(scanned.ps,dashed.genenames.ps)

## June 2008: 4714 of these:
## September 2009: 4713 of these:
expressed.scanned.ps <- scanned.ps[scanned.ps %in% lps.full.ps.sig]
expressed.scanned.egid <- as.character(ncbiID[expressed.scanned.ps])
## June 2008 :4869 of these
## Sep 2009 :4742 of these. Sampled some differences and they are consistent with ensembl website
expressed.scanned.ensembl <- as.character(unlist(ensemblIDsOfEntrezIDs[expressed.scanned.egid]))

save(expressed.scanned.ensembl, file=paste(out.dir,"expressed.scanned.ensembl.RData",sep="/"))

