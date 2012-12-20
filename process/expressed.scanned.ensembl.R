##
## load annotations, mappings, probe set lists
##
load(paste(Sys.getenv("EXP_DIR"),"representativeProbes.RData",sep="/"))
load(paste(Sys.getenv("EXP_DIR"),"annotation.objects.RData",sep="/"))
load(paste(Sys.getenv("EXP_DIR"),"all.ps.list.objects.RData",sep="/"))
## Read Ensembl-Entrez Mappings 
## Read entrezIDofEnsemblID, ensemblIDsOfEntrezIDs
load(paste(Sys.getenv("SEQ_DIR"),"allMouseMappings.RData",sep="/"))
load(paste(Sys.getenv("SEQ_DIR"),"Parsed.Scan.Results.dat",sep="/"))

eids.on.array <- names(repProbes.ncbiID)
ensids.scanned <- names(TRE.OUT)
eids.scanned <- unique(as.character(entrezIDofEnsemblID[ensids.scanned]))
eids.on.both <- intersect(eids.on.array,eids.scanned)
scanned.ps <- as.character(repProbes.ncbiID[eids.on.both])
## these last ones are unique: probe sets that are not non-repProbes are not included

## June 2008. Hyphenated gene names were not incorporated correctly in scanning objects
## Sep 2009: May mean that they are not in TRE.OUT, e.g. entrez:Krtap5-5,114666=ensembl:"73785"=
## Sep 2009: TRE.OUT generator (read.parsed.files.R) now accomodates hyphenated names 
##hasCharacter <- function( instring, querychar ){querychar %in% strsplit(instring,split="")[[1]]} ## does R really not have this?
##dashed.genenames.ps <- names(which(sapply(cname.compare[scanned.ps],hasCharacter,"-")))
##scanned.ps <- setdiff(scanned.ps,dashed.genenames.ps)


## Restrict to differentially expressed genes
## June 2008: 4714 of these:
## September 2009: 4713 of these:
## October 2010: 4943 of these 
expressed.scanned.ps <- scanned.ps[scanned.ps %in% lps.full.ps.sig]
expressed.scanned.egid <- as.character(ncbiID[expressed.scanned.ps])


## June 2008 :4869 of these
## Sep 2009 :4742 of these. Sampled some differences and they are consistent with ensembl website
## October 2010: 4944 of these 
expressed.scanned.ensembl <- as.character(unlist(ensemblIDsOfEntrezIDs[expressed.scanned.egid]))
TRE.OUT.eset <- TRE.OUT[expressed.scanned.ensembl]

save(expressed.scanned.ensembl, file=paste(Sys.getenv("PO_DIR"),"expressed.scanned.ensembl.RData",sep="/"))
save(TRE.OUT.eset,file=paste(Sys.getenv("PO_DIR"),"Parsed.Scan.eset.RData",sep="/"))



