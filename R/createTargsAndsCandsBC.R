
##
## January 2008
##

## Updated to included GSEA version of enrichment analysis for enrp and enrs.
## all other pdna models remain the same.
## Jan 16, 2009 was as follows:
## Ran only enrp and enrs portions
## Saved them as pdna.enrp.01.RData, pdna.enrp.05.RData, pdna.enrs.001.RData

## Outside this code, 
## Loaded pdnaModels.18June2008.RData, and did a replacement of the above modules, then did the save at the end of this code

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
interact.dir <- file.path(Sys.getenv("TFINF"),"interaction_data")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
r.dir <- file.path(Sys.getenv("TFINF"),"R")

source(paste(r.dir,"utilitiesMeta.R",sep="/"))
source(paste(r.dir,"utilitiesInteractions.R",sep="/"))

load(paste(annot.dir,"TFcategories.RData",sep="/"))
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"allMouseMappings.August2009.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))

load(paste(interact.dir,"pdnaModels.10May2010.RData",sep="/"))
load(paste(seq.dir,"sigPairedSitesBC.RData",sep="/"))

### Enriched pairs
### Choices: enrichment p-value, tfs to use, whether to use interactome
##

pc1 <- list()
pc1$metampairs <- metampairs.bc
pc1$metapvals <- metapvals.bc
pc1$metams <- metams.bc

### Feb 2007 
pair.cutoff <- 0.01
pc2.01 <- filterPCbyPval ( pc1, pair.cutoff )

tfsubset <- setdiff(c(bigUps,midUps,bigDowns,midDowns),repProbes.cname[c("E2f5","Gabpa")])
## E2f5 comes up only at 24 hrs
## Gabpa has a moderate, late change . Can be included, but not the cream-of-the-crop predictions
pc2.05 <- filterPCbyPval ( pc1, 0.05 )

##Load TF network distance
load("NetworkDistance/TFNetDist.RData") ### this loads tf.dist.cn 
noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information
noConnections <- c(noConnections,"Fosl2") ## May 2008. Can't quite figure out: Is this a "new" one?
noConnections <- c(noConnections,"Stat2") ## Sep 2009. Can't quite figure out: Is this a "new" one? probably from ISRE
noConnections <- c(noConnections,"Atf1") ## May 2010

## May 2008. Repairs needed, as we indexed by gene symbol and not entrez ID
tf.names.old <- rownames(tf.dist.cn)
tf.names.new <- replace(tf.names.old,match(c("Isgf3g","Tgif","Zfpn1a1"),tf.names.old),c("Irf9","Tgif1","Ikzf1"))
rownames(tf.dist.cn) <- tf.names.new
colnames(tf.dist.cn) <- tf.names.new

## January 2008
pdna.enrp.01 <- createTFsetFromPairs(pc2.01,tfsubset=tfsubset)
pdna.enrp.05 <- createTFsetFromPairs(pc2.05,tfsubset=tfsubset)

##save(pdna.enrp.01,file="pdna.enrp.01.RData")
##save(pdna.enrp.05,file="pdna.enrp.05.RData")

## Enriched singles
## Choices: enrichemnt p-value, tfs to use
singles.cutoff <- 0.001
sc1 <- list()
sc1$metamsingles <- metamsingles.bc
sc1$metapvals.singles <- metapvals.singles.bc
sc1$metams.singles <- metams.singles.bc
sc2 <-  filterSCbyPval ( sc1, singles.cutoff )
pdna.enrs.001 <- createTFsetFromSingles(sc2,tfsubset=tfsubset)
##save(pdna.enrs.001,file="pdna.enrs.001.RData")

ofile <- paste(interact.dir,"pdnaModels.2Aug2010.RData",sep="/")
save(pdna.curated, pdna.enrs.001, pdna.enrp.01, pdna.enrp.05, pdna.hs.et.5,pdna.hs.et.1,pdna.hs.et.01,pdna.hs.et.05,pdna.hs.et.001, file=ofile)
