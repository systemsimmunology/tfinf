
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

source(file.path(Sys.getenv("UTILS_DIR"),"utilitiesMeta.R"))
source(file.path(Sys.getenv("UTILS_DIR"),"utilitiesInteractions.R"))
source(file.path(Sys.getenv("UTILS_DIR"),"utilitiesHitMat.R"))

load(file.path(Sys.getenv("PO_DIR"),"TFcategories.RData"))
load(file.path(Sys.getenv("EXP_DIR"),"representativeProbes.RData"))
load(file.path(Sys.getenv("EXP_DIR"),"annotation.objects.RData"))
load(file.path(Sys.getenv("SEQ_DIR"),"allMouseMappings.RData"))
load(file.path(Sys.getenv("PO_DIR"),"tteMaps.RData"))
load(file.path(Sys.getenv("AUX_DIR"),"equalizingThresholds.RData"))
load(file.path(Sys.getenv("PO_DIR"),"sigPairedSitesBC.RData"))
load(file.path(Sys.getenv("PO_DIR"),"pdna.curated.RData"))
load(file.path(Sys.getenv("PO_DIR"),"Parsed.Scan.eset.RData"))
targScans <- TRE.OUT.eset

tfsubset <- setdiff(c(bigUps,midUps,bigDowns,midDowns),repProbes.cname[c("E2f5","Gabpa")])
## E2f5 comes up only at 24 hrs
## Gabpa has a moderate, late change . Can be included, but not the cream-of-the-crop predictions

### High scoring individual hits
targScans.et.0.5 <- getScoreFilteredTE( targScans, et.0.5 )
targScans.et.0.1 <- getScoreFilteredTE( targScans, et.0.1 )
targScans.et.0.05 <- getScoreFilteredTE( targScans, et.0.05 )
targScans.et.0.01 <- getScoreFilteredTE( targScans, et.0.01 )
targScans.et.0.001 <- getScoreFilteredTE( targScans, et.0.001 )
pdna.hs.et.5 <- getTFsForTE(targScans.et.0.5,tfsubset=tfsubset)
pdna.hs.et.1 <- getTFsForTE(targScans.et.0.1,tfsubset=tfsubset )
pdna.hs.et.05 <- getTFsForTE(targScans.et.0.05,tfsubset=tfsubset)
pdna.hs.et.01 <- getTFsForTE(targScans.et.0.01,tfsubset=tfsubset )
pdna.hs.et.001 <- getTFsForTE(targScans.et.0.001,tfsubset=tfsubset )

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
pc2.05 <- filterPCbyPval ( pc1, 0.05 )

##Load TF network distance
load(file.path(Sys.getenv("AUX_DIR"),"TFNetDist.RData")) ### this loads tf.dist.cn 
noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information
noConnections <- c(noConnections,"Fosl2") ## May 2008. Can't quite figure out: Is this a "new" one?
noConnections <- c(noConnections,"Stat2") ## Sep 2009. Can't quite figure out: Is this a "new" one? probably from ISRE
noConnections <- c(noConnections,"Atf1") ## May 2010
noConnections <- c(noConnections,"Smad2","Bach1") ## Nov 2010

## May 2008. Repairs needed, as we indexed by gene symbol and not entrez ID
tf.names.old <- rownames(tf.dist.cn)
tf.names.new <- replace(tf.names.old,match(c("Isgf3g","Tgif","Zfpn1a1"),tf.names.old),c("Irf9","Tgif1","Ikzf1"))
rownames(tf.dist.cn) <- tf.names.new
colnames(tf.dist.cn) <- tf.names.new

## January 2008
pdna.enrp.01 <- createTFsetFromPairs(pc2.01[1:3],tfsubset=tfsubset,noConnections=noConnections)
pdna.enrp.05 <- createTFsetFromPairs(pc2.05,tfsubset=tfsubset,noConnections=noConnections)

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

ofile <- file.path(Sys.getenv("PO_DIR"),"pdnaModels.RData")
save(file=ofile,pdna.curated, pdna.enrs.001, pdna.enrp.01, pdna.enrp.05, pdna.hs.et.5,pdna.hs.et.1,pdna.hs.et.01,pdna.hs.et.05,pdna.hs.et.001)

