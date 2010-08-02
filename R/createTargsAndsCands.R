
##
## September 2009
##

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
interact.dir <- file.path(Sys.getenv("TFINF"),"interaction_data")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")

source("./utilitiesMeta.R")
source("./utilitiesInteractions.R")
source("./utilitiesHitMat.R")
 
load(paste(annot.dir,"TFcategories.RData",sep="/"))
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"allMouseMappings.August2009.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))
load(paste(annot.dir,"equalizingThresholds.RData",sep="/"))
load(paste(annot.dir,"expressed.scanned.ensembl.RData",sep="/"))
load(paste(annot.dir,"TFcategories.RData",sep="/"))
load(paste(interact.dir,"pdna.curated.RData",sep="/"))
load(paste(seq.dir,"sigPairedSites.28Sep2009.RData",sep="/"))
load(paste(seq.dir,"Parsed.Scan.eset.allMouseSeptember2009.RData",sep="/"))

targScans <- TRE.OUT.eset

### Enriched pairs
### Choices: enrichment p-value, tfs to use, whether to use interactome
##
ensids.sansSSS <- intersect(names(metampairs),expressed.scanned.ensembl)
## the last set is equal to ensemblIDs.expressed.sansSSS
## September 2009: Appears this is no longer necessary
## However, the intersect removes ensids that didn't make it into metampairs
pc1 <- list()
pc1$metampairs <- metampairs[ensids.sansSSS]
pc1$metapvals <- metapvals[ensids.sansSSS]
pc1$metams <- metams[ensids.sansSSS]

### Feb 2007 
pair.cutoff <- 0.01
pc2.01 <- filterPCbyPval ( pc1, pair.cutoff )
tfsubset <- setdiff(c(bigUps,midUps,bigDowns,midDowns),repProbes.cname[c("E2f5","Gabpa")])
## E2f5 comes up only at 24 hrs
## Gabpa has a moderate, late change . Can be included, but not the cream-of-the-crop predictions
pc2.05 <- filterPCbyPval ( pc1, 0.05 )


##Load TF network distance
load(paste(interact.dir,"TFNetDist.RData",sep="/"))
           
noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information
noConnections <- c(noConnections,"Fosl2") ## May 2008. Can't quite figure out: Is this a "new" one?
noConnections <- c(noConnections,"Stat2") ## Sep 2009. Can't quite figure out: Is this a "new" one? probably from ISRE
noConnections <- c(noConnections,"Atf1") ## May 2010

## May 2008. Repairs needed, as we indexed by gene symbol and not entrez ID

tf.names.old <- rownames(tf.dist.cn)
tf.names.new <- replace(tf.names.old,match(c("Isgf3g","Tgif","Zfpn1a1","Tcfe3"),tf.names.old),c("Irf9","Tgif1","Ikzf1","Tcf3"))
rownames(tf.dist.cn) <- tf.names.new
colnames(tf.dist.cn) <- tf.names.new

pdna.enrp.01 <- createTFsetFromPairs(pc2.01,tfsubset=tfsubset)
pdna.enrp.05 <- createTFsetFromPairs(pc2.05,tfsubset=tfsubset)

## Enriched singles
## Choices: enrichemnt p-value, tfs to use
singles.cutoff <- 0.001
ensids.sansSSS <- intersect(names(metamsingles),expressed.scanned.ensembl)
sc1 <- list()
sc1$metamsingles <- metamsingles[ensids.sansSSS]
sc1$metapvals.singles <- metapvals.singles[ensids.sansSSS]
sc1$metams.singles <- metams.singles[ensids.sansSSS]
sc2 <-  filterSCbyPval ( sc1, singles.cutoff )
pdna.enrs.001 <- createTFsetFromSingles(sc2,tfsubset=tfsubset)


##These may be commented out because they are time-consuming
targScans.et.0.5 <- getScoreFilteredTE( targScans, et.0.5 )
targScans.et.0.1 <- getScoreFilteredTE( targScans, et.0.1 )
targScans.et.0.05 <- getScoreFilteredTE( targScans, et.0.05 )
targScans.et.0.01 <- getScoreFilteredTE( targScans, et.0.01 )
targScans.et.0.001 <- getScoreFilteredTE( targScans, et.0.001 )

save(targScans.et.0.5,file=paste(seq.dir,"targScans.et.0.5.RData",sep="/"))
save(targScans.et.0.1,file=paste(seq.dir,"targScans.et.0.1.RData",sep="/"))
save(targScans.et.0.05,file=paste(seq.dir,"targScans.et.0.05.RData",sep="/"))
save(targScans.et.0.01,file=paste(seq.dir,"targScans.et.0.01.RData",sep="/"))
save(targScans.et.0.001,file=paste(seq.dir,"targScans.et.0.001.RData",sep="/"))

baddies <- as.character(repProbes.cname[c("E2f5","Gabpa")])

tfsubset <- setdiff( c(bigUps,midUps,bigDowns,midDowns), baddies ) 

pdna.hs.et.5 <- getTFsForTE(targScans.et.0.5,tfsubset=tfsubset)
pdna.hs.et.1 <- getTFsForTE(targScans.et.0.1,tfsubset=tfsubset )
pdna.hs.et.05 <- getTFsForTE(targScans.et.0.05,tfsubset=tfsubset)
pdna.hs.et.01 <- getTFsForTE(targScans.et.0.01,tfsubset=tfsubset )
pdna.hs.et.001 <- getTFsForTE(targScans.et.0.001,tfsubset=tfsubset )

ofile <- paste(interact.dir,"pdnaModels.10May2010.RData",sep="/")
save(pdna.curated, pdna.enrs.001, pdna.enrp.01, pdna.enrp.05, pdna.hs.et.5,pdna.hs.et.1,pdna.hs.et.01,pdna.hs.et.05,pdna.hs.et.001, file=ofile )
