
##
## September 2009
##

source("utilitiesMeta.R")
source("utilitiesHitMat.R")
source("utilities.R")
load("TFcategories.RData")
load("~/macrophage/AffyArray/newAffy/representativeProbes.RData")
load("~/macrophage/AffyArray/newAffy/annotation.objects.RData")
load("~/macrophage/RShared/allMouseMappings.August2009.RData")
load("~/macrophage/TFinfluence/R/tteMaps.RData")
load("./pdna.curated.RData")


load("sigPairedSites.18Sep2009.RData")
load("~/data/TransScan/RcodeFromDan/Parsed.Scan.eset.allMouseSeptember2009.RData")
##load("/proj/ilyalab/Vesteinn/data/TransScan/RcodeFromDan/Parsed.Scan.eset.allMouseSeptember2009.RData")

targScans <- TRE.OUT.eset

### Enriched pairs
### Choices: enrichment p-value, tfs to use, whether to use interactome
##
expressed.scanned.ensembl <- as.character(read.table("expressed.scanned.ensembl")$V1)
ensids.sansSSS <- intersect(names(metampairs),expressed.scanned.ensembl)
## the last set is equal to ensemblIDs.expressed.sansSSS
## September 2009: Appears this is no longer necessary
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
load("NetworkDistance/TFNetDist.RData") ### this loads tf.dist.cn 
noConnections <- c("Egr3","Mafb","Foxj2","Zfp161") ## TFs for which we have no connection information
noConnections <- c(noConnections,"Fosl2") ## May 2008. Can't quite figure out: Is this a "new" one?
noConnections <- c(noConnections,"Stat2") ## Sep 2009. Can't quite figure out: Is this a "new" one? probably from ISRE

## May 2008. Repairs needed, as we indexed by gene symbol and not entrez ID

tf.names.old <- rownames(tf.dist.cn)
tf.names.new <- replace(tf.names.old,match(c("Isgf3g","Tgif","Zfpn1a1"),tf.names.old),c("Irf9","Tgif1","Ikzf1"))
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

load("~/macrophage/RShared/equalizingThresholds.RData")
##These may be commented out because they are time-consuming
targScans.et.0.5 <- getScoreFilteredTE( targScans, et.0.5 )
targScans.et.0.1 <- getScoreFilteredTE( targScans, et.0.1 )
targScans.et.0.05 <- getScoreFilteredTE( targScans, et.0.05 )
targScans.et.0.01 <- getScoreFilteredTE( targScans, et.0.01 )
targScans.et.0.001 <- getScoreFilteredTE( targScans, et.0.001 )

save(targScans.et.0.5,file="targScans.et.0.5.RData")
save(targScans.et.0.1,file="targScans.et.0.1.RData")
save(targScans.et.0.05,file="targScans.et.0.05.RData")
save(targScans.et.0.01,file="targScans.et.0.01.RData")
save(targScans.et.0.001,file="targScans.et.0.001.RData")

baddies <- as.character(repProbes.cname[c("E2f5","Gabpa")])

tfsubset <- setdiff( c(bigUps,midUps,bigDowns,midDowns), baddies ) 

pdna.hs.et.5 <- getTFsForTE(targScans.et.0.5,tfsubset=tfsubset)
pdna.hs.et.1 <- getTFsForTE(targScans.et.0.1,tfsubset=tfsubset )
pdna.hs.et.05 <- getTFsForTE(targScans.et.0.05,tfsubset=tfsubset)
pdna.hs.et.01 <- getTFsForTE(targScans.et.0.01,tfsubset=tfsubset )
pdna.hs.et.001 <- getTFsForTE(targScans.et.0.001,tfsubset=tfsubset )

save(pdna.curated, pdna.enrs.001, pdna.enrp.01, pdna.enrp.05, pdna.hs.et.5,pdna.hs.et.1,pdna.hs.et.01,pdna.hs.et.05,pdna.hs.et.001, file="pdnaModels.22September2009.RData" )
