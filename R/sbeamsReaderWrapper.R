

sbeams.data.directory  <- file.path(Sys.getenv("TFINF"),"expression_data")
## Formerly
## sbeams.data.directory <- "/users/thorsson/macrophage/AffyArray/AllSBEAMS"

## Handy quick way to get column headers. Use inline.
## colheaders <- names(read.table(sbeams.file,header=TRUE,sep='\t',nrow=1))
##  colSuffs <- mapply(tailEnd,colheaders,rep(nsuff,length(colheaders)))
##  needed.conds <- names(which(colSuffs==desiredSuffix))

## Read in CpG data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsCpG.tsv",sep="/")
sbeamsCondOrder <- c(2,3,4,5,1)
cond.names.short <- c("min20","min40","min60","min80","min120")
desiredSuffixes <- c("log10_ratio","mu_x","false_discovery_rate")

cpg.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
cpg.ratios <- cpg.mats[[desiredSuffixes[1]]]
cpg.mus <- cpg.mats[[desiredSuffixes[2]]]
cpg.fdr <- cpg.mats[[desiredSuffixes[3]]]


## Read in R848 data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsR848.tsv",sep="/")
sbeamsCondOrder <- c(1,2,3,4,5)
cond.names.short <- c("min20","min40","min60","min80","min120")
desiredSuffixes <- c("log10_ratio","mu_y","false_discovery_rate")

r848.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
r848.ratios <- r848.mats[[desiredSuffixes[1]]]
r848.mus <- r848.mats[[desiredSuffixes[2]]]
r848.fdr <- r848.mats[[desiredSuffixes[3]]]


## Read in PolyIC data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsPolyIC.tsv",sep="/")
sbeamsCondOrder <- c(1,2,3,4,5)
cond.names.short <- c("min20","min40","min60","min80","min120")
desiredSuffixes <- c("log10_ratio","mu_y","false_discovery_rate")

polyIC.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
polyIC.ratios <- polyIC.mats[[desiredSuffixes[1]]]
polyIC.mus <- polyIC.mats[[desiredSuffixes[2]]]
polyIC.fdr <- polyIC.mats[[desiredSuffixes[3]]]

## Read in PAM2 data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsPAM2.tsv",sep="/")
sbeamsCondOrder <- c(1,2,3,4,5)
cond.names.short <- c("min20","min40","min60","min80","min120")
desiredSuffixes <- c("log10_ratio","mu_y","false_discovery_rate")

pam2.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
pam2.ratios <- pam2.mats[[desiredSuffixes[1]]]
pam2.mus <- pam2.mats[[desiredSuffixes[2]]]
pam2.fdr <- pam2.mats[[desiredSuffixes[3]]]

## Read in PAM3 data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsPAM3.tsv",sep="/")
sbeamsCondOrder <- c(1,2,3,4,5)
cond.names.short <- c("min20","min40","min60","min80","min120")
desiredSuffixes <- c("log10_ratio","mu_y","false_discovery_rate")

pam3.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
pam3.ratios <- pam3.mats[[desiredSuffixes[1]]]
pam3.mus <- pam3.mats[[desiredSuffixes[2]]]
pam3.fdr <- pam3.mats[[desiredSuffixes[3]]]

## Read in MyD88 KO data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsMyD88KO.tsv",sep="/")

sbeamsCondOrder <- c(5,6,4) ## myD88 KO LPS 
cond.names.short <- c("min20","min60","min120")
desiredSuffixes <- c("log10_ratio","mu_x","false_discovery_rate")
myd88KOlps.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
myd88KOlps.ratios <- myd88KOlps.mats[[desiredSuffixes[1]]]
myd88KOlps.mus <- myd88KOlps.mats[[desiredSuffixes[2]]]
myd88KOlps.fdr <- myd88KOlps.mats[[desiredSuffixes[3]]]

sbeamsCondOrder <- c(8,9,7) ## MyD88 KO PAM3 
cond.names.short <- c("min20","min60","min120")
desiredSuffixes <- c("log10_ratio","mu_x","false_discovery_rate")
myd88KOpam3.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
myd88KOpam3.ratios <- myd88KOpam3.mats[[desiredSuffixes[1]]]
myd88KOpam3.mus <- myd88KOpam3.mats[[desiredSuffixes[2]]]
myd88KOpam3.fdr <- myd88KOpam3.mats[[desiredSuffixes[3]]]

sbeamsCondOrder <- c(11,12,10) ## MyD88 KO POLYIC 
cond.names.short <- c("min20","min60","min120")
desiredSuffixes <- c("log10_ratio","mu_x","false_discovery_rate")
myd88KOpolyIC.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
myd88KOpolyIC.ratios <- myd88KOpolyIC.mats[[desiredSuffixes[1]]]
myd88KOpolyIC.mus <- myd88KOpolyIC.mats[[desiredSuffixes[2]]]
myd88KOpolyIC.fdr <- myd88KOpolyIC.mats[[desiredSuffixes[3]]]



## Read in PAM3 PolyIC  data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsPAM3PolyIC.tsv",sep="/")
names(read.table(sbeams.file,header=TRUE,sep='\t',nrow=1))

sbeamsCondOrder <- c(5,6,3) 
cond.names.short <- c("min20","min60","min120")
desiredSuffixes <- c("log10_ratio","mu_x","false_discovery_rate")
pam3polyIC.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)
pam3polyIC.ratios <- pam3polyIC.mats[[desiredSuffixes[1]]]
pam3polyIC.mus <- pam3polyIC.mats[[desiredSuffixes[2]]]
pam3polyIC.fdr <- pam3polyIC.mats[[desiredSuffixes[3]]]

## Read in TRIF  KO data from SBEAMS file 
sbeams.file <- paste(sbeams.data.directory,"sbeamsTRIFKO.tsv",sep="/")

sbeamsCondOrder <- c(2,1) ## trif KO LPS 
cond.names.short <- c("min60","min120")
desiredSuffixes <- c("log10_ratio","mu_x")
trifKOlps.mats <- getMatricesFromSbeamsFile( sbeams.file, sbeamsCondOrder,desiredSuffixes,cond.names.short)

trifKOlps.ratios <- trifKOlps.mats[[desiredSuffixes[1]]]
trifKOlps.mus <- trifKOlps.mats[[desiredSuffixes[2]]]

#### Ratio cube for wildtype

ps.max <- 45037

ratioCube <- rep(NA,6*ps.max*5)
dim(ratioCube) <- c(6,ps.max,5)
dimnames(ratioCube)[[3]] <- c("min20","min40","min60","min80","min120")
dimnames(ratioCube)[[2]] <- probeset[1:ps.max]
dimnames(ratioCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

ratioCube["LPS",,] <- lps.ratios[1:ps.max,1:5]
ratioCube["PAM2",,] <- pam2.ratios[1:ps.max,1:5]
ratioCube["PAM3",,] <- pam3.ratios[1:ps.max,1:5]
ratioCube["PolyIC",,] <- polyIC.ratios[1:ps.max,1:5]
ratioCube["R848",,] <- r848.ratios[1:ps.max,1:5]
ratioCube["CpG",,] <- cpg.ratios[1:ps.max,1:5]

