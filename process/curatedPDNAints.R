## Transfac interactions, documented

source(paste(Sys.getenv("UTILS_DIR"),"selectionUtilities.R",sep="/"))
source(paste(Sys.getenv("UTILS_DIR"),"readInteractions.R",sep="/"))

load(paste(Sys.getenv("EXP_DIR"),"annotation.objects.RData",sep="/"))
load(paste(Sys.getenv("EXP_DIR"),"representativeProbes.RData",sep="/"))
load(paste(Sys.getenv("EXP_DIR"),"all.ps.list.objects.RData",sep="/"))
load(paste(Sys.getenv("PO_DIR"),"tteMaps.RData",sep="/"))

result <- readInteractions (file.path(Sys.getenv("AUX_DIR"),"n28.sif"))
pd.indices <- which(result$interactType=="Protein-DNA")

tfs.np <- result$partner1[pd.indices]
targets.np <- result$partner2[pd.indices]

tfs <- repProbes.np[tfs.np]
targets <- repProbes.np[targets.np] 

## lambda
expressed.set <- lps.full.ps.sig

tfs.exp <- as.character( tfs[(tfs %in% expressed.set )  &  (targets %in% expressed.set)])
targets.exp <- as.character( targets[(tfs %in% expressed.set )  &  (targets %in% expressed.set)])
## these are still paired, equal in length

## There are a few discrepancies revealed in nonzero setdiff(us(tfs.exp),tte)
## Stat1 and Stat3 have different repProbes plus we have 'new' TFs
## A few have a decent expression level eg. Bcl6
##cc[setdiff(us(tfs.exp),tte)]
##           1416224_at          1416244_a_at            1418114_at
##             "Zbtb17"               "Cnbp1"              "Rbpsuh"
##         1418275_a_at            1421818_at            1422742_at
##               "Elf2"                "Bcl6"              "Hivep1"
##         1425996_a_at          1426587_a_at            1426744_at
##            "Smarca3"               "Stat3"              "Srebf2"
##           1439336_at          1441956_s_at          1449947_s_at
##               "Tcf4"                 "---"               "Atbf1"
##         1450033_a_at          1451012_a_at            1459360_at
##              "Stat1"                "Csda"                 "Sp3"
##         1460279_a_at
##"Gtf2i /// LOC669007"

tfs.exp <- replace(tfs.exp,which(tfs.exp=="1450033_a_at"),"1440481_at") ## Stat1
tfs.exp <- replace(tfs.exp,which(tfs.exp=="1426587_a_at"),"1459961_a_at") ## Stat3
targets.exp <- replace(targets.exp,which(targets.exp=="1450033_a_at"),"1440481_at") ## Stat1
targets.exp <- replace(targets.exp,which(targets.exp=="1426587_a_at"),"1459961_a_at") ## Stat3

## Further weeding, so we only work with transfac.tfs.expressed
tfs.exp.2 <- as.character(tfs.exp[(tfs.exp %in% transfac.tfs.expressed)&(targets.exp %in% expressed.set)])
targets.exp.2 <- as.character(targets.exp[(tfs.exp %in% transfac.tfs.expressed )&(targets.exp %in% expressed.set)])
tfs.exp <- tfs.exp.2 ; rm(tfs.exp.2)
targets.exp <- targets.exp.2 ; rm(targets.exp.2)
## these are still paired, equal in length

 

tfs.total.exp <- tfs.exp
targets.total.exp <- targets.exp
 
sc.iws <- getIncomingWires(tfs.total.exp,targets.total.exp)
n.sc.iws <- length(sc.iws)
targs.total.exp.unique <- unique(sort(targets.total.exp))
n.max.preds <- numeric(length=length(targs.total.exp.unique))
names(n.max.preds) <- targs.total.exp.unique

## Count number of incoming
for ( targ in targs.total.exp.unique ){
  invars <- sc.iws[[targ]]
  n.max.predictors <- length(invars)
  n.max.preds[targ] <- n.max.predictors
}

## Convention with:
## targets: psois 
## tfs: cnames

pdna.curated <- list()
for ( targ in targs.total.exp.unique ){
  invars <- sc.iws[[targ]]
  invars.c <- as.character(cname.compare[invars])
  invars.c <- setdiff(invars.c,"---")
  pdna.curated[[targ]] <- invars.c
}

save(pdna.curated,file=file.path(Sys.getenv("PO_DIR"),"pdna.curated.RData"))
