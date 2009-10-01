

coi <- "Atf3"

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
interact.dir <- file.path(Sys.getenv("TFINF"),"interaction_data")
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
ddata.dir <- file.path(Sys.getenv("TFINF"),"derived_data")

source("./multiAssayPlot.R")
source("./newFunctions.R")

load(paste(exp.dir,"scaled.mus.objects.RData",sep="/"))
load(paste(exp.dir,"all.mus.objects.RData",sep="/"))
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"annotation.objects.RData",sep="/"))

t.array <- c(0,20,40,60,80,120,240,360,480,1080,1440)
t.realtime <- c(0,60,120,240,360,1080)
t.chip <- c(0,60,120,240,360)
t.western <- c(0,30,60,120,240,360)
t.all <- unique(sort(c(t.western,t.chip,t.realtime,t.array)))

t1 <- read.table('TFwesternScoredCytoplasmic.tsv',as.is=TRUE,header=TRUE,row.names=1)
colnames(t1) <- c("min0","min30","min60","min120","min240","min360")
protein.cytoplasmic.western <- as.matrix(t1)
t1 <- read.table('TFwesternScoredNuclear.tsv',as.is=TRUE,header=TRUE,row.names=1)
colnames(t1) <- c("min0","min30","min60","min120","min240","min360")
protein.nuclear.western <- as.matrix(t1)
rm(t1)

## Chromatin IP
t.chip <- c(0,60,120,240,360)
atf3.chip <- c(1, 18, 51, 102, 110)
rel.chip <- c(3, 42, 22, 13, 4)
control.chip  <- c(143, 124, 126, 137, 140 )
norm.control <- control.chip/mean(control.chip)
atf3.chip <- c(1, 18, 51, 102, 110)/ norm.control
rel.chip <- c(3, 42, 22, 13, 4)/ norm.control
protein.il6.bound <- rbind(atf3.chip,rel.chip)
colnames(protein.il6.bound) <- paste("min",as.character(t.chip),sep="")
rownames(protein.il6.bound) <- c("Atf3","Rel")


## mRNA realtime
t.realtime <- c(0,60,120,240,360,1080)
atf3.realtime <- c(0.041235,0.486327,0.08362,0.190782,0.108067,0.070805)
rel.realtime <- c(0.00018,0.006479,0.002613,0.001456,0.001236,0.000284)
mRNA.realtime <- rbind(rel.realtime,atf3.realtime)
colnames(mRNA.realtime) <- paste("min",as.character(t.realtime),sep="")
rownames(mRNA.realtime) <- c("Rel","Atf3")
  
## mRNA array 
t.array <- c(0,20,40,60,80,120,240,360,480,1080,1440)
cnois <- rownames(protein.cytoplasmic.western)
psois <- repProbes.cname[cnois]
mRNA.array <- lps.mus[psois,]
rownames(mRNA.array) <- cnois

psoi <- repProbes.cname[coi]

zz <- rbind(lps.mat.max1[psoi,1:6],
      pam2.mat.max1[psoi,1:6],
      pam3.mat.max1[psoi,1:6],
      polyIC.mat.max1[psoi,1:6],
      r848.mat.max1[psoi,1:6],
      cpg.mat.max1[psoi,1:6])

x11()
par(mfrow=c(3,1))
paiWT(psoi)

multiAssayPlot(coi,120)
      
pNLA(psoi)
