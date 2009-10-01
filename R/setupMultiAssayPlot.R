## February 25, 2006

load("../annotations/representativeProbes.RData")
load("../expression_data/all.mus.objects.RData")


## Protein Westerns. Cytoplasmic and Nuclear
t1 <- read.table('TFwesternScoredCytoplasmic.tsv',as.is=TRUE,header=TRUE,row.names=1)
colnames(t1) <- c("min0","min30","min60","min120","min240","min360")
protein.cytoplasmic.western <- as.matrix(t1)
t1 <- read.table('TFwesternScoredNuclear.tsv',as.is=TRUE,header=TRUE,row.names=1)
colnames(t1) <- c("min0","min30","min60","min120","min240","min360")
protein.nuclear.western <- as.matrix(t1)
rm(t1)
t.western <- c(0,30,60,120,240,360)

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

##
## Summary: We now have
## protein.cytoplasmic.western ; t.western
## protein.nuclear.western ; t.western
## protein.il6.bound ; t.chip
## mRNA.realtime ; t.realtime
## mRNA.array ; t.array

t.all <- unique(sort(c(t.western,t.chip,t.realtime,t.array)))

## Scale each to its maximum? 
## Westerns, are already scaled

protein.il6.bound.max <- apply(protein.il6.bound,1,max)
protein.il6.bound <- protein.il6.bound/protein.il6.bound.max

mRNA.realtime.max <- apply(mRNA.realtime,1,max)
mRNA.realtime <- mRNA.realtime/mRNA.realtime.max

mRNA.array.max <- apply(mRNA.array,1,max)
mRNA.array <- mRNA.array/mRNA.array.max
