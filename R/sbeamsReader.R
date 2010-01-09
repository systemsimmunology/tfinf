
fpe <- file.path(Sys.getenv("TFINF"),"expression_data")

data.mus <- read.table(paste(fpe,"sbeams.mus",sep="/"),header=TRUE,sep='\t')
data.mus.min0 <- read.table(paste(fpe,"sbeams.minute0.mus",sep="/"),header=TRUE,sep='\t')
data.fdr <- read.table(paste(fpe,"sbeams.fdr",sep="/"),header=TRUE,sep='\t')
data.ratios <- read.table(paste(fpe,"sbeams.log10_ratio",sep="/"),header=TRUE,sep='\t')
data.lambdas <- read.table(paste(fpe,"LPS.lambdas",sep="/"),header=TRUE,sep='\t')

if ( length(which(as.character(data.ratios[["reporter_name"]]) != as.character(data.fdr[["reporter_name"]]))) != 0 ){
  cat("Error: Ratio and Fdr Rows do not match\n")
}
       
if ( length(which(as.character(data.ratios[["reporter_name"]]) != as.character(data.mus[["reporter_name"]]))) != 0 ){
  cat("Error: Ratios and mus rows do not match \n")
}


if ( length(which(as.character(data.ratios[["reporter_name"]]) != as.character(data.mus.min0[["reporter_name"]]))) != 0 ){
  cat("Error: Ratios and mus time 0  rows do not match \n")
}


### Identifiers

attach(data.mus)
gene.ids.all <- as.character(reporter_name)
probeset <- as.character(reporter_name)
commonName <- as.character(common_name)
np <- as.character(canonical_name)
ncbiID <- as.character(second_name)
detach(data.mus)

names(np) <- probeset
names(commonName) <- probeset
names(ncbiID) <- probeset


gene.names.all <- commonName



## Slice out the AFFY- controls

keep.indexes <- 1:45037

probeset <- probeset[keep.indexes]
commonName <- commonName[keep.indexes]
np <- np[keep.indexes]
ncbiID <- ncbiID[keep.indexes]

gene.ids.all <- gene.ids.all[keep.indexes]
gene.names.all <- gene.names.all[keep.indexes]



## test again

if ( length(which(as.character(data.lambdas[["Probesets"]]) != probeset)) != 0 ){
  cat("Error: Lambda matrix incosistent with Probeset order  \n")
}


lps.lambdas <- as.matrix(data.lambdas[,3:10])
rownames(lps.lambdas) <- probeset

## Condition names

cond.names.long <- c("LPS_020_min_vs_LPS_0_min__mu_x",
"LPS_040_min_vs_LPS_0_min__mu_x",
"LPS_060_min_vs_LPS_0_min__mu_x",
"LPS_080_min_vs_LPS_0_min__mu_x",
"LPS_120_min_vs_LPS_0_min__mu_x",
"LPS_24_hr_vs_LPS_0_min__mu_x",
"LPS_240_min_vs_LPS_0_min__mu_x",
"LPS_8_hr_vs_LPS_0_min__mu_x")

cond.names.long.timeordered <- c("LPS_020_min_vs_LPS_0_min__mu_x",
"LPS_040_min_vs_LPS_0_min__mu_x",
"LPS_060_min_vs_LPS_0_min__mu_x",
"LPS_080_min_vs_LPS_0_min__mu_x",
"LPS_120_min_vs_LPS_0_min__mu_x",
"LPS_240_min_vs_LPS_0_min__mu_x",
"LPS_8_hr_vs_LPS_0_min__mu_x",
"LPS_24_hr_vs_LPS_0_min__mu_x")

cond.names.long.timeordered <- c("LPS_020_min_vs_LPS_0_min",
"LPS_040_min_vs_LPS_0_min",
"LPS_060_min_vs_LPS_0_min",
"LPS_080_min_vs_LPS_0_min",
"LPS_120_min_vs_LPS_0_min",
"LPS_240_min_vs_LPS_0_min",
"LPS_8_hr_vs_LPS_0_min",
"LPS_24_hr_vs_LPS_0_min")

cond.names.short <- c("min20","min40","min60","min80","min120","hr4","hr8","hr24")

## Numerical data

lps.mus <- apply((as.matrix(data.mus))[,paste(cond.names.long.timeordered,"__mu_x",sep="")],2,as.numeric)
lps.mus.0 <- as.numeric(data.mus.min0[["LPS_020_min_vs_LPS_0_min__mu_y"]])
lps.mus <- cbind(lps.mus.0,lps.mus)                                      

colnames(lps.mus) <- c("min0",cond.names.short)
lps.mus <- lps.mus[keep.indexes,]
rownames(lps.mus) <- probeset

lps.fdr <- apply((as.matrix(data.fdr))[,paste(cond.names.long.timeordered,"__false_discovery_rate",sep="")],2,as.numeric)
colnames(lps.fdr) <- cond.names.short
lps.fdr <- lps.fdr[keep.indexes,]
rownames(lps.fdr) <- probeset

lps.ratios <- apply((as.matrix(data.ratios))[,paste(cond.names.long.timeordered,"__log10_ratio",sep="")],2,as.numeric)
colnames(lps.ratios) <- cond.names.short
lps.ratios <- lps.ratios[keep.indexes,]
rownames(lps.ratios) <- probeset


rm(data.mus,data.mus.min0,data.fdr,data.ratios,data.lambdas)
