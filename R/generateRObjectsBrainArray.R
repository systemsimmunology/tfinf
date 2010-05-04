 
source("./utilitiesExpression.R")
expression.directory <- file.path(Sys.getenv("AA"),"data/20100426.curated.3prime")

conds <- c("min0","min20","min40","min60","min80","min120")

## handy utility function - depends on globals, so don't move
miniRead <- function(filelabel,internallabel,conds){
  res <- readVERAandSAMfiles("PAM2",expression.directory,niceheaders=conds)
  assign(paste(c(internallabel,".mus"),collapse=""),res$mus)
  assign(paste(c(internallabel,".ratios"),collapse=""),res$ratios)
  assign(paste(c(internallabel,".lambdas"),collapse=""),res$lambdas)
  assign(paste(c(internallabel,".muStdErrs"),collapse=""),res$muStdErrs)
## ugh can't pass these out without creating another structure
}  

conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
res <- readVERAandSAMfiles("LPS",expression.directory,niceheaders=conds)
lps.mus <- res$mus
lps.ratios <- res$ratios
lps.lambdas <- res$ratios
lps.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
res <- readVERAandSAMfiles("PAM2",expression.directory,niceheaders=conds)
pam2.mus <- res$mus
pam2.ratios <- res$ratios
pam2.lambdas <- res$ratios
pam2.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
#original (male - female mix)
#conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("PAM3",expression.directory,niceheaders=conds)
pam3.mus <- res$mus
pam3.ratios <- res$ratios
pam3.lambdas <- res$ratios
pam3.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
#original (male - female mix)
#conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("R848",expression.directory,niceheaders=conds)
r848.mus <- res$mus
r848.ratios <- res$ratios
r848.lambdas <- res$ratios
r848.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
res <- readVERAandSAMfiles("CpG",expression.directory,niceheaders=conds)
cpg.mus <- res$mus
cpg.ratios <- res$ratios
cpg.lambdas <- res$ratios
cpg.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("Poly-IC",expression.directory,niceheaders=conds)
polyIC.mus <- res$mus
polyIC.ratios <- res$ratios
polyIC.lambdas <- res$ratios
polyIC.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
res <- readVERAandSAMfiles("MyD88-LPS",expression.directory,niceheaders=conds)
myd88KOlps.mus <- res$mus
myd88KOlps.ratios <- res$ratios
myd88KOlps.lambdas <- res$ratios
myd88KOlps.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
res <- readVERAandSAMfiles("MyD88-PolyIC",expression.directory,niceheaders=conds)
myd88KOpolyIC.mus <- res$mus
myd88KOpolyIC.ratios <- res$ratios
myd88KOpolyIC.lambdas <- res$ratios
myd88KOpolyIC.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
myd88KOpam3.mus <- read.matrix(paste(c(expression.directory,"/","MyD88-PAM3.mus"),collapse=""))
colnames(myd88KOpam3.mus) <- conds

conds <- c("min0","min60","min120")
trifKOpam2.mus <- read.matrix(paste(c(expression.directory,"/","TRIF-PAM2.mus"),collapse=""))
colnames(trifKOpam2.mus) <- conds

conds <- c("min0","min60","min120")
trifKOlps.mus <- read.matrix(paste(c(expression.directory,"/","TRIF-LPS.mus"),collapse=""))
colnames(trifKOlps.mus) <- conds

all.ratios.objects <- ls()[grep("ratios$",ls(),extended=TRUE)]
save ( list=all.ratios.objects, file=paste(expression.directory,"all.ratios.objects.RData",sep="/"))
all.mus.objects <- ls()[grep("mus$",ls(),extended=TRUE)]
save ( list=all.mus.objects, file=paste(expression.directory,"all.mus.objects.RData",sep="/"))
all.lambdas.objects <- ls()[grep("lambdas$",ls(),extended=TRUE)]
save ( list=all.lambdas.objects, file=paste(expression.directory,"all.lambdas.objects.RData",sep="/"))
all.muStdErrs.objects <- ls()[grep("muStdErrs$",ls(),extended=TRUE)]
save ( list=all.muStdErrs.objects, file=paste(expression.directory,"all.muStdErrs.objects.RData",sep="/"))

annotation.directory <- file.path(Sys.getenv("AA"),"annotation")

rt <- read.table(file.path(Sys.getenv("AA"),"annotation/MoGene4302_Mm_ENTREZG_mapfile"),sep="\t",as.is=TRUE,header=TRUE)
probeset <- rt$ProbesetID
ncbiID <- as.character(rt$EntrezID)
names(ncbiID) <- probeset
cname.compare <- rt$GeneSymbol
names(cname.compare) <- probeset
## aliases
gene.names.all <- cname.compare
commonName <- cname.compare
gene.ids.all <- probeset

load("~/data/ncbi/np.RData")
names(np) <- paste(names(np),"_at",sep="")

annotation.objects <- c("probeset",
                        "gene.ids.all",
                        "commonName",
                        "np",
                        "ncbiID",
                        "gene.names.all",
                        "cname.compare")

# identical: commonName, gene.names.all
# identical: probeset, gene.ids.all
# keep both, as functions may rely on them

save (list=annotation.objects, file=paste(annotation.directory,"annotation.objects.RData",sep="/"))

repProbes.cname <- probeset
names(repProbes.cname) <- cname.compare
repProbes.ncbiID <- probeset
names(repProbes.ncbiID) <- ncbiID
repProbes.np <- names(np)
names(repProbes.np) <- np

save(list=c("repProbes.cname","repProbes.np","repProbes.ncbiID"),file=paste(annotation.directory,"representativeProbes.RData",sep="/"))

source("countSigGenes.R")

## Save everythingq

##save.image(file="allExpression.RData")

## Extract relevant portions and save 

all.Cube.objects <- ls()[grep("Cube$",ls(),extended=TRUE)]
save( list=all.Cube.objects, file=paste(expression.directory,"all.Cube.objects.RData",sep="/"))

save(list=c("repProbes.cname","repProbes.np","repProbes.ncbiID"),file=paste(annotation.directory,"representativeProbes.RData",sep="/"))

##save(list=c("representativeProbes","repProbes.cname","repProbes.np","repProbes.ncbiID","repProbes.lfps.np","repProbes.lfps.cname","repProbes.lfps.ncbiID"),file="representativeProbes.RData")

## Older verson of GO analysis. Pull it out.
##source("compileGoListAll.R")
##save(list=c("goListAll.bp5","goListAll.mf5"),file="goListAll.RData")

save(all.mic,file=paste(expression.directory,"all.mic.RData",sep="/"))

all.ps.lists <- ls()[grep("ps.sig$",ls(),extended=TRUE)]
all.ps.lists <- c(all.ps.lists,"lps.full.ps.sig.rep")
save( list=all.ps.lists, file=paste(annotation.directory,"all.ps.list.objects.RData",sep="/"))

annotation.objects <- c("probeset",
                        "gene.ids.all",
                        "commonName",
                        "np",
                        "ncbiID",
                        "gene.names.all",
                        "cname.compare")

# identical: commonName, gene.names.all
# identical: probeset, gene.ids.all
# keep both, as functions may rely on them

save (list=annotation.objects, file=paste(annotation.directory,"annotation.objects.RData",sep="/"))

