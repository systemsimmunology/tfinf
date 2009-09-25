

## Read and save all expression data 

source("./newFunctions.R")

source("./sbeamsReaderFunction.R")
source("./sbeamsReader.R") ## reads LPS
source("./sbeamsReaderWrapper.R") ## read other than LPS
source("./readMus.R")
source("./readLambdas.R")
source("./readStdErrs.R")
source("./readMyD88Ratios.R")
source("./readLPSratios.R")

expression.directory <- file.path(Sys.getenv("HOME"),"tfinf/expression_data")
annotation.directory <- file.path(Sys.getenv("HOME"),"tfinf/annotations")

###
affies.annotations <- read.table(paste(annotation.directory,"Mouse430_2.na25.annot_ncbi_mappings.tsv",sep="/"), header=TRUE, sep='\t', strip.white=TRUE,as.is=TRUE)
##affies.annotations <- read.table(paste(annotation.directory,"Mouse430_2.na29.annot_ncbi_mappings.tsv",sep="/"), header=TRUE, sep='\t', strip.white=TRUE,as.is=TRUE)
probeset <- affies.annotations[["Probe.Set.ID"]][1:45037]
cname.compare <- affies.annotations[["Gene.Symbol"]][1:45037]
names(cname.compare) <- probeset[1:45037]
ncbiID <- affies.annotations[["Entrez.Gene"]][1:45037]
names(ncbiID) <- probeset[1:45037]

##source("repairBadMaps.R")
##August 2009: Skip this, as for earlier cases, not needed or not helpful.

source("repairExpressionObjectNames.R") ## removes the .2 extension

source("getRepresentativeProbes.R")

source("countSigGenes.R")

## Save everything

##save.image(file="allExpression.RData")

## Extract relevant portions and save 
 
all.ratios.objects <- ls()[grep("ratios$",ls(),extended=TRUE)]
save ( list=all.ratios.objects, file=paste(expression.directory,"all.ratios.objects.RData",sep="/"))

all.mus.objects <- ls()[grep("mus$",ls(),extended=TRUE)]
save ( list=all.mus.objects, file=paste(expression.directory,"all.mus.objects.RData",sep="/"))

all.lambdas.objects <- ls()[grep("lambdas$",ls(),extended=TRUE)]
save ( list=all.lambdas.objects, file=paste(expression.directory,"all.lambdas.objects.RData",sep="/"))

all.muStdErrs.objects <- ls()[grep("muStdErrs$",ls(),extended=TRUE)]
save ( list=all.muStdErrs.objects, file=paste(expression.directory,"all.muStdErrs.objects.RData",sep="/"))

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

