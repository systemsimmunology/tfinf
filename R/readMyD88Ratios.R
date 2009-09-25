
muLambda.data.directory   <- file.path(Sys.getenv("HOME"),"tfinf/expression_data")
## was
## muLambda.data.directory <- "/users/thorsson/macrophage/AffyArray/AllLambdaMu"


### MyD88 KO LPS
ratios.file <- paste(muLambda.data.directory,"MyD88KO_LPS.log10_ratios",sep="/")
data.ratios <- read.table(ratios.file,header=TRUE,sep='\t')
myd88KOlps.ratios.2 <- apply(as.matrix(data.ratios[3:5]),2,as.numeric)
rownames(myd88KOlps.ratios.2) <- as.character(data.ratios$Probesets)

### MyD88 KO polyIC
ratios.file <- paste(muLambda.data.directory,"MyD88KO_PIC.log10_ratios",sep="/")
data.ratios <- read.table(ratios.file,header=TRUE,sep='\t')
myd88KOpolyIC.ratios.2 <- apply(as.matrix(data.ratios[3:5]),2,as.numeric)
rownames(myd88KOpolyIC.ratios.2) <- as.character(data.ratios$Probesets)

