muLambda.data.directory   <- file.path(Sys.getenv("HOME"),"tfinf/expression_data")
## Was
## muLambda.data.directory <- "/users/thorsson/macrophage/AffyArray/AllLambdaMu"

ratios.file <- paste(muLambda.data.directory,"LPS.log10_ratios",sep="/")
data.ratios <- read.table(ratios.file,header=TRUE,sep='\t')
lps.ratios.2 <- apply(as.matrix(data.ratios[3:12]),2,as.numeric)
rownames(lps.ratios.2) <- as.character(data.ratios$Probesets)

plot(lps.ratios[,"min120"],lps.ratios.2[,"min120"])

lps.original.ratios <- lps.ratios
lps.ratios <- lps.ratios.2
