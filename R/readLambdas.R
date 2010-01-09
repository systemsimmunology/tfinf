muLambda.data.directory   <- file.path(Sys.getenv("TFINF"),"expression_data")
muLambda.data.directory.2009   <- file.path(Sys.getenv("TFINF"),"expression_data")

## Formerly 
## muLambda.data.directory <- "/users/thorsson/macrophage/AffyArray/AllLambdaMu"
## muLambda.data.directory.2009 <- "/users/thorsson/macrophage/AffyArray/AllLambdaMuApril2009"

## LPS 
##data.lambda.file <- read.table("LPS.lambdas",header=TRUE,sep='\t')
##lps.original.lambdas <- lps.lambdas 
lambda.file <- paste(muLambda.data.directory,"LPS.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
lps.lambdas <- apply(as.matrix(data.lambda[3:12]),2,as.numeric)
rownames(lps.lambdas) <- as.character(data.lambda$Probesets)
lps.lambdas <- lps.lambdas[1:45037,]

### Poly IC 
lambda.file <- paste(muLambda.data.directory,"PIC.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
polyIC.lambdas <- apply(as.matrix(data.lambda[3:7]),2,as.numeric)
rownames(polyIC.lambdas) <- as.character(data.lambda$Probesets)

lambda.file <- paste(muLambda.data.directory.2009,"PIC.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
polyIC.lambdas.3 <- apply(as.matrix(data.lambda[3:10]),2,as.numeric)
rownames(polyIC.lambdas.3) <- as.character(data.lambda$Probesets)

### R 848
lambda.file <- paste(muLambda.data.directory.2009,"R848.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
r848.lambdas <- apply(as.matrix(data.lambda[3:7]),2,as.numeric)
rownames(r848.lambdas) <- as.character(data.lambda$Probesets)

### CpG
lambda.file <- paste(muLambda.data.directory,"CPG.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
cpg.lambdas <- apply(as.matrix(data.lambda[3:7]),2,as.numeric)
rownames(cpg.lambdas) <- as.character(data.lambda$Probesets)

### PAM2
lambda.file <- paste(muLambda.data.directory,"PAM2.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
pam2.lambdas <- apply(as.matrix(data.lambda[3:7]),2,as.numeric)
rownames(pam2.lambdas) <- as.character(data.lambda$Probesets)

### PAM3
lambda.file <- paste(muLambda.data.directory.2009,"PAM3.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
pam3.lambdas <- apply(as.matrix(data.lambda[3:7]),2,as.numeric)
rownames(pam3.lambdas) <- as.character(data.lambda$Probesets)

### MyD88 KO LPS
lambda.file <- paste(muLambda.data.directory,"MyD88KO_LPS.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
myd88KOlps.lambdas <- apply(as.matrix(data.lambda[3:5]),2,as.numeric)
rownames(myd88KOlps.lambdas) <- as.character(data.lambda$Probesets)

### MyD88 KO polyIC
lambda.file <- paste(muLambda.data.directory,"MyD88KO_PIC.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
myd88KOpolyIC.lambdas <- apply(as.matrix(data.lambda[3:5]),2,as.numeric)
rownames(myd88KOpolyIC.lambdas) <- as.character(data.lambda$Probesets)

### PAM3 PolyIC
lambda.file <- paste(muLambda.data.directory,"PAM3PIC.lambdas",sep="/")
data.lambda <- read.table(lambda.file,header=TRUE,sep='\t')
PAM3polyIC.lambdas <- apply(as.matrix(data.lambda[3:5]),2,as.numeric)
rownames(PAM3polyIC.lambdas) <- as.character(data.lambda$Probesets)

pam3polyIC.lambdas <- PAM3polyIC.lambdas
rm(PAM3polyIC.lambdas)



############ Cube, for wild-type

ps.max <- 45037

lambdaCube <- rep(NA,6*ps.max*5)
dim(lambdaCube) <- c(6,ps.max,5)
dimnames(lambdaCube)[[3]] <- c("min20","min40","min60","min80","min120")
dimnames(lambdaCube)[[2]] <- probeset[1:ps.max]
dimnames(lambdaCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

lambdaCube["LPS",,] <- lps.lambdas[1:ps.max,1:5]
lambdaCube["PAM2",,] <- pam2.lambdas[1:ps.max,1:5]
lambdaCube["PAM3",,] <- pam3.lambdas[1:ps.max,1:5]
lambdaCube["PolyIC",,] <- polyIC.lambdas[1:ps.max,1:5]
lambdaCube["R848",,] <- r848.lambdas[1:ps.max,1:5]
lambdaCube["CpG",,] <- cpg.lambdas[1:ps.max,1:5]


## 
polyIC.lambdas <- polyIC.lambdas.3
