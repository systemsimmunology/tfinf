
muLambda.data.directory   <- file.path(Sys.getenv("TFINF"),"expression_data")
muLambda.data.directory.2009   <- file.path(Sys.getenv("TFINF"),"expression_data")

## Formerly
##muLambda.data.directory <- "/users/thorsson/macrophage/AffyArray/AllLambdaMu"
##muLambda.data.directory.2009 <- "/users/thorsson/macrophage/AffyArray/AllLambdaMuApril2009"

## LPS 
##data.muStdErr.file <- read.table("LPS.muStandardErrors",header=TRUE,sep='\t')
##lps.muStdErrs <- apply(as.matrix(data.muStdErr.file[3:8]),2,as.numeric)
##rownames(lps.muStdErrs) <- as.character(data.muStdErr.file$Probesets)
##lps.muStdErrs <- lps.muStdErrs[1:45037,]

stderr.file <- paste(muLambda.data.directory,"LPS.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
lps.muStdErrs <- apply(as.matrix(data.muStdErr[3:13]),2,as.numeric)
rownames(lps.muStdErrs) <- as.character(data.muStdErr$Probesets)


### Poly IC 
stderr.file <- paste(muLambda.data.directory,"PIC.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
polyIC.muStdErrs <- apply(as.matrix(data.muStdErr[3:8]),2,as.numeric)
rownames(polyIC.muStdErrs) <- as.character(data.muStdErr$Probesets)

stderr.file <- paste(muLambda.data.directory.2009,"PIC.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
polyIC.muStdErrs.3 <- apply(as.matrix(data.muStdErr[3:11]),2,as.numeric)
rownames(polyIC.muStdErrs.3) <- as.character(data.muStdErr$Probesets)

### R 848
stderr.file <- paste(muLambda.data.directory.2009,"R848.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
r848.muStdErrs <- apply(as.matrix(data.muStdErr[3:8]),2,as.numeric)
rownames(r848.muStdErrs) <- as.character(data.muStdErr$Probesets)

### CpG
stderr.file <- paste(muLambda.data.directory,"CPG.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
cpg.muStdErrs <- apply(as.matrix(data.muStdErr[3:8]),2,as.numeric)
rownames(cpg.muStdErrs) <- as.character(data.muStdErr$Probesets)

### PAM2
stderr.file <- paste(muLambda.data.directory,"PAM2.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
pam2.muStdErrs <- apply(as.matrix(data.muStdErr[3:8]),2,as.numeric)
rownames(pam2.muStdErrs) <- as.character(data.muStdErr$Probesets)

### PAM3
stderr.file <- paste(muLambda.data.directory.2009,"PAM3.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
pam3.muStdErrs <- apply(as.matrix(data.muStdErr[3:8]),2,as.numeric)
rownames(pam3.muStdErrs) <- as.character(data.muStdErr$Probesets)

### MyD88 KO LPS
stderr.file <- paste(muLambda.data.directory,"MyD88KO_LPS.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
myd88KOlps.muStdErrs <- apply(as.matrix(data.muStdErr[3:5]),2,as.numeric)
rownames(myd88KOlps.muStdErrs) <- as.character(data.muStdErr$Probesets)

### MyD88 KO polyIC
stderr.file <- paste(muLambda.data.directory,"MyD88KO_PIC.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
myd88KOpolyIC.muStdErrs <- apply(as.matrix(data.muStdErr[3:5]),2,as.numeric)
rownames(myd88KOpolyIC.muStdErrs) <- as.character(data.muStdErr$Probesets)

### PAM3 PolyIC
stderr.file <- paste(muLambda.data.directory,"PAM3PIC.muStandardErrors",sep="/")
data.muStdErr <- read.table(stderr.file,header=TRUE,sep='\t')
PAM3polyIC.muStdErrs <- apply(as.matrix(data.muStdErr[3:5]),2,as.numeric)
rownames(PAM3polyIC.muStdErrs) <- as.character(data.muStdErr$Probesets)

pam3polyIC.muStdErrs <- PAM3polyIC.muStdErrs
rm(PAM3polyIC.muStdErrs)



############ Cube, for wild-type

ps.max <- 45037

muStdErrCube <- rep(NA,6*ps.max*6)
dim(muStdErrCube) <- c(6,ps.max,6)
dimnames(muStdErrCube)[[3]] <- c("min0","min20","min40","min60","min80","min120")
dimnames(muStdErrCube)[[2]] <- probeset[1:ps.max]
dimnames(muStdErrCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

muStdErrCube["LPS",,] <- lps.muStdErrs[1:ps.max,1:6]
muStdErrCube["PAM2",,] <- pam2.muStdErrs[1:ps.max,1:6]
muStdErrCube["PAM3",,] <- pam3.muStdErrs[1:ps.max,1:6]
muStdErrCube["PolyIC",,] <- polyIC.muStdErrs[1:ps.max,1:6]
muStdErrCube["R848",,] <- r848.muStdErrs[1:ps.max,1:6]
muStdErrCube["CpG",,] <- cpg.muStdErrs[1:ps.max,1:6]

## 
polyIC.muStdErrs <- polyIC.muStdErrs.3

