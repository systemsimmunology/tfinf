muLambda.data.directory   <- file.path(Sys.getenv("TFINF"),"expression_data")
muLambda.data.directory.2009   <- file.path(Sys.getenv("TFINF"),"expression_data")
## formerly
## muLambda.data.directory <- "/users/thorsson/macrophage/AffyArray/AllLambdaMu"
## muLambda.data.directory.2009 <- "/users/thorsson/macrophage/AffyArray/AllLambdaMuApril2009"

## LPS
data.mu.file <- read.table(paste(muLambda.data.directory,"LPS.alpha.mus",sep="/"),header=TRUE,sep='\t') ## March 2008: Not sure why this was here
lps.mus.2 <- apply(as.matrix(data.mu.file[3:8]),2,as.numeric)
rownames(lps.mus.2) <- as.character(data.mu.file$Probesets)
lps.mus.2 <- lps.mus.2[1:45037,]

mu.file <- paste(muLambda.data.directory,"LPS.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
lps.mus.3 <- apply(as.matrix(data.mu[3:13]),2,as.numeric)
rownames(lps.mus.3) <- as.character(data.mu$Probesets)
lps.mus.3 <- lps.mus.3[1:45037,]

### Poly IC 
mu.file <- paste(muLambda.data.directory,"PIC.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
polyIC.mus.2 <- apply(as.matrix(data.mu[3:8]),2,as.numeric)
rownames(polyIC.mus.2) <- as.character(data.mu$Probesets)
plot(polyIC.mus[1:45037,"min120"],polyIC.mus.2[,"min120"])

mu.file <- paste(muLambda.data.directory.2009,"PIC.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
polyIC.mus.3 <- apply(as.matrix(data.mu[3:11]),2,as.numeric)
rownames(polyIC.mus.3) <- as.character(data.mu$Probesets)

plot(polyIC.mus.2[1:45037,"min120"],polyIC.mus.3[,"min120"])

### R848
mu.file <- paste(muLambda.data.directory,"R848.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
r848.mus.2 <- apply(as.matrix(data.mu[3:8]),2,as.numeric)
rownames(r848.mus.2) <- as.character(data.mu$Probesets)
plot(r848.mus[1:45037,"min120"],r848.mus.2[,"min120"])

mu.file <- paste(muLambda.data.directory.2009,"R848.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
r848.mus.3 <- apply(as.matrix(data.mu[3:11]),2,as.numeric)
rownames(r848.mus.3) <- as.character(data.mu$Probesets)
plot(r848.mus.2[1:45037,"min120"],r848.mus.3[,"min120"])

### CpG
mu.file <- paste(muLambda.data.directory,"CPG.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
cpg.mus.2 <- apply(as.matrix(data.mu[3:8]),2,as.numeric)
rownames(cpg.mus.2) <- as.character(data.mu$Probesets)
plot(cpg.mus[1:45037,"min120"],cpg.mus.2[,"min120"])

### PAM2
mu.file <- paste(muLambda.data.directory,"PAM2.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
pam2.mus.2 <- apply(as.matrix(data.mu[3:8]),2,as.numeric)
rownames(pam2.mus.2) <- as.character(data.mu$Probesets)
plot(pam2.mus[1:45037,"min120"],pam2.mus.2[,"min120"])

### PAM3
mu.file <- paste(muLambda.data.directory,"PAM3.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
pam3.mus.2 <- apply(as.matrix(data.mu[3:8]),2,as.numeric)
rownames(pam3.mus.2) <- as.character(data.mu$Probesets)
plot(pam3.mus[1:45037,"min120"],pam3.mus.2[,"min120"])

mu.file <- paste(muLambda.data.directory.2009,"PAM3.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
pam3.mus.3 <- apply(as.matrix(data.mu[3:11]),2,as.numeric)
rownames(pam3.mus.3) <- as.character(data.mu$Probesets)

plot(pam3.mus.3[1:45037,"min120"],pam3.mus.2[,"min120"])

### MyD88 KO LPS
mu.file <- paste(muLambda.data.directory,"MyD88KO_LPS.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
myd88KOlps.mus.2 <- apply(as.matrix(data.mu[3:6]),2,as.numeric)
rownames(myd88KOlps.mus.2) <- as.character(data.mu$Probesets)
plot(myd88KOlps.mus[1:45037,"min120"],myd88KOlps.mus.2[,"min120"])

## The do not agree as well. The mu lamba files were obtained from a clone
## file I generated Aug31. The getExpression uploads were done earlier

### MyD88 KO PolyIC
mu.file <- paste(muLambda.data.directory,"MyD88KO_PIC.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
myd88KOpolyIC.mus.2 <- apply(as.matrix(data.mu[3:6]),2,as.numeric)
rownames(myd88KOpolyIC.mus.2) <- as.character(data.mu$Probesets)
plot(myd88KOpolyIC.mus[1:45037,"min120"],myd88KOpolyIC.mus.2[,"min120"])

### MyD88 KO PAM3
mu.file <- paste(muLambda.data.directory,"MyD88KO_PAM3.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
myd88KOpam3.mus.2 <- apply(as.matrix(data.mu[3:6]),2,as.numeric)
rownames(myd88KOpam3.mus.2) <- as.character(data.mu$Probesets)
plot(myd88KOpam3.mus[1:45037,"min120"],myd88KOpam3.mus.2[,"min120"])

### TRIF  KO LPS
mu.file <- paste(muLambda.data.directory,"TRIFKO_LPS.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
trifKOlps.mus.2 <- apply(as.matrix(data.mu[3:5]),2,as.numeric)
rownames(trifKOlps.mus.2) <- as.character(data.mu$Probesets)
plot(trifKOlps.mus[1:45037,"min120"],trifKOlps.mus.2[,"min120"])

mu.file <- paste(muLambda.data.directory,"TRIFKO_PAM2.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
trifKOpam2.mus <- apply(as.matrix(data.mu[3:5]),2,as.numeric)
rownames(trifKOpam2.mus) <- as.character(data.mu$Probesets)

### PAM3 PolyIC
mu.file <- paste(muLambda.data.directory,"PAM3PIC.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
PAM3polyIC.mus.2 <- apply(as.matrix(data.mu[3:6]),2,as.numeric)
rownames(PAM3polyIC.mus.2) <- as.character(data.mu$Probesets)
##plot(PAM3polyIC.mus[1:45037,"min120"],PAM3polyIC.mus.2[,"min120"])

### ATF3 KO
mu.file <- paste(muLambda.data.directory,"ATF3KO_LPS.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
atf3KOlps.mus <- apply(as.matrix(data.mu[3:6]),2,as.numeric)
rownames(atf3KOlps.mus) <- as.character(data.mu$Probesets)

mu.file <- paste(muLambda.data.directory,"ATF3KO_PAM2.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
atf3KOpam2.mus <- apply(as.matrix(data.mu[3:5]),2,as.numeric)
rownames(atf3KOpam2.mus) <- as.character(data.mu$Probesets)


mu.file <- paste(muLambda.data.directory,"ATF3KO_POLYIC.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
atf3KOpolyIC.mus <- apply(as.matrix(data.mu[3:5]),2,as.numeric)
rownames(atf3KOpolyIC.mus) <- as.character(data.mu$Probesets)


mu.file <- paste(muLambda.data.directory,"ATF3KO_CPG.mus",sep="/")
data.mu <- read.table(mu.file,header=TRUE,sep='\t')
atf3KOcpg.mus <- apply(as.matrix(data.mu[3:5]),2,as.numeric)
rownames(atf3KOcpg.mus) <- as.character(data.mu$Probesets)


######### We need the zero values ###########
lps.original.mus <- lps.mus ## before including hour 6, 8
lps.mus <- lps.mus.3 ## after including hour 6, 8
pam2.mus <- pam2.mus.2
pam3.mus <- pam3.mus.3 ## after including hour 4, 8, 12
polyIC.mus <- polyIC.mus.3 ## after including hour 4, 8, 12
r848.mus <- r848.mus.3 ## after including hour 4, 8, 12
cpg.mus <- cpg.mus.2


############ Cube, for wild-type

ps.max <- 45037

muCube <- rep(NA,6*ps.max*6)
dim(muCube) <- c(6,ps.max,6)
dimnames(muCube)[[3]] <- c("min0","min20","min40","min60","min80","min120")
dimnames(muCube)[[2]] <- probeset[1:ps.max]
dimnames(muCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

muCube["LPS",,] <- lps.mus[1:ps.max,1:6]
muCube["PAM2",,] <- pam2.mus[1:ps.max,1:6]
muCube["PAM3",,] <- pam3.mus[1:ps.max,1:6]
muCube["PolyIC",,] <- polyIC.mus[1:ps.max,1:6]
muCube["R848",,] <- r848.mus[1:ps.max,1:6]
muCube["CpG",,] <- cpg.mus[1:ps.max,1:6]

