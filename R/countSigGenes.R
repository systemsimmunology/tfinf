
source("./utilitiesExpression.R")

mu.cutoff <- 100
##mu.cutoff <- 0
lambda.cutoff <- 57.2
#lambda.cutoff <- 150

## May 2010
lambda.cutoff <- 26.61275

imax <- 6

stimulus <- "LPS"
lps.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,lps.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
lps.ps.sig <- setdiff(lps.ps.sig,low.expressors)
lps.cname.sig <- sort(unique(cname.compare[lps.ps.sig]))

stimulus <- "PAM2"
pam2.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,pam2.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
pam2.ps.sig <- setdiff(pam2.ps.sig,low.expressors)
pam2.cname.sig <- sort(unique(cname.compare[pam2.ps.sig]))

stimulus <- "PAM3"
pam3.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,pam3.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
pam3.ps.sig <- setdiff(pam3.ps.sig,low.expressors)
pam3.cname.sig <- sort(unique(cname.compare[pam3.ps.sig]))

stimulus <- "PolyIC"
polyIC.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,polyIC.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
polyIC.ps.sig <- setdiff(polyIC.ps.sig,low.expressors)
polyIC.cname.sig <- sort(unique(cname.compare[polyIC.ps.sig]))

stimulus <- "R848"
r848.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,r848.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
r848.ps.sig <- setdiff(r848.ps.sig,low.expressors)
r848.cname.sig <- sort(unique(cname.compare[r848.ps.sig]))

stimulus <- "CpG"
cpg.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,cpg.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
cpg.ps.sig <- setdiff(cpg.ps.sig,low.expressors)
cpg.cname.sig <- sort(unique(cname.compare[cpg.ps.sig]))

all.ps.sig <- union(lps.ps.sig,union(pam2.ps.sig,union(pam3.ps.sig,union(polyIC.ps.sig,union(r848.ps.sig,cpg.ps.sig)))))
all.cname.sig <- sort(unique(cname.compare[all.ps.sig]))
n.all.ps.sig <- length(all.ps.sig)

## Subtract one for the dash dash dash gene name
cat("No. LPS genes: ",length(lps.cname.sig)-1,"\n")
cat("No. PAM2 genes: ",length(pam2.cname.sig)-1,"\n")
cat("No. PAM3 genes: ",length(pam3.cname.sig)-1,"\n")
cat("No. PolyIC genes: ",length(polyIC.cname.sig)-1,"\n")
cat("No. R848 genes: ",length(r848.cname.sig)-1,"\n")
cat("No. CpG genes: ",length(cpg.cname.sig)-1,"\n")
cat("No. All genes: ",length(all.cname.sig)-1,"\n")


#compareSets(lps.cname.sig,pam2.cname.sig)


sigCompCube <- rep(NA,6*6*6)
dim(sigCompCube) <- c(6,6,6)
dimnames(sigCompCube)[[3]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
dimnames(sigCompCube)[[2]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
dimnames(sigCompCube)[[1]] <- c("a","b","a^b","a-b","b-a","a+b")

sigCompCube[,"LPS","PAM2"] <- compareSets(lps.cname.sig,pam2.cname.sig)
sigCompCube[,"LPS","PAM3"] <- compareSets(lps.cname.sig,pam3.cname.sig)
sigCompCube[,"LPS","PolyIC"] <- compareSets(lps.cname.sig,polyIC.cname.sig)
sigCompCube[,"LPS","R848"] <- compareSets(lps.cname.sig,r848.cname.sig)
sigCompCube[,"LPS","CpG"] <- compareSets(lps.cname.sig,cpg.cname.sig)

sigCompCube[,"PAM2","PAM3"] <- compareSets(pam2.cname.sig,pam3.cname.sig)
sigCompCube[,"PAM2","PolyIC"] <- compareSets(pam2.cname.sig,polyIC.cname.sig)
sigCompCube[,"PAM2","R848"] <- compareSets(pam2.cname.sig,r848.cname.sig)
sigCompCube[,"PAM2","CpG"] <- compareSets(pam2.cname.sig,cpg.cname.sig)

sigCompCube[,"PAM3","PolyIC"] <- compareSets(pam3.cname.sig,polyIC.cname.sig)
sigCompCube[,"PAM3","R848"] <- compareSets(pam3.cname.sig,r848.cname.sig)
sigCompCube[,"PAM3","CpG"] <- compareSets(pam3.cname.sig,cpg.cname.sig)

sigCompCube[,"PolyIC","R848"] <- compareSets(polyIC.cname.sig,r848.cname.sig)
sigCompCube[,"PolyIC","CpG"] <- compareSets(polyIC.cname.sig,cpg.cname.sig)

sigCompCube[,"R848","CpG"] <- compareSets(r848.cname.sig,cpg.cname.sig)

sigCompCube["a^b",,]/sigCompCube["a+b",,]


outMat <- matrix(NA,nrow=6,ncol=6)
rownames(outMat) <-  c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
colnames(outMat) <-  c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
m1 <- sigCompCube["a^b",,]/sigCompCube["a",,]
m2 <- t(sigCompCube["a^b",,]/sigCompCube["b",,])
outMat[upper.tri(outMat)] <- m1[upper.tri(m1)]
outMat[lower.tri(outMat)] <- m2[lower.tri(m2)]

## write.table(100*outMat,file="sigCompCube.tsv",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)



mu.cutoff <- 100
##mu.cutoff <- 0
lambda.cutoff <- 57.2
#lambda.cutoff <- 150


imax <- 11
lps.full.ps.sig <- rownames(sigSlice(lambda.cutoff,lps.mus[,1:imax],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[,1:imax]<mu.cutoff,1,sum)==imax))
lps.full.ps.sig <- setdiff(lps.full.ps.sig,low.expressors)

imax <- 6
lps.ps.sig <- rownames(sigSlice(lambda.cutoff,lps.mus[,1:imax],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[,1:imax]<mu.cutoff,1,sum)==imax))
lps.ps.sig <- setdiff(lps.ps.sig,low.expressors)

lps.full.ps.sig.rep <- sort(repProbes.cname[unique(sort(cname.compare[lps.full.ps.sig]))])

outMat1 <- cbind(lps.full.ps.sig.rep,cname.compare[lps.full.ps.sig.rep],ncbiID[lps.full.ps.sig.rep],lps.ratios[lps.full.ps.sig.rep,])
colnames(outMat1) <- c("Probeset","GeneName","GeneID",colnames(lps.ratios))
##write.table(outMat1,file="LPS.repProbe.ratios",quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

outMat1 <- cbind(lps.full.ps.sig.rep,cname.compare[lps.full.ps.sig.rep],ncbiID[lps.full.ps.sig.rep],lps.mus[lps.full.ps.sig.rep,])
colnames(outMat1) <- c("Probeset","GeneName","GeneID",colnames(lps.mus))
##write.table(outMat1,file="LPS.repProbe.mus",quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

outMat1 <- cbind(lps.full.ps.sig.rep,cname.compare[lps.full.ps.sig.rep],ncbiID[lps.full.ps.sig.rep],lps.lambdas[lps.full.ps.sig.rep,])
colnames(outMat1) <- c("Probeset","GeneName","GeneID",colnames(lps.lambdas))
##write.table(outMat1,file="LPS.repProbe.lambdas",quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

outMat1 <- cbind(lps.full.ps.sig.rep,cname.compare[lps.full.ps.sig.rep],ncbiID[lps.full.ps.sig.rep],lps.muStdErrs[lps.full.ps.sig.rep,])
colnames(outMat1) <- c("Probeset","GeneName","GeneID",colnames(lps.muStdErrs))
##write.table(outMat1,file="LPS.repProbe.muStdErrs",quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

##########################################################
# Overall signal integral
##########################################################

imax <- 6
m <- lps.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
lps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#lps.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1


imax <- 6
m <- pam2.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
pam2.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#pam2.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1


imax <- 6
m <- pam3.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
pam3.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#pam3.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- polyIC.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
polyIC.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#polyIC.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- r848.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
r848.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#r848.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- cpg.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
cpg.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#CpG.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

all.mic <- cbind(lps.mic,pam2.mic,pam3.mic,polyIC.mic,r848.mic,cpg.mic)
colnames(all.mic) <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
