## Initialize with initializeBC.R
##
load(paste(seq.dir,"mdof.RData",sep="/"))

expressedFamilies <- names(expressedFamilyMats)
n.fams <- length(expressedFamilies)
n.pairs <- n.fams*(n.fams-1)/2
 
fampairs <- c()
for ( i in 1:(n.fams-1)){
  for ( j in (i+1):n.fams ){
    streeng <- paste(expressedFamilies[i],expressedFamilies[j],sep="-")
    fampairs <- c(fampairs,streeng)
  }
}
 
featureMatrix <- array(FALSE,dim=c(n.pairs,length(eids.all)))
rownames(featureMatrix) <- fampairs
colnames(featureMatrix) <- eids.all

for ( i in 1:(n.fams-1)){
  for ( j in (i+1):n.fams ){
    fam1 <- expressedFamilies[i]
    fam2 <- expressedFamilies[j]
    streeng <- paste(fam1,fam2,sep="-")
    featureMatrix[streeng,] <- unlist(lapply(mdof[good.ensids],"[[",fam1,fam2))
  }
}

Mpair <- apply(featureMatrix,1,sum)


### Can we get the randomizations and scoring done fast?

### NN


t.index.min <- 1
t.index.max <- 11
cordistH <- 1-cor(t(lps.mus[psois.all,t.index.min:t.index.max]))

 
system.time( dist.mat.sorted <- getDistMatSorted(psois=psois.all,n.cloud=length(rownames(cordistH)), distmat=cordistH))

nns <- numeric(length=length(psois.all))
names(nns) <- psois.all

maxn <-  600 
mincor <- 0.1

for ( psoi in psois.all){
  if ( dist.mat.sorted[psoi,1] < mincor ) {
    n.neighbors <- max(which(dist.mat.sorted[psoi,]<mincor))
    if ( n.neighbors <= maxn ){
      nns[psoi] <- n.neighbors
    } else {
      nns[psoi] <- maxn
    }
  } else {
    nns[psoi] <- 0
  }
}



##### Cytokines
cytokines.es.psois <- read.table(file="cytokines.es.psois",as.is=TRUE)$V1
cytokinebinding.es.psois <- read.table(file="cytokinebinding.es.psois",as.is=TRUE)$V1
car.psois <- c(cytokines.es.psois,cytokinebinding.es.psois) 
psois <- car.psois

## Remember to exclude "the gene itself"

n.reps <- 10

### Can we create a random feature vector for
## a given gene, given it's nn
 
psoi <- car.psois[1]
nn <- nns[psoi]
eid <- ncbiID[psoi]

##eids.univ <- setdiff(eids.all,eid)
eids.univ <- eids.all ## include one more for generality of the randomization

N <- length(eids.univ)


### Single draw
frand <- featureMatrix[,sample(eids.univ,nn)]
hmm <- apply(frand,2,cumsum)


mmat <- t(hmm[1:10,1:3])

matn <- t(matrix(rep(1:10,3),ncol=3))

TPmat <- mmat
FPmat <- matn - mmat

## this seems to work, who knows why
FNmat <- Mpair[1:3] -mmat

##FNmat <- M - mmat
TNmat <- N - TPmat - FPmat - FNmat
FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 

ESmat <- TPRmat - FPRmat 

ES <- ESmat

apply(ES,1,max)
## gives ESmax per matrix for ths single draw.
## Can we do multiple draws?

## Eg 4 draws
## 3 site pairs
## nn=10

ptm <- proc.time()

##eids.univ <- setdiff(eids.all,eid)
eids.univ <- eids.all ## include one more for generality of the randomization

N <- length(eids.univ)

nn <- 200
npairs <- n.pairs
ndraws <- 100

matn <- t(matrix(rep(1:nn,npairs),ncol=npairs))
mcube <- array(0,dim=c(ndraws,npairs,nn))
ncube <- array(0,dim=c(ndraws,npairs,nn))
bigMcube <- array(0,dim=c(ndraws,npairs,nn))

for ( i in 1:ndraws ){
  cloudsamp <- sample(eids.univ,nn)
  frand <- featureMatrix[1:npairs,cloudsamp]
  hmm <- t(apply(frand,1,cumsum))
  ##  frand <- featureMatrix[,sample(eids.univ,nn)]
  mcube[i,,] <- hmm
  ncube[i,,] <- matn
  bigMcube[i,,] <- matrix(rep(Mpair[1:npairs],nn),ncol=nn)
}
TPcube <- mcube
FPcube <- ncube - mcube
FNcube <- bigMcube - mcube
rm(mcube,ncube,bigMcube)

##FNcube <- M - mcube
TNcube <- N - TPcube - FPcube - FNcube
FPRcube <- FPcube/(FPcube+TNcube) ## False Positive Rate
TPRcube <- TPcube/(TPcube+FNcube) ## True Positive Rate = Recall = Sensitivity ! 
EScube <- TPRcube - FPRcube 
rm(TPcube,FPcube,FNcube,TNcube,TPRcube,FPRcube)

### This has the right dims nsamps, npair 
proc.time() - ptm


ESmaxcube <- apply(EScube,c(1,2),cummax)


## This works (empirically!)
##ESmaxmat <- apply(EScube,c(1,2),max)


proc.time() - ptm

##
## Back to the original calculation: ESmax, from the *ordered set* per gene !!
##

## For every row of feature matrix, need ESmax for found features

## 

## Get neighbor set


## For each found feature
### Calculate ESvector
## Maximize 

## More matrixy:
### Calculate ESvalues for all neighbors and matrices
### Maximize over neighbor direction

maxn <-  600 
mincor <- 0.1
system.time( nbrs <- getNbrs(psois=psois.all,cloud=rownames(cordistH), distmat=cordistH, threshold=mincor))

ptm <- proc.time()
 
ffeature.list <- list()
##for ( psoi in psois.all  ){
for ( psoi in rc["Il12b"] ){
  nbs <- nbrs[[psoi]]
  nn <- length(nbs)

  if ( nn>5 ){

  nbs.eid <- as.character(ncbiID[nbs])
  fm.sub <- featureMatrix[,nbs.eid]
  mmat <- t(apply(fm.sub,1,cumsum))
  matn <- t(matrix(rep(1:nn,n.pairs),ncol=n.pairs))
  TPmat <- mmat
  FPmat <- matn - mmat
  FNmat <- Mpair- mmat
  TNmat <- N - TPmat - FPmat - FNmat
  FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
  TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
  ESmat <- TPRmat - FPRmat 
  ES <- ESmat
  Esmax <- apply(ES,1,max)

  randmat <- ESmaxcube[nn,,]

  pair.keep <- character()
  score.keep <- numeric()
  pval.keep <- numeric()
  
  for ( i in 1:n.pairs ){
    pval <- getpval(Esmax[i],allMs[i],nn)
    ##Fn <- ecdf(randmat[,i])
    ##pval <- 1-Fn(Esmax[i])
    if ( pval <= 0.05 ){
      pair.keep <- c(pair.keep,fampairs[i])
      score.keep <- c(score.keep,Esmax[i])
      names(pval) <- fampairs[i]
      pval.keep <- c(pval.keep,pval)
    }
  }
  feature.save <- rbind(score.keep,pval.keep)
  ffeature.list[[psoi]] <- feature.save 

}
  
} ## end loop over psois

##save(ffeature.list,file="ffeature.list.RData")

proc.time() - ptm



