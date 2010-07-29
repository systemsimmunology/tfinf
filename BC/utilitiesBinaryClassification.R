
### Input: a square distance matrix, named in both dimensions
###       assume ranging from "close"=low values to "far"=high values
###       
binaryClassificationMats <- function (psois, cloud, featurevec, distmat, pgsea=1e-8 ) {

  checkBinaryClassPars(psois, cloud, featurevec, distmat)
  
  n.psois <- length(psois)
  n.cloud <- length(cloud)

  logmat <- getLogMat(psois, cloud, featurevec, distmat)
  mmat <- t(apply(logmat,1,cumsum)) ## matrix of "m" values  
  matn <- t(matrix(rep((1:(n.cloud-1)),n.psois),ncol=n.psois))

  TPmat <- mmat
  FPmat <- matn - mmat
  FNmat <- M - mmat
  TNmat <- N - TPmat - FPmat - FNmat
  FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
  TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 
  Precmat <- TPmat/(TPmat+FPmat) ## Precision. Limits to M/N
  Recmat <- TPRmat
  
  rm(mmat,matn) ## cleanup, may not be needed

  ## PrecRec AUC
  xs <- Recmat
  ys <- Precmat
  del.xs <- xs[,2:(n.cloud-1)]-xs[,1:(n.cloud-2)]
  avg.ys <- (ys[,2:(n.cloud-1)]+ys[,1:(n.cloud-2)])/2.
  precrec.auc <- apply(del.xs*avg.ys,1,sum)
  precrec.auc.vs.n <- t(apply(del.xs*avg.ys,1,cumsum))
  
  ## ROC AUC   
  xs <- FPRmat
  ys <- TPRmat
  del.xs <- xs[2:(n.cloud-1)]-xs[1:(n.cloud-2)]
  avg.ys <- (ys[2:(n.cloud-1)]+ys[1:(n.cloud-2)])/2.
  del.xs <- xs[,2:(n.cloud-1)]-xs[,1:(n.cloud-2)]
  avg.ys <- (ys[,2:(n.cloud-1)]+ys[,1:(n.cloud-2)])/2.
  roc.auc <- apply(del.xs*avg.ys,1,sum)
  roc.auc.vs.n <- t(apply(del.xs*avg.ys,1,cumsum))

  ## Cleanup, not sure if it is needed
  rm(xs,ys,del.xs,avg.ys,roc.auc,precrec.auc)

  ## Distances go from 0 (close) to 2 (far)
  dist.mat.sorted <- t(apply(distmat[psois,],1,sort))
  dist.mat.sorted <- dist.mat.sorted[,2:n.cloud]
  
  returnList <- list(FPR=FPRmat,TPR=TPRmat,Prec=Precmat,Rec=Recmat,
                     TP=TPmat, FP=FPmat, FN=FNmat, TN=TNmat, 
                     precrec.auc.vs.n=precrec.auc.vs.n, roc.auc.vs.n=roc.auc.vs.n,
                     dist.mat.sorted=dist.mat.sorted)

  list.gsea <- gsea(psois=psois,cloud=cloud,featurevec=featurevec,distmat=distmat,pgsea=pgsea) 
  returnList <- c(returnList,list.gsea)
  
  return (returnList )
}


### Input: a square distance matrix, named in both dimensions
###       assume ranging from "close"=low values to "far"=high values
###

getDistMatSorted <- function ( distmat, psois, n.cloud ){
  ## Distances go from 0 (close) to 2 (far)
  dist.mat.sorted <- t(apply(distmat[psois,],1,sort))
  dist.mat.sorted <- dist.mat.sorted[,2:n.cloud]
  return(dist.mat.sorted)
}

## for set of genes, find neighbors in a cloud
## Neighbor has distance below a given distance threshold
getNbrs <- function ( psois, cloud, distmat, threshold ){
  outList <- list()
  for ( psoi in psois ){
    cloudnoself <- setdiff(cloud,psoi)
    keepers <- which( distmat[psoi,cloudnoself] <= threshold )
    if ( length(keepers) >  0 ) {
      keepers.ps <- cloudnoself[keepers]
      keepers.ps <- names(sort(distmat[psoi,keepers.ps]))
      outList[[psoi]] <- keepers.ps
    } else {
      outList[[psoi]] <- NULL
    }
  }    
  return(outList)
}

### This is a pgsea=0 version
binaryClassificationMatsLean <- function (psois, cloud, featurevec, distmat ) {
  
  n.psois <- length(psois)
  n.cloud <- length(cloud)

  ptm1 <- proc.time() 
  
  logmat <- getLogMat(psois, cloud, featurevec, distmat)

  ptm2 <- proc.time()
  cat("getlogmat",ptm2-ptm1,"\n")
  
  mmat <- t(apply(logmat,1,cumsum)) ## matrix of "m" values  
  
  ptm3 <- proc.time()
  cat("mmat",ptm3-ptm2,"\n")

  matn <- t(matrix(rep((1:(n.cloud-1)),n.psois),ncol=n.psois))

  ptm4 <- proc.time()
  cat("matn",ptm4-ptm3,"\n")
  
  ### TEMPORARY FOR TESTING
  ##mmat <- mmat[,1:100]
  ##matn <- matn[,1:100]
  
  TPmat <- mmat
  FPmat <- matn - mmat
  FNmat <- M - mmat
  TNmat <- N - TPmat - FPmat - FNmat
  FPRmat <- FPmat/(FPmat+TNmat) ## False Positive Rate
  TPRmat <- TPmat/(TPmat+FNmat) ## True Positive Rate = Recall = Sensitivity ! 

  ESmat <- TPRmat - FPRmat 
  
  ES <- ESmat

  ptm5 <- proc.time()
  cat("rest",ptm5-ptm4,"\n")

  return ( ES )
}





## ES(S), KS and all tha
## Input:
##    pgsea: GSEA parameter p (to get correct p=0 limit, use a tiny p )
##    abs.flag: set according to whether or not you want to use absolute value 
##
gsea <- function (psois, cloud, featurevec, distmat,pgsea=1e-8, abs.flag=FALSE) {

  checkBinaryClassPars(psois, cloud, featurevec, distmat)

  n.psois <- length(psois)
  n.cloud <- length(cloud)

  logmat <- getLogMat(psois, cloud, featurevec, distmat) ## this could maybe be passed in, instead

  ## Distances go from 0 (close) to 2 (far)
  dist.mat.sorted <- t(apply(distmat[psois,],1,sort))
  dist.mat.sorted <- dist.mat.sorted[,2:n.cloud]
  ## Correlations go from 1 (close ) to -1 (far)
  cor.mat <- 1 - dist.mat.sorted
  ## Correlations squished to 1 (close) to 0 (far)
  scaled.cormat <- (cor.mat+1)/2.

  ### Phit
  if ( !abs.flag ){
    scaled.cormat.trueonly <- logmat*scaled.cormat
    Phit.genecontrib <- scaled.cormat.trueonly^pgsea
  }
  else if ( abs.flag ){
    cormat.trueonly.abs <- abs(logmat*cor.mat)
    Phit.genecontrib <- cormat.trueonly.abs^pgsea
  } 
  Phit.genecontrib.cumsum <-  t(apply(Phit.genecontrib,1,cumsum))
  Phitmat <- Phit.genecontrib.cumsum/Phit.genecontrib.cumsum[,n.cloud-2]


  ## Pmiss
  if ( abs.flag ){
    cormat.falseonly.abs <- abs((!logmat)*cor.mat)
    Pmiss.genecontrib <- cormat.falseonly.abs^pgsea
  } else if ( !abs.flag ){
    scaled.cormat.falseonly  <- (!logmat)*scaled.cormat
    Pmiss.genecontrib <- scaled.cormat.falseonly^pgsea
  }
  Pmiss.genecontrib.cumsum <-  t(apply(Pmiss.genecontrib,1,cumsum))
  Pmissmat <- Pmiss.genecontrib.cumsum/Pmiss.genecontrib.cumsum[,n.cloud-2] 
  
  ## Subramanian et al GSEA statistic ES(S,i) prior to maximizatoin over i
  ESmat <- Phitmat-Pmissmat

  returnList <- list(pgsea=pgsea,ES=ESmat,Phit=Phitmat,Pmiss=Pmissmat,Phit.genecontrib=Phit.genecontrib,Pmiss.genecontrib=Pmiss.genecontrib)
  return (returnList )
  
}  


## For a given feature and genes, return a matrix with
## rows = genes
## row vector: order by neighbors, indicating wether neighbor has feature
getLogMat <- function(psois, cloud, featurevec, distmat){
  sortindex <- function(x){sort(x,index.return=TRUE)$ix} ## this could be placed "outside"
  ind.mat <- t(apply(distmat[psois,cloud],1,sortindex))
  ind.mat <- ind.mat[,2:ncol(ind.mat)]
  logmat <- matrix(featurevec[ind.mat],nrow=nrow(ind.mat) )
  rownames(logmat) <- psois
  return (logmat)
}




checkBinaryClassPars <- function (psois, cloud, featurevec, distmat) {
  if ( length(setdiff(psois,rownames(distmat))) != 0 ){
    stop("Error: Some psois are not in the distance matrix\n")
  }
  if ( length(setdiff(cloud,colnames(distmat))) != 0 ){
    stop("Error: Some cloud variables are not in the distance matrix\n")
  }
  if ( length(setdiff(cloud,names(has.feature.by.ps))) != 0 ){
    stop("Error: Some cloud variables are not in the feature vector\n")
  }
}  


##
## processFeatureVector - Process feature vector by Ensembl ID to one by Entrez ID
## 
## Input: feature vector by EnsemblID (logical vector)
## Output: feature vector by EntrezID (logical vector)
## Processing:
##   Map ensemblIDs to entrezIDs
##   Where the first mapping is multiple to one,
##   retain those for which the calls agree
##   toss the ones that the promoters "disagree" e.g. 1=found 2=not found

processEnsidFeatureVector <- function ( featurevec.ensids ) { 

  redundancy <- unlist(lapply(ensemblIDsOfEntrezIDs[expressed.scanned.egid],length))
  good.eids <- names(which(redundancy==1)) ## Gene IDs with non-redundant ensembl IDs
  good.ensids <- as.character(unlist(ensemblIDsOfEntrezIDs[good.eids])) 
  baddies <- names(which(redundancy>1)) ## Gene IDs with redundant ensembl IDs
  ## The above could be precomputed, outside this function
  
  ## First treat non-redundant set  
  has.feature <- featurevec.ensids[good.ensids]
  names(has.feature) <- good.eids

  ##salvageable if promoters agree
  salvageable <- character()
  has.to.append <- logical()
  for ( baddy in baddies ){
    logvec <- featurevec.ensids[ensemblIDsOfEntrezIDs[[baddy]]]
    crit1 <- (sum(logvec)==length(logvec)) ## all true
    crit2 <- (sum(logvec)==0) ## all false
    if ( crit1 | crit2 ) {
      salvageable <- c(salvageable, baddy )
      logical.value <- c(FALSE,TRUE)[sum(logvec)/length(logvec)+1]
      has.to.append <- c(has.to.append,logical.value)
    }
  }
  
  ## these had conflicts
  troublemaker.eids  <-  setdiff(baddies,salvageable)
  
  ## Append the ones that are consistent over mutliple ensids
  names(has.to.append) <- salvageable
  has.feature.by.eid <- c(has.feature,has.to.append)

  return (has.feature.by.eid) 

}
  
##
## Plot ROC curve 
## (uses matrix representations FPR and TPR )

plotROC <- function( psoi,FPRmat,TPRmat,minval=0,maxval=0.1,maxind=(N-1)){
  main <- as.character(cc[psoi])
  if ( FALSE ){
    score <- roc.auc.vs.n[psoi,maxind]/((FPRmat[psoi,maxind])^2/2)
    main <- paste("ROC:",as.character(cc[psoi]),"Score:",round(score,3))
  }
  if ( maxind == N){
    plot(FPRmat[psoi,],TPRmat[psoi,],xlim=c(minval,maxval),ylim=c(minval,maxval),xlab="False Positive Rate",ylab="True Positive Rate",main=main)
    abline(0,1)
  } else {
    plot(FPRmat[psoi,1:maxind],TPRmat[psoi,1:maxind],xlab="False Positive Rate",ylab="True Positive Rate",main=main)
    abline(0,1)
  }
}

##
## Plot Precison-Recall curve 
## (uses matrix representations Precmat and Recmat )

plotPrecRec <- function( psoi, Precmat, Recmat, prec.min=0., prec.max=0.01, maxind=(N-1), main=main){
  xminval <- prec.min
  xmaxval <- prec.max
  yminval <- 0.
  ymaxval <- 1.
  main <- as.character(cc[psoi])
  if ( FALSE ){
    score <- precrec.auc.vs.n[psoi,maxind]/Recmat[psoi,maxind]/(M/N)
    main <- paste("PrecRec:",as.character(cc[psoi]),"Score:",round(score,3))
  }
  if ( maxind == N ){
    plot(Recmat[psoi,],Precmat[psoi,],xlim=c(xminval,xmaxval),ylim=c(yminval,ymaxval),xlab="Recall",ylab="Precision",main=main)
  } else {
    plot(Recmat[psoi,1:maxind],Precmat[psoi,1:maxind],ylim=c(yminval,ymaxval),xlab="Recall",ylab="Precision",main=main)
  }
  abline(M/N,0)

}

## 
plotES <- function( psoi, pgsea,ESmat, maxind=N ){
  plotvec <- ESmat[psoi,1:maxind]
  maxval <- max(plotvec)
  mvs <- paste(round(maxval*1000,digits=1),"e-3",sep="")
  main <- paste("p=",pgsea,",ES(S)=",mvs,sep="")
  ##main <- paste("ES,p=",pgsea,",",as.character(cc[psoi]),",ES(S)=",mvs,sep="")
  plot(ESmat[psoi,1:maxind],xlab="i",ylab="ES(S,i)",type='l',main=main)
  abline(0,0)
}


## 
computeES <- function( psoi, pgsea,ESmat, maxind=N ){
  vec <- ESmat[psoi,1:maxind]
  maxval <- max(vec)
  return(maxval)
}


panelPlot <- function ( psoi, bco ,cordist, mincor=0.2 ) {

  if ( bco$dist.mat.sorted[psoi,1] < mincor ) {
    par(mfrow=c(2,3),oma=c(0,0,1,0))
    n.neighbors <- max(which(bco$dist.mat.sorted[psoi,]<mincor)) ## this used to be dist.mat
    plotROC(psoi,bco$FPR,bco$TPR,maxind=n.neighbors) ; ## roc.auc[psoi]-M/N
    plotPrecRec(psoi,bco$Prec,bco$Rec,maxind=n.neighbors) ; ##precrec.auc[psoi]-0.5
    plot(1:n.neighbors,llvec(psoi,bco,n.neighbors),type='l',main="Log Likelihood\n")
    tt <- names(sort(cordist[psoi,psois.full])[2:(n.neighbors+1)])
    ##tt <- c(psoi,names(corvals[[psoi]][1:(max(which(dist.mat[psoi,]<mincor)))]))
    profileplotRatio(lps.ratios[tt,],cc[tt],main="Expression Profiles")
    plotES(psoi,bco$pgsea,bco$ES,maxind=n.neighbors)
    main <- paste(cc[psoi],", Number of neighbors closer than ",mincor,":",n.neighbors)
    title(main, outer=TRUE,)
  }

}

## Log likelihood as function of neighbors
## can be slow !!
llvec <- function( psoi, bco, n ){
  llvals <- numeric()
  for ( i in 1:n ){
    test.mat <- rbind(c(bco$TP[psoi,i],bco$FP[psoi,i]),c(bco$FN[psoi,i],bco$TN[psoi,i]))
    llvals <- c(llvals,-log10(fisher.test(test.mat)$p.value))
  }
  llvals
}


## Estimate p-value of ESmax value
## Uses "cube" of ESmaxvalues derived from random:
## ESmaxcubeMgridNew
ESmaxPval <- function (esmaxval,em, nn,ESmaxcube ){
  mind <- nearestMpoint(em)
  Mval <- Msteps[mind]
  ESmaxdistrib <- ESmaxcube[mind,nn,]
  Fn <- ecdf(ESmaxdistrib)
  pval <- 1-Fn (esmaxval)
  return(pval)
}
nearestMpoint <- function ( Mval ){ ## helper function for the ESmaxPval. Needs generalization.
  index <- round(Mval/20)
  if ( index== 0 ){ index=1}
  return(index)
}

write.feature.vector <- function ( mm, filename ){
  matdim <- nrow(mm)
  linelength <- 80
  total.counter <- 0
  completed.line.counter <- 0
   
  for ( i in 1:(matdim-1)){
    for ( j in (i+1):matdim ){
      bitval <- mm[i,j]+0
      cat(as.character(bitval),file=filename,append=TRUE)
      total.counter <- total.counter + 1
      if ( total.counter/linelength==round(total.counter/linelength)){
        completed.line.counter <- completed.line.counter + 1
        cat( "\n",file=filename,append=TRUE)
      }
    }
  }
  cat("\n",file=filename,append=TRUE)
}

