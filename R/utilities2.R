
us <- function( vec ){ unique(sort(vec)) }

s2o <- function( string ) { eval(parse(text=string)) }

quickGroupPlot <- function ( psois, main="" ){
  profileplotRatio(lps.ratios[psois,],cname.compare[psois],main=main)
##  profileplotIntensity(lps.ratios[psois,],cname.compare[psois],main=main)
}
quickGP <- quickGroupPlot
qgp <- quickGroupPlot

## Do you have a variable containing inString?
gls <- function( inString ){
  return( ls(.GlobalEnv)[grep(inString,ls(.GlobalEnv))])
}


"lsfunc" <- function ( resp, time, krat, kd,deltaT ){
  pred <- krat * (1-exp(-kd*(time-deltaT)))
  (resp-pred)
}

#krat.external is an externally imposed krat
"lsfunc2" <- function ( resp, time, kd,deltaT ){
  pred <- krat.external * (1-exp(-kd*(time-deltaT)))
  (resp-pred)
}


# ls as a function of parameters only
# data.vals and t.samples need to be known externally
"lsfuncPars" <- function ( parvec ){
  resp <- data.vals
  time <- t.samples
  krat <- parvec[1]
  kd <- parvec[2]
  deltaT <- parvec[3]
  pred <- krat * (1-exp(-kd*(time-deltaT)))
  sum((resp-pred)^2)
}


# ls as a function of parameters only
# data.vals and t.samples need to be known externally
"lsfuncPars2" <- function ( parvec ){
  resp <- data.vals
  time <- t.samples
  krat <- krat.external
  kd <- parvec[1]
  deltaT <- parvec[2]
  pred <- krat * (1-exp(-kd*(time-deltaT)))
  sum((resp-pred)^2)
}


# without time shift
#"lsfunc" <- function ( resp, time, krat, kd ){
#  pred <- krat * (1-exp(-kd*time))
#  (resp-pred)/sqrt(resp)
#}


# pure exponential
"lsfuncexp" <- function ( resp, time, krat, kd ){
  pred <- krat * exp(kd*time)
  (resp-pred)/sqrt(resp)
}



"peakIndex" <- function( jj ) {
  n.j <- length(jj)
  j.change <- jj[2:n.j]-0.9*jj[1:(n.j-1)] ## allow some slop
  if ( sum(j.change<0) ){
    return.index <- min(which(j.change < 0))
  } else {
    return.index <- n.j
  }
  return.index
}

"plotFit" <- function(observed,predicted,main.title="",time.constant=999999,deltaT=0){

  times <- c(0,20,40,60,80,120) ## fixed for now, could use as an argument

  # ymax should bound the data above by integer * power of ten
  yy <- max(c(data.vals,pred))
  dec <- ceiling(log10(yy))-1
  ymax <- 10^dec*ceiling(yy/10^dec)
  
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  plot(times,observed, ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity Gain',axes=TRUE,main=main.title,lwd=3)
  lines(times,predicted, col='green',lwd=3)

  lines(c(deltaT,deltaT),c(0,ymax),col='magenta',lwd=3) 
  
  if ( time.constant != 999999 ){
    lines(c(time.constant+deltaT,time.constant+deltaT),c(0,ymax),col='red',lwd=3)
    legend(80,ymax/4., legend=c("data","prediction","start time","start+induction times"), col=c('blue','green','magenta','red'), lwd=3)
  } else {
    legend(80,ymax/4., legend=c("data","prediction","start time"), col=c('blue','green','magenta'), lwd=3)
  }

  
  par(op)

}


  



## returns pairs
"grabPairs" <- function(soi.tfs,near.tfs){

  ns <- length(soi.tfs)
  nt <- length(near.tfs)
  
  if ( ns ==0 || nt == 0){
    cat("Error: Missing TFs as input to grabPairs function\n")
  }

  if ( (FALSE %in% (soi.tfs %in% colnames(tf.dist.cn))) || (FALSE %in% (near.tfs %in% colnames(tf.dist.cn)))){
      cat("Error: No distance available ( in grabPairs function )\n")
      cat(setdiff(soi.tfs,colnames(tf.dist.cn)),"\n")
      cat(setdiff(near.tfs,colnames(tf.dist.cn)),"\n")
  }

  distMat <- tf.dist.cn[soi.tfs,near.tfs] ## Can be a vector in some cases
  collection <- character()
  ## Unfortunately, matrices collapsing to vectors (ns=1 and/or nt=1)
  ## have posed the usual R headaches with losing names
  ## and nrow/ncol evaluating to NaN and so on
  ## ( using data.frame didn't solve this completely )
  ## Hence all the case statements and cumbersome code
  
  if ( ns==1 && nt==1 ){
    partner1 <- soi.tfs
    partner2 <- near.tfs
    if ( distMat==1 | distMat==2 ){
      collection <- rbind(character(),c(partner1,partner2))
    }
  }## end ns==1, nt==1 case

  if ( ns>1 && nt==1){
    partner1 <- near.tfs
    one.indices <- which(distMat==1)
    two.indices <- which(distMat==2)
    keepers <- c(one.indices,two.indices)
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(near.tfs,length(keepers)),soi.tfs[keepers]) 
      rownames(collection) <- NULL
    }
    #if ( length(one.indices) != 0 ){
    #  for ( one.index  in 1:length(one.indices) ){
    #    partner2 <- soi.tfs[one.index]
    #    collection <- rbind(collection,c(partner1,partner2))
    #  }
    #}
    #if ( length(two.indices) != 0 ){
    #  for ( two.index  in 1:length(two.indices) ){
    #    partner2 <- soi.tfs[two.index]
    #    collection <- rbind(collection,c(partner1,partner2))
    #  }
    #}
  }## end ns>1, nt==1 case

  if ( nt>1 && ns==1){
    partner1 <- soi.tfs
    one.indices <- which(distMat==1)
    two.indices <- which(distMat==2)
    keepers <- c(one.indices,two.indices)
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(soi.tfs,length(keepers)),near.tfs[keepers]) 
      rownames(collection) <- NULL
    }
    ##if ( length(one.indices) != 0 ){
    ##  for ( one.index  in 1:length(one.indices) ){
    ##    partner2 <- near.tfs[one.index]
    ##    collection <- rbind(collection,c(partner1,partner2))
    ##  }
    ##}
    ##if ( length(two.indices) != 0 ){
    ##  for ( two.index  in 1:length(two.indices) ){
    ##    partner2 <- near.tfs[two.index]
    ##    collection <- rbind(collection,c(partner1,partner2))
    ##  }
    ##}
  }## end nt>1, ns==1 case

  if ( nt>1 &&  ns> 1 ){
    one.indices <- which(distMat==1,arr.ind=TRUE) ## 
    two.indices <- which(distMat==2,arr.ind=TRUE) 

    if ( length(one.indices) != 0 ){
      for ( one.index.row  in 1:nrow(one.indices) ){
        partner1 <- soi.tfs[one.indices[one.index.row,1]]
        partner2 <- near.tfs[one.indices[one.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
  
    if ( length(two.indices) != 0 ){
      for ( two.index.row  in 1:nrow(two.indices) ){
        partner1 <- soi.tfs[two.indices[two.index.row,1]]
        partner2 <- near.tfs[two.indices[two.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
    
  }## end nt>1, ns>1 case

  colnames(collection) <- NULL
  return (collection) 
}
## returns pairs


## Grab pairs, nearest neighbors only
"grabPairsNN" <- function(soi.tfs,near.tfs){

  ns <- length(soi.tfs)
  nt <- length(near.tfs)
  
  if ( ns ==0 || nt == 0){
    cat("Error: Missing TFs as input to grabPairs function\n")
  }

  if ( (FALSE %in% (soi.tfs %in% colnames(tf.dist.cn))) || (FALSE %in% (near.tfs %in% colnames(tf.dist.cn)))){
      cat("Error: No distance available ( in grabPairs function )\n")
  }

  distMat <- tf.dist.cn[soi.tfs,near.tfs] ## Can be a vector in some cases
  collection <- character()
  ## Unfortunately, matrices collapsing to vectors (ns=1 and/or nt=1)
  ## have posed the usual R headaches with losing names
  ## and nrow/ncol evaluating to NaN and so on
  ## ( using data.frame didn't solve this completely )
  ## Hence all the case statements and cumbersome code
  
  if ( ns==1 && nt==1 ){
    partner1 <- soi.tfs
    partner2 <- near.tfs
    if ( distMat==1 ){
      collection <- rbind(character(),c(partner1,partner2))
    }
  }## end ns==1, nt==1 case

  if ( ns>1 && nt==1){
    partner1 <- near.tfs
    one.indices <- which(distMat==1)
    keepers <- one.indices
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(near.tfs,length(keepers)),soi.tfs[keepers]) 
      rownames(collection) <- NULL
    }
  }## end ns>1, nt==1 case

  if ( nt>1 && ns==1){
    partner1 <- soi.tfs
    one.indices <- which(distMat==1)
    keepers <- one.indices
    if ( length(keepers) > 0 ){
      collection <- cbind(rep(soi.tfs,length(keepers)),near.tfs[keepers]) 
      rownames(collection) <- NULL
    }
  }## end nt>1, ns==1 case

  if ( nt>1 &&  ns> 1 ){
    one.indices <- which(distMat==1,arr.ind=TRUE)

    if ( length(one.indices) != 0 ){
      for ( one.index.row  in 1:nrow(one.indices) ){
        partner1 <- soi.tfs[one.indices[one.index.row,1]]
        partner2 <- near.tfs[one.indices[one.index.row,2]]
        collection <- rbind(collection,c(partner1,partner2))
      }
    }
  
  }## end nt>1, ns>1 case

  colnames(collection) <- NULL
  return (collection) 
}



##
## Get nearest neighbors in a distance matrix
## Set either N: number of neighbors OR radius: maximum radius
##

"getExpressionNeighbors" <- function(probeset,distmat,N=NULL,cutoff=NULL){ 
  if ( is.null(N) & is.null(cutoff) ){ stop("No constraint on neighbors")}
  nomos <- rownames(distmat)
  others <- nomos[-which(nomos==probeset)]
  distvec <- distmat[probeset,others]
  rank.indices <- sort(distvec,index.return=TRUE)$ix
  if ( is.null(cutoff) & !is.null(N) ){
    winners <- others[rank.indices[1:N]]
  }
  if ( is.null(N) & !is.null(cutoff) ){
    winners <- (others[rank.indices])[which(distvec[others[rank.indices]] <= cutoff )]
  }
  ## if both are specified, return below cutoff, but not exceeding N
  if ( !is.null(N) & !is.null(cutoff) ){ 
    winners <- (others[rank.indices])[which(distvec[others[rank.indices]] <= cutoff )]
    if ( length(winners) >= N ){
      winners <- others[rank.indices[1:N]]
    }
  }
  return (winners)
}


split1 <- function(splitme,splitchar='_'){ #parse binding site names by dropping "_" and everything behind it
 strsplit(splitme,split=splitchar,fixed=TRUE,extended=FALSE)[[1]][1];
}



## Miniminum distance between elements of a vector
##
## For ordered mode, think of
## "the distance to the closest a upstream of b"
##
## distance returned as a positive number
## fla, flb, Feature length associated with a and b, respectively
minPairDistance <- function(a,b,fla,flb,keep.order=TRUE){
  na <- length(a)
  nb <- length(b)
  min.dist <- NULL
  if ( na==1 | nb== 1){
    mm <- a - b
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(abs(mm))
    } else if ( keep.order & (TRUE %in% (mm <0) )) { # need to test that there actually is negative number
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  } else {
    a.repped <- t(matrix(rep(a,nb),nrow=na,ncol=nb))
    b.repped <- matrix(rep(b,na),nrow=nb,ncol=na)
    mm <- a.repped - b.repped
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(min(abs(mm)))
    } else if ( keep.order & (TRUE %in% (mm <0)) ){
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  }
  return(min.dist)
}



## Miniminum distance between elements of a vector
## Same as above, but fla and flb are vectors
##
## For ordered mode, think of
## "the distance to the closest a upstream of b"
##
## distance returned as a positive number
## fla, flb, Feature length associated with a and b, respectively
minPairDistanceV <- function(a,b,fla,flb,keep.order=TRUE){
  na <- length(a)
  nb <- length(b)
  min.dist <- NULL
  if ( na==1 | nb== 1){
    mm <- a - b
    mm[which(mm>=0 & mm <= fla )] <- Inf
    mm[which(mm<=0 & mm >= -flb )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(abs(mm))
    } else if ( keep.order & (TRUE %in% (mm <0) )) { # need to test that there actually is negative number
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  } else {
    a.repped <- t(matrix(rep(a,nb),nrow=na,ncol=nb))
    b.repped <- matrix(rep(b,na),nrow=nb,ncol=na)
    fla.repped <- t(matrix(rep(fla,nb),nrow=na,ncol=nb))
    flb.repped <- matrix(rep(flb,na),nrow=nb,ncol=na)

    mm <- a.repped - b.repped
    mm[which(mm>=0 & mm <= fla.repped )] <- Inf
    mm[which(mm<=0 & mm >= -flb.repped )] <- -Inf
    if ( !keep.order ){
      min.dist <- min(min(abs(mm)))
    } else if ( keep.order & (TRUE %in% (mm <0)) ){
      min.dist <- abs(max(mm[which(mm<0)]))
    }
  }
  return(min.dist)
}


estimatePrefixMatrixLengths <- function(matrixLength){
  pfxs <- as.character(sapply(names(matrixLength),split1))
  ml <- matrixLength
  names(ml) <- pfxs
  el <- round(tapply(ml,as.factor(names(ml)),mean))
  return(el)
}




countInterLap <- function ( alist ){
  nmems <- length(alist)
  membs <- names(alist)
  olap <- matrix(nrow=nmems,ncol=nmems)
  rownames(olap) <- membs
  colnames(olap) <- membs
  for ( mem1 in membs ){
    for ( mem2 in membs ){
      olap[mem1,mem2] <- length(intersect(alist[[mem1]],alist[[mem2]]))
    }
  }
  return(olap)
}
    
