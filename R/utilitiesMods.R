
## filter by a correlation threshold
## unrolled=TRUE means we apply to full unrolled values
## unrolled=FALSE means we apply to the derivative

filterModsByCorrelation <- function( mods, cor.threshold, unrolled=TRUE, n.cands=1){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    if ( !unrolled ){
      keep.indices <- which(mod$cors >= cor.threshold )
    } else {
      keep.indices <- which(mod$cors.unrolled >= cor.threshold )
    }
    if (length(keep.indices)>0){
      returnList[[psoi]] <- NULL
      if( n.cands==1 ) {
        returnList[[psoi]]$tfs <- mod$tfs[keep.indices]
        returnList[[psoi]]$tfs.cname <- mod$tfs.cname[keep.indices]
        returnList[[psoi]]$betas <- mod$betas[keep.indices]      
      } else if ( n.cands==2 ) {
        ## matrix construct has no effect on matrices, but ensures we get one-row matrices
        returnList[[psoi]]$tfs <- matrix(mod$tfs[keep.indices,],ncol=2)
        returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[keep.indices,],ncol=2)
        returnList[[psoi]]$betas <- matrix(mod$betas[keep.indices,],ncol=2)      
      }
      returnList[[psoi]]$offsets <- mod$offsets[keep.indices]
      returnList[[psoi]]$rmsds <- mod$rmsds[keep.indices]
      returnList[[psoi]]$cors <- mod$cors[keep.indices]
      returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled[keep.indices]
      returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled[keep.indices]
      returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter[keep.indices]
      returnList[[psoi]]$gcv <- mod$gcv[keep.indices]    
      

    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}

  
## filter by an RMSD threshold
## unrolled=TRUE means we apply to full unrolled values
## unrolled=FALSE means we apply to the derivative
## scale=FALSE means we use RMSD 'as is'
## scale=TRUE means we scale RMSD by magnitude of expression. Currently works only for unrolled=TRUE

## new :march 2007 
## rmsd.type should be one of "orig", "unrolled", "fullode"

## unrolled flag to be phased away

filterModsByRMSD <- function( mods, rmsd.threshold, rmsd.type="orig",unrolled=TRUE, scale=FALSE,n.cands=1){

  ## we'll allow a threshold per gene
  ## simple (simplistic) test for how/whether to use that 
  use.gene.specific <- FALSE
  if ( length(rmsd.threshold) > 1 ){
    use.gene.specific <- TRUE
  }
  
  returnList <- list()
  psois <- names(mods)

  if ( use.gene.specific ){
    if ( FALSE %in% ( psois %in% names(rmsd.threshold) ) ){
      stop("You are missing a threshold for a gene\n")
    }
  }
  
  for ( psoi in psois ){

    if ( use.gene.specific ){
      rmsd.thr <- rmsd.threshold[psoi]
    } else {
      rmsd.thr <- rmsd.threshold
    }    
    
    mod <- mods[[psoi]]
    if ( !scale ){
      scalefac <- 1.
    } else {
      scalefac <- max.intensity[psoi]
    }

    if ( rmsd.type == "orig" ){
      keep.indices <- which(mod$rmsds/scalefac < rmsd.thr )
    } else if ( rmsd.type == "unrolled" ){
      keep.indices <- which(mod$rmsds.unrolled/scalefac < rmsd.thr )
    } else if ( rmsd.type == "fullode" ){
      keep.indices <- which(mod$rmsd.ode/scalefac < rmsd.thr )
    }
    
    ##if ( !unrolled ){
    ##  keep.indices <- which(mod$rmsds/scalefac < rmsd.thr )
    ##} else {
    ##  keep.indices <- which(mod$rmsds.unrolled/scalefac < rmsd.thr )
    ##}

    if (length(keep.indices)>0){
      returnList[[psoi]] <- NULL
      if( n.cands==1 ) {
        returnList[[psoi]]$tfs <- mod$tfs[keep.indices]
        returnList[[psoi]]$tfs.cname <- mod$tfs.cname[keep.indices]
        returnList[[psoi]]$betas <- mod$betas[keep.indices]      
      } else if ( n.cands==2 ) {
        returnList[[psoi]]$tfs <- matrix(mod$tfs[keep.indices,],ncol=2)
        returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[keep.indices,],ncol=2)
        returnList[[psoi]]$betas <- matrix(mod$betas[keep.indices,],ncol=2)      
      }
      returnList[[psoi]]$offsets <- mod$offsets[keep.indices]
      returnList[[psoi]]$rmsds <- mod$rmsds[keep.indices]
      returnList[[psoi]]$cors <- mod$cors[keep.indices]    
      returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled[keep.indices]
      returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled[keep.indices]
      returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter[keep.indices]
      returnList[[psoi]]$gcv <- mod$gcv[keep.indices]
      if ( !is.null(mod$binding.partner) ) {
        returnList[[psoi]]$binding.partner <- mod$binding.partner[keep.indices]
      }
      if( n.cands==1 ) {
        returnList[[psoi]]$betas.ode <- mod$betas.ode[keep.indices]      
      } else if ( n.cands==2 ) {
        returnList[[psoi]]$betas.ode <- matrix(mod$betas.ode[keep.indices,],ncol=2)      
      }
      returnList[[psoi]]$rmsd.ode <- mod$rmsd.ode[keep.indices]
      returnList[[psoi]]$tau.ode <- mod$tau.ode[keep.indices]
      

    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}

## filter by an Offset threshold
filterModsByOffset <- function( mods, offset.threshold, n.cands=1){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    keep.indices <- which(abs(mod$offsets) < offset.threshold )
    if (length(keep.indices)>0){
      returnList[[psoi]] <- NULL
      if( n.cands==1 ) {
        returnList[[psoi]]$tfs <- mod$tfs[keep.indices]
        returnList[[psoi]]$tfs.cname <- mod$tfs.cname[keep.indices]
        returnList[[psoi]]$betas <- mod$betas[keep.indices]      
      } else if ( n.cands==2 ) {
        returnList[[psoi]]$tfs <- matrix(mod$tfs[keep.indices,],ncol=2)
        returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[keep.indices,],ncol=2)
        returnList[[psoi]]$betas <- matrix(mod$betas[keep.indices,],ncol=2)      
      }
      returnList[[psoi]]$offsets <- mod$offsets[keep.indices]
      returnList[[psoi]]$rmsds <- mod$rmsds[keep.indices]
      returnList[[psoi]]$cors <- mod$cors[keep.indices]    
    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}
  

## filter to get highest correlation for each target
getHighestCorrelation <- function( mods, n.cands=1){
  returnList <- list()
  psois <- names(mods)
  
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    keep.index <- which.max(mod$cors)
    
    returnList[[psoi]] <- NULL
    if( n.cands==1 ) {
      returnList[[psoi]]$tfs <- mod$tfs[keep.index]
      returnList[[psoi]]$tfs.cname <- mod$tfs.cname[keep.index]
      returnList[[psoi]]$betas <- mod$betas[keep.index]      
    } else if ( n.cands==2 ) {
      returnList[[psoi]]$tfs <- matrix(mod$tfs[keep.index,],ncol=2)
      returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[keep.index,],ncol=2)
      returnList[[psoi]]$betas <- matrix(mod$betas[keep.index,],ncol=2)      
    }
    returnList[[psoi]]$offsets <- mod$offsets[keep.index]
    returnList[[psoi]]$rmsds <- mod$rmsds[keep.index]
    returnList[[psoi]]$cors <- mod$cors[keep.index]    
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}

## For a certain pairs create mods with single inputs
## Input: mods: a pairwise mod structure
## Input: pairs: a matrix of pairs.

pairToSingleMods <- function( mods, pairs){ 

  returnListP <- list()
  returnListS <- list()

  pairs.sorted <- unique(t(apply(pairs, 1, sort))) ## sort first within a row, then finds unique rows
  strpairs <- paste(pairs.sorted[,1],pairs.sorted[,2]) ## each pair (now sorted), represented as single string "tf1 tf2"
  
  psois <- names(mods)
  
  for ( targ in psois ){
    mod <- mods[[targ]]

    n.regpairs <- nrow(mod$tfs)

    for ( i in 1:n.regpairs ){

      tfs <- mod$tfs[i,]
      tfs.as.string <- paste(sort(tfs),collapse=" ")
      
      if ( tfs.as.string %in% strpairs ){

        ## Note!! Keep in mind that the structure *looks* somewhat like a pair structure,
        ## but the interpretation is different
        returnListS[[targ]]$tfs <- c(returnListS[[targ]]$tfs, tfs )
        returnListS[[targ]]$tfs.cname <- c(returnListS[[targ]]$tfs.cname, as.character(cname.compare[tfs]) )

        ## fit parameters will need to be recomputed for each of the two tfs
        ##returnListS[[targ]]$betas <- c(returnListS[[targ]]$betas, as.numeric(coefs.ris) )
        ##returnListS[[targ]]$offsets <- c(returnListS[[targ]]$offsets, as.numeric(ris[1]))
        ##returnListS[[targ]]$rmsds <- c(returnListS[[targ]]$rmsds, rmsd.val)
        ##returnListS[[targ]]$cors <- c(returnListS[[targ]]$cors, cor.val)

        ## Consider whether to retain knowledge of partner
        ## returnListS[[targ]]$binding.partner <- c(returnListS[[targ]]$binding.partner,setdiff(ivars,ivars.keepers))        
      } else { ## if not in set to be treated, just copy results 
        returnListP[[targ]]$tfs <- rbind(returnListP[[targ]]$tfs,matrix(mod$tfs[i,],ncol=2))
        returnListP[[targ]]$tfs.cname <- rbind(returnListP[[targ]]$tfs.cname,matrix(mod$tfs.cname[i,],ncol=2))
        returnListP[[targ]]$betas <- rbind(returnListP[[targ]]$betas,matrix(mod$betas[i,],ncol=2))
        returnListP[[targ]]$offsets <- c(returnListP[[targ]]$offsets,mod$offsets[i])
        returnListP[[targ]]$rmsds <- c(returnListP[[targ]]$rmsds,mod$rmsds[i])
        returnListP[[targ]]$cors <- c(returnListP[[targ]]$cors,mod$cors[i])
        
        if ( !is.null(mod$cors.unrolled[i]) ) {
          returnListP[[targ]]$cors.unrolled <- c(returnListP[[targ]]$cors.unrolled,mod$cors.unrolled[i])
        }
            
        if ( !is.null(mod$rmsds.unrolled[i]) ) {
          returnListP[[targ]]$rmsds.unrolled <- c(returnListP[[targ]]$rmsds.unrolled,mod$rmsds.unrolled[i])
        }
            
        if ( !is.null(mod$taus[i]) )  {
          returnListP[[targ]]$taus <- c(returnListP[[targ]]$taus,mod$taus[i])
        }
        
      }
    } ## end loop over regpairs

    ## Looks like we need to handle multiple occurences in returnListS

    if ( length(returnListS[[targ]]$tfs) > 0 ){
      returnListS[[targ]]$tfs <- unique(sort(returnListS[[targ]]$tfs))
      returnListS[[targ]]$tfs.cname <- as.character(cname.compare[returnListS[[targ]]$tfs])
    }
    ## Cludgy mod
    ## For singles, we really want the correct format for subsequent model fitting. Thus:
    if ( length(returnListS[[targ]]$tfs) > 0 ){
      tempvec <- returnListS[[targ]]$tfs.cname
      returnListS[[targ]]$tfs <- NULL
      returnListS[[targ]]$tfs.cname <- NULL
      returnListS[[targ]] <- NULL
      returnListS[[targ]] <- tempvec
      rm(tempvec)
    }
    
    
    
  } ## end loop over targets


  ## should be checking each for non-zero length, really
  return(list(mods.pairs=returnListP,targsAndCands.singles=returnListS))

}



## prefix: file prefix
## mode:
## np - creates NP file
## cname - creates cname file

## ppints - specify whether to show protein-protein interactions amongst TFs


##For node attribute files
## LOCAL
##tc.dir.linux <- "/users/thorsson/macrophage/TFinfluence/R/REGPREDS/"
##tc.dir.windows <- "\\\\Isb-1\\thorsson\\macrophage\\TFinfluence\\R\\REGPREDS\\"

##tc.dir.linux <- paste("/users/thorsson/macrophage/TFinfluence/R/",outputdir,"/",sep="")

## NET/PUBLIC
##tc.dir.linux <- "/net/public/Vesteinn/RegulationHypotheses/nw1/"
##tc.dir.windows <- "\\\\Isb-2\\PUBLIC\\Vesteinn\\RegulationHypotheses\\nw1\\"


##For node attribute files
## PRION
##tc.dir.linux <- "/users/thorsson/Prion/REGPREDS/"
##tc.dir.windows <- "\\\\Isb-1\\thorsson\\Prion\\REGPREDS\\"


##  outFile1 <- paste(prefix,".sif",sep="")
##  outFile2 <- paste(prefix,".eda",sep="")

##  outFile1 <- paste(tc.dir.linux,prefix,".sif",sep="")
##  outFile2 <- paste(tc.dir.linux,prefix,".eda",sep="")
  
##  outFile3.linux <- paste(tc.dir.linux,"timecourseLinuxPath.noa",sep="") ## repair this fixed reference
##  outFile3.windows <- "Networks/timecourseWindowsPath.noa" ## repair this fixed reference

##  outFile3.linux <- "Networks/timecourseLinuxPath.noa" ## repair this fixed reference
##  outFile3.windows <- "Networks/timecourseWindowsPath.noa" ## repair this fixed reference


##
## prefix: prefix for the .noa, .eda files
## linuxdir: linux representation of directory, full path where model figures will be found
## windir: windows representation of directory, full path where model figures will be found
## both are the same physical space
## full path must include slashes 

## if multimod is TRUE then we plot the multiple possible models for each target
## In this selection, we will NOT have access to the model plots at the nodes

createCytoscapeFiles <- function(mod,prefix,linuxdir=NULL,windir=NULL,ppints=TRUE,multimod=FALSE){

  dir.create(linuxdir) ## Assume that we are running linux, for now
  
  sifFile <- paste(linuxdir,prefix,".sif",sep="") ## Assume that we are running linux, for now
  weightsFile <- paste(linuxdir,prefix,".eda",sep="") ## Assume that we are running linux, for now
  ## remove any pre-existing files with these names
  file.remove(sifFile)
  file.remove(weightsFile)
  if ( !multimod ){
    modelFigFileLinux <- paste(linuxdir,prefix,"ModelFigsLinux.noa",sep="") ## model figure locations
    modelFigFileWindows <- paste(linuxdir,prefix,"ModelFigsWindows.noa",sep="") ## model figure locations
    ## remove any pre-existing files with these names
    file.remove(modelFigFileLinux)
    file.remove(modelFigFileWindows)
  }
  



  ##
  ## protein-protein interactions
  ##
  if ( ppints ){
    tfs.ps <- unique(sort(as.character(unlist(lapply(mod,"[[","tfs"))))) ## tfs, as psois  
    index.matrix <- which(tf.dist.cn[cname.compare[tfs.ps],cname.compare[tfs.ps]]==1,arr.ind=TRUE) ## which of these interact directly
    if ( (length(index.matrix)!=0) & (length(index.matrix)!= 2) ){
      index.matrix <- index.matrix[(index.matrix[,1] > index.matrix[,2]),] ## eliminate duplicated interactions ( tf.dist.cn is symmetric )
    }
    if ( length(index.matrix)!=0 ){
      for ( i in 1:nrow(index.matrix) ) {
        tf1 <- tfs.ps[index.matrix[i,1]]
        tf2 <- tfs.ps[index.matrix[i,2]]
        outstring1 <- paste(np[tf1],"pp",np[tf2],sep=" ")
        cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
      }
    }
  }

  ##
  ## protein-DNA interactions
  ## 

  
  cat("weight\n",file=weightsFile)

  if ( !(multimod) ) { ## These commands are for the one model per target case

    cat("TimeCourseLinuxPath (class=java.net.URL)\n",file=modelFigFileLinux)
    cat("TimeCourseWindowsPath (class=java.net.URL)\n",file=modelFigFileWindows)
    
    for ( targ in names(mod) ) {
      regs <- drop( mod[[targ]]$tfs)
      targ.name <- np[targ]
      ## Write location of model figure
      outstringTC.linux <- paste(targ.name," = file://",linuxdir,targ,"_profile.png",sep="")
      outstringTC.windows <- paste(targ.name," = file:\\\\",windir,targ,"_profile.png",sep="")
      cat(outstringTC.linux,"\n",sep="",file=modelFigFileLinux,append=TRUE)
      cat(outstringTC.windows,"\n",sep="",file=modelFigFileWindows,append=TRUE)
      ## write interactions and weights
      if ( length(regs) > 1 ){
        for ( i in 1:length(regs) ){
          outstring1 <- paste(np[regs[i]],"Protein-DNA",targ.name,sep=" ")
          cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
          outstring1.1 <- paste(np[regs[i]],"(Protein-DNA)",targ.name,sep=" ")
          outstring2 <- paste(outstring1.1,"=",mod[[targ]]$betas[i],sep="")
          cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
        }
      } else if ( length(regs)==1 ) {
        outstring1 <- paste(np[regs],"Protein-DNA",targ.name,sep=" ")
        cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
        outstring1.1 <- paste(np[regs],"(Protein-DNA)",targ.name,sep=" ")
        outstring2 <- paste(outstring1.1,"=",mod[[targ]]$betas,sep="")
        cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
      }
    } ## end loop over targets
  } else { ## These commands are for the one model per target case 
    
    for ( targ in names(mod) ) {
      
      targ.name <- np[targ]
      
      tfs <- mod[[targ]]$tfs
      n.tfs <- nrow(tfs) ############### WON'T WORK FOR SINGLE REGULATORS ########################
      
      for ( i in 1:n.tfs ){
        
        regs <- tfs[i,]
 
        ## write interactions and weights
        if ( length(regs) > 1 ){
          for ( i in 1:length(regs) ){
            outstring1 <- paste(np[regs[i]],"Protein-DNA",targ.name,sep=" ")
            cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
            outstring1.1 <- paste(np[regs[i]],"(Protein-DNA)",targ.name,sep=" ")
            outstring2 <- paste(outstring1.1,"=",mod[[targ]]$betas[i],sep="")
            cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
          }
        } else if ( length(regs)==1 ) {
          outstring1 <- paste(np[regs],"Protein-DNA",targ.name,sep=" ")
          cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
          outstring1.1 <- paste(np[regs],"(Protein-DNA)",targ.name,sep=" ")
          outstring2 <- paste(outstring1.1,"=",mod[[targ]]$betas,sep="")
          cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
        }
      } ## end loop over tfs for this target
        
    } ## end loop over targets

  } ## end if over multimod case

}




## psois are assumed to be known
## as is boost.vec,tau, s, t.index.max, t.index.min


## savedir should be full path

predPlotWrapper <- function( psois, mods,savedir="REGPREDS" ) {

  ### KEEP AN EYE ON THIS VARIABLE ##########
  zero.offset <- TRUE
  
  psois <- psois[sort(as.character(cname.compare[psois]),index.return=TRUE)$ix]
  ##savedir <- "REGPREDS"
  dir.create(savedir)
  ##outMat <- cbind(psois,cname.compare[psois],ncbiID[psois])
  ##colnames(outMat) <- c("ProbeSetID","MGNC","GeneID")
  ##outMatFileName <-  paste("REGPREDS/outMat.tsv",sep="")
  ##write.table(outMat,file=outMatFileName,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

  for ( psoi in psois ) {
    ##cat("psoi:",psoi,"\n")
    ivars <- drop(mods[[psoi]]$tfs)
    betas <- drop(mods[[psoi]]$betas)
    names(betas) <- ivars
    wt.predicted.fullscale <- calculateTimeEvolPrediction.3(psoi,ivars,tau,s=1,boost.vec,zero.offset)
    imageFileName <- paste(savedir,"/",psoi,"_profile.png",sep="")
    png(imageFileName,height=800,width=1100)
    plotPrediction.3(psoi,betas,wt.predicted.fullscale)
    dev.off()
  }

}

## A version of plotPrediction.3 with some minor adjustments to
## plot parameters needed to do multiple plots on a single page
plotPrediction.hack <- function(target.pid,pred.coefs,pred.out,main=""){

  pred.pids <- names(pred.coefs)
  coefs <- as.numeric(pred.coefs)
  
  array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)
  
  profile.true <- lps.mus[target.pid,1:t.index.max] ## we'll plot the initial values, even if they were not use for the fit

  ## rescale profile of predicted TFs (scale determined by target intensity ) for plotting purposes
  profile.pred <- 0.5 * max.intensity[target.pid] * lps.mat.max1[pred.pids,1:t.index.max]

  ymax <- max(c(profile.true,pred.out))
  ##xvals <- array.times[1:t.index.max]/60.
  xvals <- 1:t.index.max
  
  ##op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='',ylab='',main=main,axes=FALSE)
  ##plot(xvals,profile.true[1:t.index.max], ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',main=main,lwd=3,axes=FALSE)
  
  axislabels <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
  axis(1,1:t.index.max,labels=axislabels[1:t.index.max])
  axis(2)

  lines(xvals[(t.index.min+1):t.index.max],pred.out, col='blue',lty=2)
  ##lines(xvals[(t.index.min+1):t.index.max],pred.out, col='blue',lty=2,lwd=3)

  colorvec = c('red','magenta','pink','cyan','green','black')
  col <- 1:10
  lty <- 1:6

  col <- rep(col, length.out=length(pred.pids))
  lty <- rep(lty, length.out=length(pred.pids))

  if ( length(pred.pids) > 1 ){
    for ( i in 1:length(pred.pids) ){
      lines( xvals, profile.pred[i,1:t.index.max], col=colorvec[i],lty=1)
##      lines( xvals, profile.pred[i,1:t.index.max], col=colorvec[i],lty=1,lwd=3)
    }
  }
  if ( length(pred.pids) == 1 ){
    lines( xvals, profile.pred[1:t.index.max], col=colorvec[1],lty=1)
    ##lines( xvals, profile.pred[1:t.index.max], col=colorvec[1],lty=1,lwd=3)
  }

  s1 <- gid2gname(pred.pids)
  s2 <- as.character(round(coefs,digits=1))
  s2 <- paste("(",s2,")",sep="")
  legend <- paste(s1,s2)
  
  lg <- c(1.5,0.9*max.intensity[target.pid] )
  legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1 )
  ##legend(lg[1],lg[2], legend=legend, col=colorvec, lty=1,lwd=3)
  

  ##par(op)

}


## Plot 'contact sheets' with predictions
## filename must include full path
gridPredPlotWrapper <- function(mods,file="grid.ps",nx=1,ny=1,tindex.min=1,t.index.max=11,tau=600./log(2.) ) {

  ### KEEP AN EYE ON THIS VARIABLE ##########
  zero.offset <- TRUE
  
  psois <- names(mods)

  ## could do an input parameter check
  
  ## psois follow order presented to this function 
  postscript(file,paper='letter')
  ##png(file,paper='letter',height=800,width=1100)
  
  ##def.par <- par(no.readonly = TRUE)
  ##nf <- layout(matrix(1:(nx*ny),nx,ny,byrow=TRUE),widths=lcm(rep(0.01,nx)),heights=lcm(rep(0.01,ny)))

  par(mfrow=c(nx,ny),mar=(c(3, 2, 2, 1)+0.1),mgp=c(1.75, 0.75, 0))
  par(mfrow=c(nx,ny),mar=(c(2, 2, 1, 1)+0.1),mgp=c(1.75, 0.75, 0))

  counter <- 0
  for ( psoi in psois ) {
##postscript(file,paper='letter')
    ivars <- drop(mods[[psoi]]$tfs)
    betas <- drop(mods[[psoi]]$betas)
    offset <- drop(mods[[psoi]]$offsets)
    
    names(betas) <- ivars

    dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,psoi,ivars,boost.vec,lps.mat.max1)
    outpred <- calculateOutpred(dataToFit$inputs,offset,betas)
    start.val <- lps.mat.max1[psoi,t.index.min]
    scale.fac <- max.intensity[psoi]
    wt.predicted.fullscale <- calculateTimeEvolPrediction.6(t.index.min,t.index.max,tau,outpred,start.val,scale.fac)
    ##cat(wt.predicted.fullscale,"\n")
    
    counter <- counter+1
    main=paste(counter,": ",gid2gname(psoi),":",psoi)
   
    plotPrediction.hack(psoi,betas,wt.predicted.fullscale,main=main)
  }

  dev.off()
 
  
}





## Count number of targets for each pair of regulators 
## ( May alter this so that it is sign-specific )

countTargets <- function( mods ){
  returnMatrix <- matrix(0,nrow=length(tte),ncol=length(tte ) )
  rownames(returnMatrix)=tte
  colnames(returnMatrix)=tte
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]

    n.pairs <- nrow(mod$tfs)
    for ( i in 1:n.pairs ) {

      tf1 <- mod$tfs[i,1]
      tf2 <- mod$tfs[i,2]
      returnMatrix[tf1,tf2] <-  returnMatrix[tf1,tf2]+1
    }
  }
  return(returnMatrix)
}



## Strip empty mods. (These should not really have been generated in the first place!)
removeEmpties <- function ( mods ){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    if ( length(mod) > 0 ) {
      returnList[[psoi]]$tfs <- mod$tfs
      returnList[[psoi]]$tfs.cname <- mod$tfs.cname
      returnList[[psoi]]$betas <- mod$betas
      returnList[[psoi]]$offsets <- mod$offsets
      returnList[[psoi]]$rmsds <- mod$rmsds
      returnList[[psoi]]$cors <- mod$cors
    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }  
}

## Filter by, tf (probesets)
##
## tf.sign is "" if not used
## "pos" for beta >0
## type, is either ode or fd:finite difference result
filterByTF <- function( mods, tf, tf.sign="",type="ode" ){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    n.tfs <- length(mod$tfs)
    tf.logical <- (mod$tfs == tf)
    if ( identical(tf.sign,"") ){
      vec.logical <- tf.logical
    } else {
      ##cat(mod$betas,"\n")
      if ( type == "fd" ){
        sign.string <- c("neg","pos")[sign(mod$betas)/2+3/2]
      } else if ( type == "ode" ){
        sign.string <- c("neg","pos")[sign(mod$betas.ode)/2+3/2]
      }
      sign.logical <-  sign.string == tf.sign
      ##cat(sign.string,"\n")
      ##cat(sign.logical,"\n")
      vec.logical <- tf.logical & sign.logical
    }
    keep.indices <- (1:n.tfs)[vec.logical]
    
    if (length(keep.indices)>0 ){
      
#######################################################################
## THE FOLLOWING IS A HACK UNTIL WE FIGURE OUT HOW TO HANDLE THESE
## RETAIN ONLY FIRST INSTANCE
      keep.indices <- min(keep.indices)
######################################################################

      returnList[[psoi]] <- NULL
      returnList[[psoi]]$tfs <- mod$tfs[keep.indices]
      returnList[[psoi]]$tfs.cname <- mod$tfs.cname[keep.indices]
      returnList[[psoi]]$betas <- mod$betas[keep.indices]      
      returnList[[psoi]]$offsets <- mod$offsets[keep.indices]
      returnList[[psoi]]$rmsds <- mod$rmsds[keep.indices]
      returnList[[psoi]]$cors <- mod$cors[keep.indices]    
      returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled[keep.indices]
      returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled[keep.indices]
      returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter[keep.indices]
      returnList[[psoi]]$gcv <- mod$gcv[keep.indices]
      if ( !is.null(mod$binding.partner) )  {
        returnList[[psoi]]$binding.partner <- mod$binding.partner[keep.indices]
      }
      if ( !is.null (mod$betas.ode) ){
        returnList[[psoi]]$betas.ode <- matrix(mod$betas.ode[keep.indices],ncol=2) 
      }
      if ( !is.null (mod$tau.ode) ){
        returnList[[psoi]]$tau.ode <- mod$tau.ode[keep.indices]
      }
      if ( !is.null (mod$rmsd.ode) ){
        returnList[[psoi]]$rmsd.ode <- mod$rmsd.ode[keep.indices]
      }
      if ( !is.null (mod$tau.ode) ){
        returnList[[psoi]]$tau.ode <- mod$tau.ode[keep.indices]
      }


      
    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}
  
## Filter by pair, tf1, tf2 (probesets)
##
## tf.sign is "" if not used, i.e. we don't care about sign
## c("pos","neg") for beta.tf1>0 tf2.beta<0 etc )
## type, is either ode or fd:finite difference result
filterByTFPair <- function( mods, tf1, tf2, tf.sign="",type="ode" ){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    n.pairs <- nrow(mod$tfs)
    for ( i in 1:n.pairs ) {
      teef1 <- mod$tfs[i,1]
      teef2 <- mod$tfs[i,2]
      if ( (teef1==tf1)&(teef2==tf2) |(teef1==tf2)&(teef2==tf1) ){ ## could use sort here 
        if ( identical(tf.sign,"") ){
          returnList[[psoi]]$tfs <- matrix(mod$tfs[i,],ncol=2)
          returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[i,],ncol=2)
          returnList[[psoi]]$betas <- matrix(mod$betas[i,],ncol=2)      
          returnList[[psoi]]$offsets <- mod$offsets[i]
          returnList[[psoi]]$rmsds <- mod$rmsds[i]
          returnList[[psoi]]$cors <- mod$cors[i]
          returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled[i]
          returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled[i]
          returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter[i]
          returnList[[psoi]]$gcv <- mod$gcv[i]    
        }  else {
          if ( type == "fd" ){
            b1 <- mod$betas[i,1]
            b2 <- mod$betas[i,2]
          } else if ( type == "ode" ) {
            b1 <- mod$betas.ode[i,1]
            b2 <- mod$betas.ode[i,2]
          }
          if ( b1 >= 0 ){ b1.str="pos" }
          if ( b1 < 0  ){ b1.str="neg" }
          if ( b2 >= 0 ){ b2.str="pos" }
          if ( b2 < 0  ){ b2.str="neg" }
          logic1 <- (teef1==tf1)&(teef2==tf2)&(b1.str==tf.sign[1])&(b2.str==tf.sign[2])
          logic2 <- (teef1==tf2)&(teef2==tf1)&(b1.str==tf.sign[2])&(b2.str==tf.sign[1])
          if ( logic1 | logic2 ){
            returnList[[psoi]]$tfs <- matrix(mod$tfs[i,],ncol=2)
            returnList[[psoi]]$tfs.cname <- matrix(mod$tfs.cname[i,],ncol=2)
            returnList[[psoi]]$betas <- matrix(mod$betas[i,],ncol=2)      
            returnList[[psoi]]$offsets <- mod$offsets[i]
            returnList[[psoi]]$rmsds <- mod$rmsds[i]
            returnList[[psoi]]$cors <- mod$cors[i]
            returnList[[psoi]]$rmsds.unrolled <- mod$rmsds.unrolled[i]
            returnList[[psoi]]$cors.unrolled <- mod$cors.unrolled[i]
            returnList[[psoi]]$shrinkage.parameter <- mod$shrinkage.parameter[i]
            returnList[[psoi]]$gcv <- mod$gcv[i]    

            if ( !is.null (mod$betas.ode) ){
              returnList[[psoi]]$betas.ode <- matrix(mod$betas.ode[i,],ncol=2) 
            }
            if ( !is.null (mod$tau.ode) ){
              returnList[[psoi]]$tau.ode <- mod$tau.ode[i]
            }
            if ( !is.null (mod$rmsd.ode) ){
              returnList[[psoi]]$rmsd.ode <- mod$rmsd.ode[i]
            }
            if ( !is.null (mod$tau.ode) ){
              returnList[[psoi]]$tau.ode <- mod$tau.ode[i]
            }

          }

        }
          
      }
    }
  }

  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}




## Weed out entries with TFpair, tf1, tf2 (probesets)
removeTFPair <- function( mods, tf1, tf2 ){
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){
    mod <- mods[[psoi]]
    n.pairs <- nrow(mod$tfs)
    for ( i in 1:n.pairs ) {
      teef1 <- mod$tfs[i,1]
      teef2 <- mod$tfs[i,2]
      if ( !((teef1==tf1)&(teef2==tf2) |(teef1==tf2)&(teef2==tf1)) ){ ## the not of the condition for the filterByTFPair
        returnList[[psoi]]$tfs <- rbind(returnList[[psoi]]$tfs,matrix(mod$tfs[i,],ncol=2))

        returnList[[psoi]]$tfs.cname <- rbind(returnList[[psoi]]$tfs.cname,matrix(mod$tfs.cname[i,],ncol=2))

        returnList[[psoi]]$betas <- rbind(returnList[[psoi]]$betas,matrix(mod$betas[i,],ncol=2))
        returnList[[psoi]]$offsets <- c(returnList[[psoi]]$offsets,mod$offsets[i])
        returnList[[psoi]]$rmsds <- c(returnList[[psoi]]$rmsds,mod$rmsds[i])
        returnList[[psoi]]$cors <- c(returnList[[psoi]]$cors,mod$cors[i])
      }
    }
  }
  
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}


##
## Level predictions for different ligands
##
## Currently works only for single model per psoi

modelPredictStim <- function( mods, t.index.min=4, t.index.max=6, ligand=NULL ){
  
  expmat.string <- paste(ligand,".mat.max1",sep="")
  if ( !(expmat.string %in% ls(.GlobalEnv) )) {
    stop("Illegal condition name")
  }
  expmat <- eval(parse(text=expmat.string))

  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois[1:3] ){
    mod <- mods[[psoi]]
    n.pairs <- nrow(mod$tfs)
    for ( i in 1:n.pairs ) {
      ivars <- mod$tfs[i,]
      coefs <- mod$betas[i,]
      offset <- mod$offsets[i]
      result <- calculateTimeEvolPrediction.5(t.index.min,t.index.max,coefs,offset,targ,ivars,tau,s=1.,boost.vec, expmat )
      returnList[[psoi]]$pred <- rbind(returnList[[psoi]]$pred,result)
    }
  }
  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}

## Currently works only for single model per psoi
modelPredictMultiStim <- function( mods, t.index.min=4, t.index.max=6, ligands=NULL ){

  if ( is.null(ligands) ){
    stop("You forgot to specify ligands")
  }


  for ( ligand in ligands ){
    expmat.string <- paste(ligand,".mat.max1",sep="")
    if ( !(expmat.string %in% ls(.GlobalEnv) )) {
      stop("Illegal condition name:",ligand)
    }
  }
  
  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){

    targ <- psoi
    mod <- mods[[psoi]]
    ivars <- mod$tfs[1,]
    coefs <- mod$betas[1,]
    offset <- mod$offsets

    returnList[[psoi]]$pred <- NULL
    for ( ligand in ligands ){
      expmat.string <- paste(ligand,".mat.max1",sep="")
      expmat <- eval(parse(text=expmat.string))      
      result <- calculateTimeEvolPrediction.5(t.index.min,t.index.max,coefs,offset,targ,ivars,tau,s=1.,boost.vec, expmat )
      returnList[[psoi]]$pred <- rbind(returnList[[psoi]]$pred, result  )
    }
    rownames(returnList[[psoi]]$pred) <- ligands

  }

  if ( length(returnList) > 0 ) {
    return(returnList)
  } else {
    return(NULL)
  }
}

## Jan 2007
## Perform shrinkage model selection on models
##
##  returns two lists, one for double inputs, one for single inputs, after selection

## returnListP - pairs
## returnListS - singles

## tau <- 0 as input is reserved, and taken to mean that taus are in the model

shrinkageSelectModels <- function(inmods, t.index.min, t.index.max, tau, boost.vec, mus.mat.maxnormed, n.cands=2, zero.offset=FALSE ){

  use.tau.per.model <- FALSE
  if ( tau == 0 ){
    use.tau.per.model <- TRUE
  }
  
  returnListP <- list()
  returnListS <- list()

  psois <- names( inmods )
  
  counter <- 0
  for ( targ in psois ){

    counter <- counter + 1
    write(c(targ,counter,counter/length(psois)),file='status.txt'); #report status

    ##cat("You are here one\n")
    
    ## One candidate
    if ( n.cands==1 ){
      
      ##cat("You are here two\n")
      
      tfs <- inmods[[targ]]$tfs
      n.tfs <- length(tfs) 
      for ( i in 1:n.tfs ) {
        ivars <- repProbes.cname[tfs[i]] ## ivars has length one

        if ( !(targ %in% ivars ) ){

          dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,mus.mat.maxnormed)
          result <- getBestS(dataToFit,zero.offset)
          rel.bound.min.gcv <- result$s
          gcv.value <- result$s
          
          ## Refit model to OLS. s specificies only the complexity
          ##mrf <- modReFit(dataToFit,s=rel.bound.min.gcv,zero.oeffet=zero.offset)
          mrf <- modReFit(dataToFit,s=rel.bound.min.gcv,zero.offset,targ,t.index.min,t.index.max,tau)
          ivars.keepers <- mrf$tfs

          ##cat(ivars.keepers,"2\n")
          
          if ( length(ivars.keepers) == 1 ){
            returnListS[[targ]]$tfs <- c(returnListS[[targ]]$tfs, as.character(ivars.keepers) )
            returnListS[[targ]]$tfs.cname <- c(returnListS[[targ]]$tfs.cname, as.character(cname.compare[ivars.keepers]) )
            returnListS[[targ]]$betas <- c(returnListS[[targ]]$betas, mrf$betas ) 
            returnListS[[targ]]$offsets <- c(returnListS[[targ]]$offsets, mrf$offset )
            returnListS[[targ]]$rmsds <- c(returnListS[[targ]]$rmsds, mrf$rmsd )
            returnListS[[targ]]$cors <- c(returnListS[[targ]]$cors, mrf$cor )
            returnListS[[targ]]$rmsds.unrolled <- c(returnListS[[targ]]$rmsds.unrolled, mrf$rmsd.unrolled )
            returnListS[[targ]]$cors.unrolled <- c(returnListS[[targ]]$cors.unrolled, mrf$cor.unrolled )
            returnListS[[targ]]$shrinkage.parameter <- c(returnListS[[targ]]$shrinkage.parameter,rel.bound.min.gcv)
            returnListS[[targ]]$gcv <- c(returnListS[[targ]]$gcv,gcv.value)
          } 
          
        }## complete if ( targ not in ivars ) 
      } ## complete loop over tfs
    }## close cands=1 section

    ## Two candidates-
    
    if ( n.cands==2 ){

      ##cat("You are here three\n")
      
      mod <- inmods[[targ]]
      tf.pairs <- mod$tfs
      n.pairs <- nrow(tf.pairs) ## This does seem to work for '1 row' case
      for ( i in 1:n.pairs ) { ## Loop over pairs
        ivars <- tf.pairs[i,]

        ##cat("Looking at:",cc[ivars],"\n")
 
        if ( !(targ %in% ivars ) ){ ## ignore if self-regulating

          if (  use.tau.per.model ){
            tau <- mod$taus[i]
          }

          dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,mus.mat.maxnormed)

          ##cat("You are here 3.1\n")          
          rel.bound.min.gcv <- getBestS(dataToFit,zero.offset)

          result <- getBestS(dataToFit,zero.offset)
          rel.bound.min.gcv <- result$s
          gcv.value <- result$gcv
          
          ##cat("You are here 3.2\n")
          ##cat("rel.bound.min.gcv:",rel.bound.min.gcv, "\n")
          if ( rel.bound.min.gcv > 0 ){
            ##cat("You are here 3.3\n")
            ##cat("Targ:", targ, "\n")
            ##cat(str(dataToFit),"\n")
            ##save(dataToFit,file="dataToFit.RData")
            
            ## Refit model to OLS. s specificies only the complexity
            mrf <- modReFit(dataToFit,s=rel.bound.min.gcv,zero.offset,targ,t.index.min,t.index.max,tau)
            ## cat("You are here four\n")
            ivars.keepers <- mrf$tfs
            ##cat("You are here five\n")
            ##cat(ivars.keepers,"2\n")
          } else {
            ivars.keepers <- NULL
          }
            
          ## If both variables are retained
          if ( length(ivars.keepers) == 2 ) {
            ##cat("Keeping both\n\n")
            returnListP[[targ]]$tfs <- rbind(returnListP[[targ]]$tfs, as.character(ivars.keepers) )
            returnListP[[targ]]$tfs.cname <- rbind(returnListP[[targ]]$tfs.cname, as.character(cname.compare[ivars.keepers]) )
            returnListP[[targ]]$betas <- rbind(returnListP[[targ]]$betas,mrf$betas)
            returnListP[[targ]]$offsets <- c(returnListP[[targ]]$offsets, mrf$offset)
            if (  use.tau.per.model ){
              returnListP[[targ]]$taus <- c(returnListP[[targ]]$taus, tau)
            }
            returnListP[[targ]]$rmsds <- c(returnListP[[targ]]$rmsds, mrf$rmsd)
            returnListP[[targ]]$cors <- c(returnListP[[targ]]$cors, mrf$cor )
            returnListP[[targ]]$rmsds.unrolled <- c(returnListP[[targ]]$rmsds.unrolled, mrf$rmsd.unrolled )
            returnListP[[targ]]$cors.unrolled <- c(returnListP[[targ]]$cors.unrolled, mrf$cor.unrolled )
            returnListP[[targ]]$shrinkage.parameter <- c(returnListP[[targ]]$shrinkage.parameter,rel.bound.min.gcv)
            returnListP[[targ]]$gcv <- c(returnListP[[targ]]$gcv,gcv.value)
          }

          ## If one variable is retained
          if ( length(ivars.keepers) == 1 ){
            ##cat("=========Keeping:",cc[ivars.keepers],"\n\n")
            returnListS[[targ]]$tfs <- c(returnListS[[targ]]$tfs, as.character(ivars.keepers) )
            returnListS[[targ]]$tfs.cname <- c(returnListS[[targ]]$tfs.cname, as.character(cname.compare[ivars.keepers]) )
            returnListS[[targ]]$betas <- c(returnListS[[targ]]$betas, mrf$betas ) 
            returnListS[[targ]]$offsets <- c(returnListS[[targ]]$offsets, mrf$offset )
            if (  use.tau.per.model ){
              returnListS[[targ]]$taus <- c(returnListS[[targ]]$taus, tau)
            }
            returnListS[[targ]]$rmsds <- c(returnListS[[targ]]$rmsds, mrf$rmsd )
            returnListS[[targ]]$cors <- c(returnListS[[targ]]$cors, mrf$cor )
            returnListS[[targ]]$rmsds.unrolled <- c(returnListS[[targ]]$rmsds.unrolled, mrf$rmsd.unrolled )
            returnListS[[targ]]$cors.unrolled <- c(returnListS[[targ]]$cors.unrolled, mrf$cor.unrolled )
            returnListS[[targ]]$shrinkage.parameter <- c(returnListS[[targ]]$shrinkage.parameter,rel.bound.min.gcv)
            ## just in case, keep track of the non-influencing variable
            returnListS[[targ]]$binding.partner <- c(returnListS[[targ]]$binding.partner,setdiff(ivars,ivars.keepers))
            returnListS[[targ]]$gcv <- c(returnListS[[targ]]$gcv,gcv.value)
          } 
        }## complete if ( targ not in ivars ) 
      } ## complete loop over pairs
    }## close cands=2 section
  } ## close loop over targets
  
  return(list(mods.pairs=returnListP,mods.singles=returnListS))
  
}


##
## Utilities for ModelVec structures
##
## List: One entry for each model structure
## Indexed by ordinal

sliceModVecByTargs <- function( psois, modvec ){
  inds <- which(unlist(lapply(modvec,"[[","targ")) %in% psois)
  return( modvec[inds] )
}


## slice mod vec by TFs
sliceModVecByTFs <- function( ivars, modvec){
  intersect.len <- lapply(lapply(lapply(modvec,"[[","tfs"),intersect,ivars),length)
  does.match <- which( unlist(lapply(intersect.len,"==",length(ivars))))
  return( modvec[does.match] )
}

## Similar to slicing by TFs, but need to match sign
## give signs as "pos", "neg"
## operates in two modes
## comparison.mode = "exact", number of regulators must match.
## comparison.mode = "subset", returns target whose regulation contains the specified regulators
## so far, have only implemented "exact"

sliceModVecByRegulation <- function( ivars, ivars.sign, modvec, comparison.mode="exact"){
  
  ## alphabetic sort, by TF, to be used later for comparison
  sort.order <- sort(ivars,index.return=TRUE)$ix
  ivars <- as.character(ivars[sort.order])
  ivars.sign <- as.character(ivars.sign[sort.order])

  returnList <- list()
  for ( mod in modvec ){
    tfs <- as.character(mod$tfs)
    if ( "betas.ode" %in% names(mod) ){
      betas <- as.numeric(mod$betas.ode)
    } else {
      betas <- as.numeric(mod$betas)
    }
    tfs.sign <- c("neg","pos")[( betas >= 0 ) + 1]
    ## sort as above, to facilitate comparison
    if ( length(tfs) > 1 ){
      sort.order <- sort(tfs,index.return=TRUE)$ix
      tfs <- tfs[sort.order]
      tfs.sign <- tfs.sign[sort.order]
    }
    ## Legwork done, let's do the comparison
    if ( comparison.mode == "exact" ){
      cond1 <- identical(ivars,tfs)
      cond2 <- identical(ivars.sign,tfs.sign)
      if ( cond1 & cond2 ){
        returnList[[length(returnList)+1]] <- mod
      }
    } else if ( comparison.mode == "subset" ){
      ##other stuff
    }
  } ## end loop over mods

  return( returnList )
}



##
## Create modvec from target list + fixed regulators (one or two )
##

## regs should be in psoi form. zero.offset,tau.estimate set to false

createModVecFromTargRegLists <- function( evidenceString, targs, ivars, t.index.min,t.index.max,tau,boost.vec){
  regs <- as.character(cc[ivars])
  n.regs <- length(regs)
  if ( n.regs == 1 ){ reg = regs }
  pdna <- list()
  for ( targ in targs ){
    if ( n.regs == 1 ){
      pdna[[targ]] <- reg
    }
  }
  mods <- calcModels(pdna,t.index.min,t.index.max,tau,boost.vec,lps.mat.max1,n.cands=n.regs,zero.offset=FALSE,tau.estimate=FALSE)
  return.modvec <- createModVecFromMultMod(evidenceString,mods, n.cands=n.regs)
  return(return.modvec)
}


## For sites predicted not with kinetics, we create fake models
## so that we can proceed with plotting in BioTapestry
## BioTapestry needs to know about the sign of beta
createModVecNoKinetics <- function(evidenceString,targs,ivars){
  regs <- as.character(cc[ivars])
  n.regs <- length(regs)
  if ( n.regs == 1 ){ reg = regs }
  fmods <- list()
  for ( targ in targs ){
    fmods[[targ]] <- list()
    if ( n.regs == 1 ){
      fmods[[targ]]$tfs <- as.character(rc[reg])
      fmods[[targ]]$tfs.cname <- reg
      fmods[[targ]]$betas <- 1.11111111111 ## meant to be recognizable
    }
  }
  return.modvec <- createModVecFromMultMod(evidenceString,fmods,n.cands=n.regs)
  return(return.modvec)
}
  
                       
##
## Create model vector for smods structure
## The difference lies only in an evidence code
##

createModVecFromSmods <- function (evidencelabel, smodstructure ){
  counter <- 0 
  outList <- list()
  for ( smod in smodstructure ){
    tempstruct <- smod
    tempstruct[["evidence"]] <- evidencelabel
    tempstruct[["tfs"]] <- smod[["selected"]]
    counter <- counter + 1
    outList[[counter]] <- tempstruct
  }
  return( outList )
}

##
## create ModVec from
## structure with multiple models per target
##

createModVecFromMultMod <- function( evidencelabel, multmods, n.cands=1){
  
  ## refers to the structure in the original (?) object
  label.class1 <- c("tau.ode")
  label.class2 <- c("offsets","rmsds","cors","rmsds.unrolled","cors.unrolled","shrinkage.parameter","gcv","rmsd.ode")
  label.class3 <- c("tfs","tfs.cname","betas","betas.ode")
  
  ntargs <- length(multmods)
  targs <- names(multmods)
  
  outList <- list()
  counter <- 0 
  
  for ( i  in 1:ntargs ){
    
    targ <- targs[i]
    mod <- multmods[[i]]

    if (n.cands==1){
      n.regs <- length(mod$tfs)
    } else if (n.cands==2) {
      n.regs <- nrow(mod$tfs)
    }
    labels <- names(mod) ## the various attributes of the models that are present
    
    tempList <- list()
    
    for ( k in 1:n.regs ){
      tempList[[k]] <- list()
      tempList[[k]][["targ"]] <- targ
    }
    
    if (n.cands==1){
      for ( label in labels ){
        if ( label %in% c(label.class2,label.class3) ){
          for ( k in 1:n.regs ){
            element <- mod[[label]][k]
            tempList[[k]][[label]] <- element
          }
        } else if ( label %in% label.class1 ){
          for ( k in 1:n.regs ){
            tempList[[k]][[label]] <- mod[[label]][k]
          }
        } else {
          cat("Warning: Missing label", label, "\n" )
        }
      }
    } else if (n.cands==2) {
      for ( label in labels ){
        if ( label %in% c(label.class1,label.class2) ){
          for ( k in 1:n.regs ){
            tempList[[k]][[label]] <- mod[[label]][k]
          }
        } else if ( label %in% label.class3 ) {
          for ( k in 1:n.regs ){
            tempList[[k]][[label]] <- mod[[label]][k,]
          }
        }
      }
    } ## end if cases for n.cands
    
    ## Evidence label, could perhaps be placed higher up
    for ( k in 1:n.regs ){
      tempList[[k]][["evidence"]] <- evidencelabel
    }
    outList <- c(outList,tempList)
    
  } ## end loop over targets
  
  return( outList)
  
}


## Report statistics on a modvec
modVecStats <- function( modvec ){

  statVec <- numeric() ## A named vector with statstics 
  statList <- list() ## list object for other summaries that are not numeric

  ## Target statistics
  targ.vec <- unlist(lapply(modvec,"[[","targ"))
  statVec["ModCount"] <- length(targ.vec)
  statVec["TargCount"] <- length(unique(targ.vec))

  ## TFs (flattened, no regard for network structure , etc. )
  tf.vec <- unlist(lapply(modvec,"[[","tfs"))
  statVec["EdgeCount"] <- length(tf.vec)
  statVec["RegCount"] <- length(unique(tf.vec))
  
  tf.count.table <- sort(table(tf.vec),decreasing=TRUE)
  ## Extra transformation, as table does not return a proper named numeric vector
  name.save <- cc[names(tf.count.table)] ## easier to look at! 
  tf.count.table <- as.vector(tf.count.table)
  names(tf.count.table) <- name.save
  
  ## report major regulators
  min.targs <- 10 
  tf.count.table <- tf.count.table[which(tf.count.table > min.targs )]
  
  return(list(statVec=statVec,majorTFs=tf.count.table))
  
}
##


##
## Report statistics comparing two modvecs
## 
compareModVecs <- function( modvec.A, modvec.B ){

  statVec <- numeric()

  ## Target statistics
  targ.vec.A <- unlist(lapply(modvec.A,"[[","targ"))
  targ.vec.B <- unlist(lapply(modvec.B,"[[","targ"))

  ## We'll use the signs from the ode fit if available
  if ( "betas.ode" %in% names(modvec.A[[1]]) ){
    bstring.A <- "betas.ode"
  } else {
    bstring.A <- "betas"
  }
  if ( "betas.ode" %in% names(modvec.B[[1]]) ){
    bstring.B <- "betas.ode"
  } else {
    bstring.B <- "betas"
  }
  
  mat.A <- flattenModVec(modvec.A,bstring.A)
  mat.B <- flattenModVec(modvec.B,bstring.B)

  n.edges.A <- nrow(mat.A)
  n.edges.B <- nrow(mat.B)

  ## (check order of TF and target)
  mv.strings.A.noSign <- apply(mat.A[,c("Target","TF")],1,paste,collapse="-")
  mv.strings.B.noSign <- apply(mat.B[,c("Target","TF")],1,paste,collapse="-")
  n.edges.in.common.noSign <- length(intersect( mv.strings.A.noSign,  mv.strings.B.noSign ))
  
  mv.strings.A.wSign <- apply(mat.A[,c("Target","TF","Sign")],1,paste,collapse="-")
  mv.strings.B.wSign <- apply(mat.B[,c("Target","TF","Sign")],1,paste,collapse="-")
  n.edges.in.common.wSign <- length(intersect( mv.strings.A.wSign,  mv.strings.B.wSign ))

  statVec[["A Edges"]] <-  n.edges.A
  statVec[["B Edges"]] <-  n.edges.B
  statVec["Targets in both A and B"] <- length(intersect(unique(targ.vec.A),unique(targ.vec.B)))
  statVec[["Shared Edges"]] <- n.edges.in.common.noSign
  statVec[["Shared Edges and Sign"]] <- n.edges.in.common.wSign

  return(statVec)
  
}


### Return simple matrix representation of modvec
## construct paired regulation ('sif-style') interactions
## (can't think of anything other than a loop to do this )
flattenModVec <- function( modvec, betastring ) {

  tf.vec <- unlist(lapply(modvec,"[[","tfs"))
  targlong <- character()
  signlong <- character()
  for( m in modvec ){
    targlong <- c(targlong,rep(m$targ,length(m$tfs)))
    h <- (sign(m[[betastring]])+3)/2
    b <- c('+','-')
    signlong <- c(signlong,b[h])
  }
  
  returnMat <- cbind(targlong,tf.vec,signlong)
  colnames(returnMat) <- c("Target","TF","Sign")
  return(returnMat)
}

##
##  merge two modvecs into a single modvec
##
mergeModVecs <- function( modvec.A, modvec.B, debug=FALSE ){
  returnVec <- modvec.A
  counter <- 0 
  for ( mod in modvec.B ){
    counter <- counter + 1
    if ( debug ){
      cat("Looking at:",counter, "in modvec.B\n")
    }
    returnVec <- incrementModVec( mod, returnVec, debug=debug )
  }
  return ( returnVec )
}


##
##  Increment ModVec with a newmod
##  If model does not exist, we concatenate the modvec
##  If exists already, we merge in the model evidence
incrementModVec <- function( newmod, modvec, debug=FALSE ){
 
  returnVec <- list()

  new.targ <- newmod$targ
  new.tfs <- sort(newmod$tfs) ## sorted, for better comparisons 
  
  targ.vec <- unlist(lapply(modvec,"[[","targ"))  
  merge.counter <- 0 
    
  if ( !(new.targ %in% targ.vec )){ ## if new target, increment

    returnVec <- modvec
    returnVec[[(length(returnVec)+1)]] <- newmod

  } else { ## if target already there, need to work harder

    ## modvec indices for that target
    inds <- which(unlist(lapply(modvec,"[[","targ")) %in% new.targ)

    ## models which do not have this target, will be left as is
    modvec.non.targ <- modvec[-inds]

    ## models with this target
    ## we'll go through each, and decide whether to merge
    ## newmod into it. We'll build this incrementally, as modvec.temp
    ## Add the end of this proces, we increment by newmod
    ## if not merges have been performed
    modvec.targ <- modvec[inds]
    modvec.temp <- list()
    
    for ( mod  in modvec.targ ){
      tfs.in.mod <- sort(mod$tfs) 
      if ( !identical(new.tfs,tfs.in.mod) ) { ## if tfs not the same, do a simple copy
        modvec.temp[[length(modvec.temp)+1]] <- mod
      } else { ## if tfs same, create a merged model
        ## for now, let's just increment the evidence 
        mod.temp <- mod
        mod.temp$evidence <- c(mod.temp$evidence, newmod$evidence )
        modvec.temp[[length(modvec.temp)+1]] <- mod.temp
        merge.counter <- merge.counter + 1
      }
    }
    ## if we performed no merges, the new model should be appended
    ## to modvec.temp
    if ( merge.counter == 0 ){ ## 
      modvec.temp[[length(modvec.temp)+1]] <- newmod
    }

    ## finally, put the non.targ and targ objects together
    returnVec <- modvec.non.targ
    for ( mod in modvec.temp ){
      returnVec[[length(returnVec)+1]] <- mod
    }
    
  } ## end tasks for target overlap

  if ( debug ){
    cat("Performed ",merge.counter," merges\n")
  }
  
  return( returnVec )
}



##
## prefix: prefix for the .noa, .eda files
## linuxdir: linux representation of directory, full path where model figures will be found
## windir: windows representation of directory, full path where model figures will be found
## both are the same physical space
## full path must include slashes 

## if multimod is TRUE then we plot the multiple possible models for each target
## In this selection, we will NOT have access to the model plots at the nodes

createCytoscapeFilesMV <- function(modvec,prefix,outputdir,linuxdir=NULL,windir=NULL,osxdir=NULL,ppints=TRUE ){

  dir.create(outputdir) ## This is the directory where the output files of this routine will be located
  
  sifFile <- paste(outputdir,prefix,".sif",sep="") 
  weightsFile <- paste(outputdir,prefix,".eda",sep="") 

  ## remove any pre-existing files with these names
  file.remove(sifFile)
  file.remove(weightsFile)

  modelFigFileLinux <- paste(outputdir,prefix,"ModelFigsLinux.noa",sep="") ## model figure locations
  modelFigFileWindows <- paste(outputdir,prefix,"ModelFigsWindows.noa",sep="") ## model figure locations
  modelFigFileOSX <- paste(outputdir,prefix,"ModelFigsOSX.noa",sep="") ## model figure locations
  ## remove any pre-existing files with these names
  if(file.exists(modelFigFileLinux)){file.remove(modelFigFileLinux)}
  if(file.exists(modelFigFileOSX)){file.remove(modelFigFileOSX)}
  if(file.exists(modelFigFileWindows)){file.remove(modelFigFileWindows)}

  ##
  ## protein-protein interactions
  ##
  if ( ppints ){
    tfs.ps <- unique(sort(as.character(unlist(lapply(modvec,"[[","tfs"))))) ## tfs, as psois  
    index.matrix <- which(tf.dist.cn[cname.compare[tfs.ps],cname.compare[tfs.ps]]==1,arr.ind=TRUE) ## which of these interact directly
    if ( (length(index.matrix)!=0) & (length(index.matrix)!= 2) ){
      index.matrix <- index.matrix[(index.matrix[,1] > index.matrix[,2]),] ## eliminate duplicated interactions ( tf.dist.cn is symmetric )
    }
    if ( length(index.matrix)!=0 ){
      for ( i in 1:nrow(index.matrix) ) {
        tf1 <- tfs.ps[index.matrix[i,1]]
        tf2 <- tfs.ps[index.matrix[i,2]]
        outstring1 <- paste(np[tf1],"pp",np[tf2],sep=" ")
        cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
      }
    }
  }

  ##
  ## protein-DNA interactions
  ## 

  cat("weight\n",file=weightsFile)
  cat("TimeCourseLinuxPath (class=java.net.URL)\n",file=modelFigFileLinux)
  cat("TimeCourseWindowsPath (class=java.net.URL)\n",file=modelFigFileWindows)
  cat("TimeCourseOSXPath (class=java.net.URL)\n",file=modelFigFileOSX)
  
  for ( mod in modvec ) {
    regs <- mod$tfs
    targ <- mod$targ
    targ.name <- np[targ]
    modlabel <- paste(cc[c(targ,regs)],collapse="-")

    
    ## Write location of model figure
    outstringTC.linux <- paste(targ.name," = file://",linuxdir,modlabel,".png",sep="")
    outstringTC.windows <- paste(targ.name," = file:\\\\",windir,modlabel,".png",sep="")
    outstringTC.osx <- paste(targ.name," = file://",osxdir,modlabel,".png",sep="")
    cat(outstringTC.linux,"\n",sep="",file=modelFigFileLinux,append=TRUE)
    cat(outstringTC.osx,"\n",sep="",file=modelFigFileOSX,append=TRUE)
    cat(outstringTC.windows,"\n",sep="",file=modelFigFileWindows,append=TRUE)
    ## write interactions and weights
    if ( length(regs) > 1 ){
      for ( i in 1:length(regs) ){
        outstring1 <- paste(np[regs[i]],"Protein-DNA",targ.name,sep=" ")
        cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
        outstring1.1 <- paste(np[regs[i]],"(Protein-DNA)",targ.name,sep=" ")
        if ( "betas.ode" %in% names(mod) ){
          coefs <- mod$betas.ode
        } else {
          coefs <- mod$betas
        }
        outstring2 <- paste(outstring1.1,"=",coefs[i],sep="")
        cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
      }
    } else if ( length(regs)==1 ) {
      outstring1 <- paste(np[regs],"Protein-DNA",targ.name,sep=" ")
      cat(outstring1,"\n",sep="",file=sifFile,append=TRUE)
      outstring1.1 <- paste(np[regs],"(Protein-DNA)",targ.name,sep=" ")
      if ( "betas.ode" %in% names(mod) ){
        coefs <- mod$betas.ode
      } else {
        coefs <- mod$betas
      }
      outstring2 <- paste(outstring1.1,"=",coefs,sep="")
      cat(outstring2,"\n",sep="",file=weightsFile,append=TRUE)
    }
  } ## end loop over models

}

createBioTapestryFilesMV <- function(modvec,prefix,outputdir=NULL, regionvec, sixpanel=FALSE ){

  dir.create(outputdir)
  
  sifFile <- paste(outputdir,prefix,".sif",sep="") 
  csvFile <- paste(outputdir,prefix,".csv",sep="") 
  
  ## remove any pre-existing files with these names ( we will be using file 'append' )
  file.remove(sifFile)
  file.remove(csvFile)
  
  ## Preamble

  if ( !sixpanel ){
    h1 <- "model,root\n"
    h2 <- "model,TFMac,root\n"
    h3 <- "region,TFMac,Before 2 hours,E\n"
    h4 <- "region,TFMac,4-6 hours,M\n"
    h5 <- "region,TFMac,After 6 hours,L\n"
    cat(h1,h2,h3,h4,h5,sep="",file=csvFile,append=TRUE)
  } else {
    h1 <- "model,root\n"
    h2 <- "model,TFMac,root\n"
    h3 <- "region,TFMac,TFs Before 2 hours,E\n"
    h4 <- "region,TFMac,TFs 4-6 hours,M\n"
    h5 <- "region,TFMac,TFs After 6 hours,L\n"
    h6 <- "region,TFMac,Cytokines Before 2 hours,P\n"
    h7 <- "region,TFMac,Cytokines 4-6 hours,Q\n"
    h8 <- "region,TFMac,Cytokines After 6 hours,R\n"
    cat(h1,h2,h3,h4,h5,h6,h7,h8,sep="",file=csvFile,append=TRUE)
  }
  
  for ( mod in modvec ) {
    regs <- mod$tfs
    targ <- mod$targ
    if ( "betas.ode" %in% names(mod) ){
      coefs <- mod$betas.ode
    } else {
      coefs <- mod$betas
    }
    
    for ( i  in 1:length(regs) ){

      reg <- regs[i]
      if ( coefs[i] >= 0 ) {
        pn <- "positive"
      } else {
        pn <- "negative"
      }
      
      
      s11 <- paste("general","TFMac","gene",cc[reg],"gene",cc[targ],pn,regionvec[reg],regionvec[targ],sep=",")
      cat(s11,"\n",file=csvFile,append=TRUE)
      ##general,TFMac,gene,Gene 1,gene,Gene 2,positive,A,A

      outstring1.cname <- paste(cc[reg],"(pdna)",cc[targ],sep="\t")
      cat(outstring1.cname,"\n",sep="",file=sifFile,append=TRUE)

    }
  } ## end loop over models

}

predPlotWrapperModVec <- function( modvec, savedir="REGPREDS", t.index.min=1,t.index.max=11,zero.offset=FALSE ) {

  dir.create(savedir)

  for ( mod in modvec ){

    targ <- mod$targ
    tfs <- mod$tfs
    n.tfs <- length(tfs)
    modfields <- names(mod)
    offset <- mod$offsets

    modlabel <- paste(cc[c(targ,tfs)],collapse="-")
    imageFileName <- paste(savedir,"/",modlabel,".png",sep="")
    
    ## Determine if ODE fit has been performed and plot accordingly

    if ( "betas.ode" %in% modfields ){ ## ODE-fitted plots begin here

      ivars <- tfs
      betas.ode <- mod$betas.ode
      tau.ode <- mod$tau.ode
      parsfit <- c(betas.ode,1)/tau.ode       
      scale.fac <- max.intensity[targ]
      G0 <- lps.mat.max1.exp[targ,t.index.min]
      
      if ( n.tfs==1 ){
        Gproduction <- generateGproduction(generateTFfunctions(ivars,n.cands=1),n.cands=1)    
      }  else if ( n.tfs==2 ) {        
        Gproduction <- generateGproduction(generateTFfunctions(ivars))    
      } else {
        stop("Oops, not set up for other than 1 or 2 regulator ODE solutions")
      }
      
      wt.pred.array.times <- scale.fac*lsoda(G0,array.times,Gproduction,parsfit)[,2]
      ndiv=50.
      times.dense <- steppes(ndiv=ndiv)
      xvals.dense <- seq(0,10,1/ndiv)
      wt.pred.dense <- scale.fac*lsoda(G0,times.dense,Gproduction,parsfit)[,2]
      
      pred.out <- max.intensity[targ]*wt.pred.array.times ## needs units
      pred.out.dense <- max.intensity[targ]*wt.pred.dense ## needs units ?

      #png(imageFileName,height=600,width=900)      
      png(imageFileName,height=800,width=1100)      
      names(betas.ode) <- tfs
      #plotPrediction.ode.minimal(targ,betas.ode,tau.ode,wt.pred.array.times, xvals.dense,wt.pred.dense)
      plotPrediction.ode(targ,betas.ode,tau.ode,wt.pred.array.times, xvals.dense,wt.pred.dense)
      dev.off()

    } else { ## finite-diff style fits begin here

      ivars <- tfs
      betas <- mod$betas
      names(betas) <- ivars
      ## should this be one of the newer versions of this function?
      wt.predicted.fullscale <- calculateTimeEvolPrediction.3(targ,ivars,tau,s=1,boost.vec,zero.offset)
      png(imageFileName,height=800,width=1100)
      plotPrediction.3(targ,betas,wt.predicted.fullscale)
      dev.off()
      
    } ## end finite-diff ops

  } ## end loop over modvec
      
}

## Plot 'contact sheets' with predictions
## filename must include full path
gridPredPlotWrapperMV <- function(modvec,file="grid.ps",nx=1,ny=1,t.index.min=1,t.index.max=11,tau=600./log(2.), zero.offset=TRUE) {
  
  postscript(file,paper='letter')

  par(mfrow=c(nx,ny),mar=(c(3, 2, 2, 1)+0.1),mgp=c(1.75, 0.75, 0))
  par(mfrow=c(nx,ny),mar=(c(2, 2, 1, 1)+0.1),mgp=c(1.75, 0.75, 0))

  counter <- 0

  for ( mod in modvec ){
    
    targ <- mod$targ
    tfs <- mod$tfs
    n.tfs <- length(tfs)
    modfields <- names(mod)
    offset <- mod$offsets

    if ( "betas.ode" %in% modfields ){ ## ODE-fitted plots begin here

      ivars <- tfs
      betas.ode <- mod$betas.ode
      tau.ode <- mod$tau.ode
      parsfit <- c(betas.ode,1)/tau.ode       
      scale.fac <- max.intensity[targ]
      G0 <- lps.mat.max1.exp[targ,t.index.min]
      
      if ( n.tfs==1 ){
        Gproduction <- generateGproduction(generateTFfunctions(ivars,n.cands=1),n.cands=1)    
      }  else if ( n.tfs==2 ) {        
        Gproduction <- generateGproduction(generateTFfunctions(ivars))    
      } else {
        stop("Oops, not set up for other than 1 or 2 regulator ODE solutions")
      }
      
      wt.pred.array.times <- scale.fac*lsoda(G0,array.times,Gproduction,parsfit)[,2]
      ndiv=50.
      times.dense <- steppes(ndiv=ndiv)
      xvals.dense <- seq(0,10,1/ndiv)
      wt.pred.dense <- scale.fac*lsoda(G0,times.dense,Gproduction,parsfit)[,2]
      
      pred.out <- max.intensity[targ]*wt.pred.array.times ## needs units
      pred.out.dense <- max.intensity[targ]*wt.pred.dense ## needs units ?

      names(betas.ode) <- tfs
      plotPrediction.ode.hack(targ,betas.ode,tau.ode,wt.pred.array.times, xvals.dense,wt.pred.dense)

      counter <- counter + 1

    } else {

      ivars <- tfs
      betas <- mod$betas
      names(betas) <- ivars
      ## should this be one of the newer versions of this function?
      ## tau is an input variable to the function 
      wt.predicted.fullscale <- calculateTimeEvolPrediction.3(targ,ivars,tau,s=1,boost.vec,zero.offset)
      plotPrediction.hack(targ,betas,wt.predicted.fullscale)
      counter <- counter+1
    }
    
  }

  dev.off()
   
}



## return modvec that has a single model for each gene
## based on the 
bestGeneModelsMV <- function  ( modvec ){
  unique.targs <- unique(sort(as.character(unlist(lapply(modvec,"[[","targ")))))
  returnList <- list()
  for ( targ in unique.targs ){
    modvec.for.targ <- sliceModVecByTargs(targ, modvec)
    best.mod.for.targ <- bestModMV(modvec.for.targ)
    returnList[[length(returnList)+1]] <- best.mod.for.targ
  }
  return( returnList  )
}

##
## From ModVec, which single model has the best score
## For now, we'll use rmsd error for rankning
bestModMV <- function( modvec ){
  if ( length(modvec)== 1){
    return(modvec[[1]])
  } else {
    rmsds <- unlist(lapply(modvec,"[[","rmsd.ode"))
    ## we'll need to do something else for smods, whose score is labeled "rmsd"
    best.index <- which.min(rmsds)
    return( modvec[[best.index]] )
  }
}

## Get index of a model withing a modvec
getModelIndexMV <- function( targ, regs, modvec){
  targ <- as.character(targ)  ## strip off any names
  regs <- as.character(regs) ## strip off any names
  index.return <- NULL
  for ( i in 1:length(modvec) ){
    mod <- modvec[[i]]
    criterion.1 <- mod$targ==targ
    criterion.2 <- identical(sort(regs),sort(mod$tfs))
    if ( criterion.1 & criterion.2 ){
      index.return <- i
    }
  }
  return( index.return )
}


## regs is to match any regulators in modvec
## we return the modvec indices of those matches
## exact.match=TRUE require exact match of regulator set
getModelIndexByRegsMV <- function( regs, modvec,exact.match=FALSE){
  regs <- as.character(regs) ## strip off any names
  index.vec.return <- NULL
  for ( i in 1:length(modvec) ){
    mod <- modvec[[i]]
    if ( !exact.match ){
      criterion <- length(intersect(regs,mod$tfs) > 0 )
    } else {
      criterion <- identical(sort(regs),sort(mod$tfs))
    }
    if ( criterion  ){
      index.vec.return <- c(index.vec.return,i)
    }
  }
  return( index.vec.return )
}

