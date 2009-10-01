
## n.cands=1 we consider single TF candidates
## n.cands=2 we consider double TF candidates

## targsAndCands: targs are psois, TFs are cnames

##tau <- 600/log(2)
##t.index.min <- 1
##t.index.max <- 11

calcModels <- function( targsAndCands, t.index.min, t.index.max, tau, boost.vec, mus.mat.maxnormed, n.cands=2, zero.offset=FALSE , tau.estimate=FALSE){

  returnList <- list()
  
  targsToPredict <- names( targsAndCands )

  counter <- 0
  for ( targ in targsToPredict ){

    returnList[[targ]] <- list()
    returnList[[targ]]$tfs <- NULL ## This entry will be a string matrix or a vector
    returnList[[targ]]$tfs.cname <- NULL  ## This entry will be a string matrix or a vector
    returnList[[targ]]$betas <- NULL ## This entry will be a matrix or a vector
    returnList[[targ]]$offsets <- NULL ## This entry will be a numeric vector.
    returnList[[targ]]$rmsds <- NULL ## This entry will be a numeric vector.
    returnList[[targ]]$cors <- NULL ## This entry will be a numeric vector.
    ## These declarations seem not to be necessary!
    
    counter <- counter+1 
    write(c(targ,counter,counter/length(targsToPredict)),file='status.txt'); #report status

    if ( n.cands==1 ){

      tfs <- targsAndCands[[targ]]
      n.tfs <- length(tfs) 
      for ( i in 1:n.tfs ) {
        ivars <- as.character(repProbes.cname[tfs[i]])
        if ( !(targ %in% ivars ) ){

          if ( tau.estimate ){
            tau <- estimateBestTau2(targ,ivars,boost.vec,lps.mat.max1)
          }

          dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,mus.mat.maxnormed)
          ris <- regInfSelectD2F(dataToFit,s=1,zero.offset)
          coefs <- as.numeric(ris[-1])
          offset <- as.numeric(ris[1])
          outpred <- calculateOutpred(drop(dataToFit$inputs),offset,coefs)
          cor <- calcCOR(outpred,dataToFit$outputs)
          rmsd <- calcRMSD(outpred,dataToFit$outputs)
          
          start.val <- lps.mat.max1[targ,t.index.min]
          scale.fac <- max.intensity[targ] ## The RHS is unfortunately a global variable
          wt.predicted.fullscale <- calculateTimeEvolPrediction.6(t.index.min,t.index.max,tau,outpred,start.val,scale.fac)
          cor.unrolled <- calcCOR(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max])
          rmsd.unrolled <- calcRMSD(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max])

          returnList[[targ]]$tfs <- c(returnList[[targ]]$tfs, as.character(ivars) )
          returnList[[targ]]$tfs.cname <- c(returnList[[targ]]$tfs.cname, as.character(cname.compare[ivars]) )
          returnList[[targ]]$betas <- c(returnList[[targ]]$betas, coefs )
          returnList[[targ]]$offsets <- c(returnList[[targ]]$offsets, offset )
          if ( tau.estimate ){
            returnList[[targ]]$taus <- c(returnList[[targ]]$taus, tau )
          }
          returnList[[targ]]$rmsds <- c(returnList[[targ]]$rmsds, rmsd)
          returnList[[targ]]$cors <- c(returnList[[targ]]$cors, cor)
          returnList[[targ]]$rmsds.unrolled <- c(returnList[[targ]]$rmsds.unrolled, rmsd.unrolled)
          returnList[[targ]]$cors.unrolled <- c(returnList[[targ]]$cors.unrolled, cor.unrolled)
        }
      } ## complete loop over pairs

    }## close cands=1 section
    
    if ( n.cands==2 ){

      tf.pairs <- targsAndCands[[targ]]
      n.pairs <- nrow(tf.pairs) ## This does seem to work for '1 row' case
      for ( i in 1:n.pairs ) {
        ivars <- repProbes.cname[tf.pairs[i,]]
        if ( !(targ %in% ivars ) ){          
          if ( tau.estimate ){
            tau <- estimateBestTau2(targ,ivars,boost.vec,lps.mat.max1)
          }
          dataToFit <- constructDataToFit(t.index.min,t.index.max,tau,targ,ivars,boost.vec,mus.mat.maxnormed)
          ris <- regInfSelectD2F(dataToFit,s=1,zero.offset)
          coefs <- as.numeric(ris[-1])
          offset <- as.numeric(ris[1])
          outpred <- calculateOutpred(dataToFit$inputs,offset,coefs)
          cor <- calcCOR(outpred,dataToFit$outputs)
          rmsd <- calcRMSD(outpred,dataToFit$outputs)
          
          start.val <- lps.mat.max1[targ,t.index.min]
          scale.fac <- max.intensity[targ] ## The RHS is unfortunately a global variable
          wt.predicted.fullscale <- calculateTimeEvolPrediction.6(t.index.min,t.index.max,tau,outpred,start.val,scale.fac)
          cor.unrolled <- calcCOR(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max])
          rmsd.unrolled <- calcRMSD(wt.predicted.fullscale, lps.mus[targ,(t.index.min+1):t.index.max])
          
          returnList[[targ]]$tfs <- rbind(returnList[[targ]]$tfs, ivars )
          returnList[[targ]]$tfs.cname <- rbind(returnList[[targ]]$tfs.cname, as.character(cname.compare[ivars]) )
          returnList[[targ]]$betas <- rbind(returnList[[targ]]$betas, coefs )
          returnList[[targ]]$offsets <- c(returnList[[targ]]$offsets, offset )
          if ( tau.estimate ){
            returnList[[targ]]$taus <- c(returnList[[targ]]$taus, tau )
          }
          returnList[[targ]]$rmsds <- c(returnList[[targ]]$rmsds, rmsd)
          returnList[[targ]]$cors <- c(returnList[[targ]]$cors, cor)
          returnList[[targ]]$rmsds.unrolled <- c(returnList[[targ]]$rmsds.unrolled, rmsd.unrolled)
          returnList[[targ]]$cors.unrolled <- c(returnList[[targ]]$cors.unrolled, cor.unrolled)
        }
      } ## complete loop over pairs

      ## remove spurious names
      rownames(returnList[[targ]]$betas) <- NULL
      colnames(returnList[[targ]]$betas) <- NULL
      rownames(returnList[[targ]]$tfs) <- NULL
      colnames(returnList[[targ]]$tfs) <- NULL
      rownames(returnList[[targ]]$tfs.cname) <- NULL
      colnames(returnList[[targ]]$tfs.cname) <- NULL 

    }## close cands=2 section
  } ## close loop over cands to predict
  
  returnList

}

