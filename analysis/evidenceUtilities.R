

fastaNbrs <- function ( psoi , outFile ) {
  fcollection <- character()
  fcollection <- c(fcollection,fastaFileName(psoi))
  friends <- nbrs[[psoi]]
  for ( friend in friends ){
    fcollection <- c(fcollection,fastaFileName(friend))
  }  
  fcollection
}

writeFastaNbrs  <- function ( psoi , outFile ) {
  fcollection <- character()
  fcollection <- c(fcollection,fastaFileName(psoi))
  friends <- nbrs[[psoi]]
  for ( friend in friends ){
    fcollection <- c(fcollection,fastaFileName(friend))
  }  
  write.table(fcollection,file=outFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

writeNbrsPrefixes <- function ( psoi , outFile ) {
  fcollection <- character()
  fcollection <- c(fcollection,fileNamePrefix(psoi))
  friends <- nbrs[[psoi]]
  for ( friend in friends ){
    fcollection <- c(fcollection,fileNamePrefix(friend))
  }  
  write.table(fcollection,file=outFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
}



## Assumes that nbrs is available
fileNamePrefix <- function ( psoi ){
  eid <- ncbiID[psoi]
  ensid <- ensemblIDsOfEntrezIDs[[ncbiID[psoi]]]
  filestring <- paste(cname.compare[psoi],eid,ensid,sep="-")
  filestring
}



## Assumes that nbrs is available
fastaFileName <- function ( psoi ){
  eid <- ncbiID[psoi]
  ensid <- ensemblIDsOfEntrezIDs[[ncbiID[psoi]]]
  filestring <- paste(cname.compare[psoi],eid,ensid,sep="-")
  filename <- paste(filestring,".fasta",sep="")
  filename
}
  

evidenceForMods <- function( mods, n.cands=1){

  returnList <- list()
  psois <- names(mods)
  for ( psoi in psois ){

    mod <- mods[[psoi]]

    returnList[[psoi]]$tfs <- mod$tfs
    returnList[[psoi]]$tfs.cname <- mod$tfs.cname

    ensids <- ensemblIDsOfEntrezIDs[[ncbiID[psoi]]]
    n.ensids <- length(ensids)
        
    cat("Looking at ", psoi,cname.compare[psoi],ensids,"\n")
    
    ## may need to check that this is unique

     if ( n.cands==1 ) {

      nmods <- length(mod$tfs)

      for ( n in 1:nmods ){

        tf.cname <- mod$tfs.cname[n]
        
        msingles <- metamsingles[[ensid]]
        pvals <- metapvals.singles[[ensid]]
        ms <- metams.singles[[ensid]]

        ## vectors, labeled by family name, indicating whether as specific TF is found
        
        logvec <- unlist(lapply(lapply(fam2tf.tte.gname,intersect,tf.cname),length))
        inds <- as.integer(which(logvec))
        pvs <- pvals[inds]
        ems <- ms[inds]
        
        ## Keep best p-value
        if ( length(inds) > 1 ){
          ibest <- which.min(pvs)
          pvs <- pvs[ibest]
          ems <- ems[ibest]
        }
        
        ## May need a list object incrementation here
        returnList[[psoi]]$pval <- c(returnList[[psoi]]$pval,pvs)
        returnList[[psoi]]$m <- c(returnList[[psoi]]$m,ems)
        
      }      
    } else if ( n.cands==2 ) {

      nmods <- nrow(mod$tfs)
      for ( n in 1:nmods ){ ## Consider each model

        tf1.cname <- mod$tfs.cname[n,1]
        tf2.cname <- mod$tfs.cname[n,2]          
        
        ## Extra machinery for the rare mutliple ensid for each entrez ID
 
        tempStore <- matrix(nrow=n.ensids,ncol=5,dimnames=list(ensids,c("pvs","ems","mexps","mp1","mp2") )) ## Will store results temporarily for each ensid
        for ( ensid in ensids ){
          
          mpairs <- metampairs[[ensid]]
          pvals <- metapvals[[ensid]]
          ms <- metams[[ensid]]
          mexpects <- metaexpected[[ensid]]
          
          ## vectors, labeled by family name, indicating whether as specific TF is found
          tf1.indicator <- unlist(lapply(lapply(fam2tf.tte.gname,intersect,tf1.cname),length))
          tf2.indicator <- unlist(lapply(lapply(fam2tf.tte.gname,intersect,tf2.cname),length))
          
          logvec.1 <- tf1.indicator[mpairs[,1]] & tf2.indicator[mpairs[,2]]
          logvec.2 <- tf1.indicator[mpairs[,2]] & tf2.indicator[mpairs[,1]]
          logvec <- logvec.1 | logvec.2
          
          inds <- as.integer(which(logvec))
          
          mps <- mpairs[inds,]
          pvs <- pvals[inds]
          ems <- ms[inds]
          mexps <- mexpects[inds]

          if ( length(inds)==0){
            tempStore[ensid,] <- c(1.,0.,0.,0.,0.) ## This p-value will lose
          } else if ( length(inds) == 1 ) {  
            tempStore[ensid,] <- c(pvs,ems,mexps,mps) 
          } else if ( length(inds) > 1 ){ ## Keep best p-value
            ibest <- which.min(pvs)
            pvs <- pvs[ibest]
            ems <- ems[ibest]
            mexps <- mexps[ibest]
            mps <- mps[ibest,]
            tempStore[ensid,] <- c(pvs,ems,mexps,mps) 
          }         
          
        } ## end loop over ensids

        ## Determine which ensid provides the sites with the lowest p-value

        index.best.ensid <- which.min(as.numeric(tempStore[,"pvs"]))
        ibe <- index.best.ensid 
        ensid <- rownames(tempStore)[ibe] ## consider storing this
        mps <- as.character(tempStore[ibe,c("mp1","mp2")])
        pval <- as.numeric(tempStore[ibe,"pvs"])
        ems <- as.numeric(tempStore[ibe,"ems"])
        mexps <- as.numeric(tempStore[ibe,"mexps"])
        
        ## Store
        ## returnList[[psoi]]$n <- metan[[ensid]] Need to understand what this is - probably the total N for Fisher test
        returnList[[psoi]]$mats <- rbind(returnList[[psoi]]$mats,drop(mps))
        returnList[[psoi]]$pval <- c(returnList[[psoi]]$pval,pvs)
        returnList[[psoi]]$mhits <- c(returnList[[psoi]]$mhits,ems)
        returnList[[psoi]]$mexpected <- c(returnList[[psoi]]$mexpected,mexps)
        

          
        } ## end loop over models

    } ## end if n.cand==2
  } ## end loop over psois
  
  returnList
  
}
  
      
