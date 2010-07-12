
##
## December 2008
##
## This code produces
## featureMatrix.mf.RData
## containing features per matrix family
## relies on preamble find.singles_preamble.R

efm <-  unique(as.character(familyMap[expressedMats]))

## Singles

nm <- length(mouseMatrices)
m.tally <- integer(length=nm)
names(m.tally) <- mouseMatrices

nmf <- length(familyNames)
mf.tally <- integer(length=nmf)
names(mf.tally) <- familyNames

counter <- 0
ptm <- proc.time()

##entrezIDs <- entrezIDs.car

efmats <- sort(unique(familyMap[expressedMats]))

featureMatrix.mf <- logical()
eids.keepers <- character()
  
for (entrezID in entrezIDs ){

  counter <- counter + 1
  ensmusgs <- ensemblIDsOfEntrezIDs[[entrezID]]

  m.tally.logical <- logical(length=nm)
  names(m.tally.logical) <- mouseMatrices
  mf.tally.logical <- logical(length=nmf)
  names(mf.tally.logical) <- familyNames

  for ( ensmusg in ensmusgs ){
    
    gene.info <- TRE.OUT[[ensmusg]]$seq.info
    tre <- TRE.OUT[[ensmusg]]$tre.hits; #select the particular tf binding site info for a specific gen
    
    keep.indices <-  which(tre$score > matrixThreshold[tre$tre]) ## which hits are above threshold
    tre <- tre[keep.indices,] ## keep only those
    ## AUGUST 2008. Restrict to expressed matrices
    keep.indices <- which(tre$tre %in% expressedMats )
    tre <- tre[keep.indices,] ## keep only those
    
    NT <- nrow(tre); #the number of hits for this gene
    if(length(NT)!=0 ){ #check to make sure there actually are hits before processing further
      if(NT > 0  ){ #check to make sure there actually are hits before processing further
        
        tre[,'tre.original'] <- tre[,'tre']; #save the original names for later
        tre[,'tre'] <- apply(as.data.frame(tre[,'tre']),1,split1); #just use the prefix for the matrix
        ##
        ## July 2008. The above seems like the wrong thing to. On the other hand, this column is not used!
        ##
        tre.P <- tre[,'tre']; #we'll store the individual occurrences along with the pairs
        
        found.mats <- unique(sort(tre[,'tre.original']))
        found.families <- unique(sort(as.character(familyMap[found.mats])))

        mf.logical <- efmats %in% found.families
        ##mf.logical <- familyNames %in% found.families

        featureMatrix.mf <- rbind(featureMatrix.mf,mf.logical)

        ##m.tally.logical <- m.tally.logical | ( mouseMatrices %in% found.mats )## logical will flip from FALSE to TRUE the first time, then remain on
        ##mf.tally.logical <- mf.tally.logical | ( familyNames %in% found.families )
      }    
    }
  } ## close loop over ensembl IDs for the specific entrez ID
 

  ##mdof[[ensmusg]] <- mm
    
  ##m.tally <- m.tally + m.tally.logical # increment single tally, utilizing integer + logical = integer
  ##mf.tally <- mf.tally + mf.tally.logical
  
  write(c(counter,counter/length(entrezIDs)),file='status.txt'); #report status
  eids.keepers <- c(eids.keepers,entrezID)
  
}; rm(counter);  #close loop over entrez IDs

rownames(featureMatrix.mf) <- entrezIDs
colnames(featureMatrix.mf) <- efmats

featureMatrix.mf <- t(featureMatrix.mf)

save(featureMatrix.mf,file=paste(seq.dir,"featureMatrix.mf.RData",sep="/"))

##save(m.tally,m.pair.tally,mf.tally,mdo.pair.tally,mf.pair.tally,mdof.pair.tally,file="hitTally.RData",compress=TRUE)
proc.time() - ptm
