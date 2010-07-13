
##
## Initialize tally matrices
## Retaining order of pair
##

## Singles

nm <- length(mouseMatrices)
m.tally <- integer(length=nm)
names(m.tally) <- mouseMatrices

nmf <- length(familyNames)
mf.tally <- integer(length=nmf)
names(mf.tally) <- familyNames

### Pairs, retaining order

m.pair.tally <- matrix(0,nrow = nm,ncol = nm)
colnames(m.pair.tally) <- mouseMatrices
rownames(m.pair.tally) <- mouseMatrices

mf.pair.tally <- matrix(0,nrow = nmf,ncol = nmf)
colnames(mf.pair.tally) <- familyNames
rownames(mf.pair.tally) <- familyNames

##
## Pairs, disregarding order
##

mdo.pair.tally <- matrix(0,nrow = nm,ncol = nm)
colnames(mdo.pair.tally) <- mouseMatrices
rownames(mdo.pair.tally) <- mouseMatrices

mdof.pair.tally <- matrix(0,nrow = nmf,ncol = nmf)
colnames(mdof.pair.tally) <- familyNames
rownames(mdof.pair.tally) <- familyNames

#define input arguments
maxsep <- 100; #the maximum allowable distance between hits to be retained
minsep <- 1; #don't retain pairs if they are less than this distance apart
pdist <- 1000; #don't retain hits if their average distance is more than this from the TSS

NG <- length(TRE.OUT); #the number of genes with data

counter <- 0 
for (entrezID in entrezIDs ){
##for(ss in 1:3){#loop over all genes

  counter <- counter + 1
  ensmusgs <- ensemblIDsOfEntrezIDs[[entrezID]]

  m.tally.logical <- logical(length=nm)
  names(m.tally.logical) <- mouseMatrices
  mf.tally.logical <- logical(length=nmf)
  names(mf.tally.logical) <- familyNames

  m.pair.logical <- matrix(FALSE,nrow = nm,ncol = nm)
  rownames(m.pair.logical) <- mouseMatrices
  colnames(m.pair.logical) <- mouseMatrices
  
  mf.pair.logical <- matrix(FALSE,nrow = nmf,ncol = nmf)
  rownames(mf.pair.logical) <- familyNames
  colnames(mf.pair.logical) <- familyNames
  
  for ( ensmusg in ensmusgs ){
  
    gene.info <- TRE.OUT[[ensmusg]]$seq.info
    tre <- TRE.OUT[[ensmusg]]$tre.hits; #select the particular tf binding site info for a specific gen
    
    keep.indices <-  which(tre$score > matrixThreshold[tre$tre]) ## which hits are above threshold
    tre <- tre[keep.indices,] ## keep only those
 
    NT <- nrow(tre); #the number of hits for this gene
    if(length(NT)!=0 ){ #check to make sure there actually are hits before processing further
      if(NT > 0  ){ #check to make sure there actually are hits before processing further
        
        tre[,'tre.original'] <- tre[,'tre']; #save the original names for later
        tre[,'tre'] <- apply(as.data.frame(tre[,'tre']),1,split1); #just use the prefix for the matrix
        tre.P <- tre[,'tre']; #we'll store the individual occurrences along with the pairs
        
        found.mats <- unique(sort(tre[,'tre.original']))
        found.families <- unique(sort(as.character(familyMap[found.mats])))
        
        m.tally.logical <- m.tally.logical | ( mouseMatrices %in% found.mats )## logical will flip from FALSE to TRUE the first time, then remain on
        mf.tally.logical <- mf.tally.logical | ( familyNames %in% found.families )
        
        ##we will start closest to the TSS and work backwards
        
        for(tt in 1:(NT-1)){ #loop through all hits and find close neighbors
          tre.o <- tre[tt,]; #the information for the particular TRE
          start.length <- matrixLength[tre.o[,"tre.original"]];   
          tre.p <- tre[(tt+1):NT,]; #looking backwards at all other TREs
          tre.p$tdist <- -tre.p[,'dist'] + tre.o$dist; #compute the distances between hits
          tre.p$adist <- abs(round(0.5*(tre.p[,'dist'] + tre.o$dist))); #compute the average distance from TSS for hits
          tre.p.select <- which((tre.p$tdist<=maxsep & tre.p$tdist>=(minsep+start.length)) & tre.p$adist <= pdist);
          tre.p <- tre.p[tre.p.select,]; #filter out that are too far away from each other, too close to each other, or too far from TSS
          np <- nrow(tre.p); #the number of remaining hits within mdist from the first one
          if(np>0){

            moi <- tre.o[,"tre.original"]
            partners <-  unique(sort(tre.p$tre.original))
            m.pair.logical[moi,] <-  m.pair.logical[moi,] | (mouseMatrices %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
         
            moif <- as.character(familyMap[moi])
            partners <- unique(sort(as.character(familyMap[tre.p$tre.original])))    
            mf.pair.logical[moif,] <-  mf.pair.logical[moif,] | (familyNames %in% partners) ## logical will flip from FALSE to TRUE the first time, then remain on
           };  #close test on np
        } ## close loop over NT

      } ## close test for NT > 1
    } ## close test for nonzero NT

  } ## close loop over ensembl IDs for the specific entrez ID

  m.tally <- m.tally + m.tally.logical # increment single tally, utilizing integer + logical = integer
  mf.tally <- mf.tally + mf.tally.logical
  
  m.pair.tally <- m.pair.tally + m.pair.logical ## increment pair tally, utilizing integer + logical = integer
  mdo.pair.tally <- mdo.pair.tally + (m.pair.logical | t(m.pair.logical) ) ## disregard order 
  
  mf.pair.tally <- mf.pair.tally + mf.pair.logical ## increment pair tally, utilizing integer + logical = integer
  mdof.pair.tally <- mdof.pair.tally + ( mf.pair.logical | t(mf.pair.logical) ) ## disregard order
  
  write(c(counter,counter/length(entrezIDs)),file='status.txt'); #report status

}; rm(counter);  #close loop over entrez IDs

save(m.tally,m.pair.tally,mf.tally,mdo.pair.tally,mf.pair.tally,mdof.pair.tally,file="hitTally.allMouse.RData",compress=TRUE)
