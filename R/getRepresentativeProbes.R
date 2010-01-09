
##   Jan 18, 2005
##  This function has been sped up with tapply
##

source("./utilitiesExpression.R")

mu.cutoff <- 100
##mu.cutoff <- 0
lambda.cutoff <- 57.2
#lambda.cutoff <- 150

imax <- 11
lps.full.ps.sig <- rownames(sigSlice(lambda.cutoff,lps.mus[,1:imax],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[,1:imax]<mu.cutoff,1,sum)==imax))
lps.full.ps.sig <- setdiff(lps.full.ps.sig,low.expressors)

## The amount of change is quantified as the sum of lambdas
change.value <- apply(lps.lambdas,1,sum)
nwm <- function(vals){names(which.max(vals))}
repProbes.np <- tapply(change.value,as.factor(np),nwm) 
repProbes.cname <- tapply(change.value,as.factor(cname.compare),nwm)
repProbes.ncbiID <- tapply(change.value,as.factor(ncbiID),nwm) 

## Hopefully this is obsolete ( April 2008 ) 
## representativeProbes <- repProbes.np


## In the above, some 'non-expressed' probesets were selected as repProbes
## derived from lps.full.ps.sig

change.value.diffexp <- apply(lps.lambdas[lps.full.ps.sig,],1,sum)
repProbes.np.diffexp <- tapply(change.value[lps.full.ps.sig],as.factor(np[lps.full.ps.sig]),nwm) 
repProbes.cname.diffexp <- tapply(change.value[lps.full.ps.sig],as.factor(cname.compare[lps.full.ps.sig]),nwm)
repProbes.ncbiID.diffexp <- tapply(change.value[lps.full.ps.sig],as.factor(ncbiID[lps.full.ps.sig]),nwm) 

##
## For the diffexp ones, 'replace' the repProbe mappings
##

repProbes.cname[names(repProbes.cname.diffexp)] <- repProbes.cname.diffexp
repProbes.np[names(repProbes.np.diffexp)] <- repProbes.np.diffexp
repProbes.ncbiID[names(repProbes.ncbiID.diffexp)] <- repProbes.ncbiID.diffexp

## Remove "multiple" annotator if present
if ( "multiple" %in% names(repProbes.cname) ){
  repProbes.cname <- repProbes.cname[-which(names(repProbes.cname)=="multiple")]
}
if ( "multiple" %in% names(repProbes.np) ){
  repProbes.np <- repProbes.np[-which(names(repProbes.np)=="multiple")]
}
if ( "multiple" %in% names(repProbes.ncbiID) ){
  repProbes.ncbiID <- repProbes.ncbiID[-which(names(repProbes.ncbiID)=="multiple")]
}

## Remove "---" annotator if present
if ( "---" %in% names(repProbes.cname) ){
  repProbes.cname <- repProbes.cname[-which(names(repProbes.cname)=="---")]
}
if ( "---" %in% names(repProbes.np) ){
  repProbes.np <- repProbes.np[-which(names(repProbes.np)=="---")]
}
if ( "---" %in% names(repProbes.ncbiID) ){
  repProbes.ncbiID <- repProbes.ncbiID[-which(names(repProbes.ncbiID)=="---")]
}



## Have checked that that all values are the same as from the loop algorithm below

##
## Get Representative Probes 
##
## Do not rerun unless needed. Takes forever!!!

##uniqueNPs <- unique(sort(as.character(np)))
##n.uniqueNPs <- length(uniqueNPs)

##representativeProbes <- character(length=n.uniqueNPs)
##names(representativeProbes) <-  uniqueNPs

##
## Analysis of duplcated vs. non-duplicted was not needed, in the end
##

## rep.indices <- as.numeric(as.factor(np)) 

##non.duplicated.nps <- np[!(duplicated(np))]
## note that the unique function does not apply here, as it returns duplicated NPs

##representativeProbes[non.duplicated.nps] <- names(non.duplicated.nps)
##dups <- np[duplicated(np)]
## unbelievable, doesn't seem to be possible to find *which* elements are unique and which are duplicates

##rep.count <- tapply(np,as.factor(np),length)
##non.duplicated.nps <- names(which(rep.count==1))





##for ( uniqueNP in uniqueNPs ){
##  members <-probesetFromNP(uniqueNP)
  ## The AFFX might appear!
  ##if ( length(elementsWithSubString("AFFX",members)) !=0  ){
  ##  toRemove <- elementsWithSubString("AFFX",members)
  ##  members <- setdiff(members,toRemove)
  ##}

  ##if ( length(members) == 0 ){
  ##  representativeProbes[uniqueNP] <- NA
  ##} else if ( length(members) == 1 ){
  ##  representativeProbes[uniqueNP] <- members
  ##} else {
    ##reducedMat <- lps.lambdas[members,]
    ##representativeProbes[uniqueNP] <- getMaxChangingProbeSet(reducedMat)
  ##}
##}


if ( FALSE ) { 

### Can we convert the whole set
tosser.indices <- which(is.na(representativeProbes)) ## Looks like the AFFX set alone
hh <- representativeProbes[-tosser.indices]


jj.ratios <- lps.ratios[hh,]
rownames(jj.ratios) <- npFromProbeSet(rownames(jj.ratios))

write.table(jj.ratios,file="lps.repProbe.ratios",quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t')

jj.lambdas <- lps.lambdas[hh,]
rownames(jj.lambdas) <- npFromProbeSet(rownames(jj.lambdas))

write.table(jj.lambdas,file="lps.repProbe.lambdas",quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t')

jj.mus <- lps.mus[hh,] 
rownames(jj.mus) <- npFromProbeSet(rownames(jj.mus))

write.table(jj.mus,file="lps.repProbe.mus",quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t')


}
