
## April 2008: Need to do a round of debugging for matrix families code


## Expressed Transfact TFs
##
## After correcting the mappings (in badMaps.tsv )
##
annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")

load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"all.ps.list.objects.RData",sep="/"))
load(paste(exp.dir,"all.mus.objects.RData",sep="/"))

gg <- read.table(paste(annot.dir,"factorsWithMatrices",sep="/"))
transfac.ncbiID <- as.character(gg$V2)
transfac.gname <- as.character(gg$V1)

transfac.gname <- transfac.gname[transfac.gname %in% cname.compare]
## Pax2 and Pou4f1 have no probes (Steve)
## VT: The latter *may* have probes ( 1429667_at 1429668_at  ) but there's no signal anyhow

## August 2006: Pou4f1 seems to be 'back' according to Affy,and to Adapt
## Ensembl Martview (still)list 1429668_at as well outside, and 1429667_at just touching the UTR
## Perhaps discrepancy is because the two probes are associated with NM_011143, maybe not in Ensembl. USCS agrees.

transfac.ps <- repProbes.cname[transfac.gname]

## Differentially expressed TFs

ps.exp <- lps.full.ps.sig 
transfac.tfs.expressed <- ps.exp[ps.exp %in% transfac.ps]
transfac.tfs.expressed <- sort(transfac.tfs.expressed) ## mostly for unified output order


## Factors expressed, but not differentially expressed
mu.cutoff <- 100
transfac.tfs.no.intensity <- names(which(apply(lps.mus[transfac.ps,]  > mu.cutoff,1,sum)==0 ))
transfac.tfs.constitutive.ps <- setdiff( setdiff(transfac.ps,transfac.tfs.no.intensity), transfac.tfs.expressed )

##Looking at these, some were just under the diff exp threshold, but do appear to change, eg Bach1


##
## Mappings to matrices 
##
gg <-read.table(paste(annot.dir,"matrixLengths",sep="/"),as.is=TRUE);
matrixLength <- gg$V2; names(matrixLength) <- gg$V1; rm(gg);

## Mapping of Matrix to TFs ( one-to-many )
t1 <- read.table(paste(annot.dir,"MatrixToMouseGeneIDuniqueGeneID.tsv",sep="/"),as.is=TRUE);
allMats <- as.character(read.table(paste(annot.dir,"MouseMatrices",sep="/"),as.is=TRUE)$V1)
mouseMatrices <- allMats

matVec <- substring(t1$V2,3) ## chops off the 'V$' prefix used by Transfac
tfVecGeneID <- as.character(t1$V3)
tfVecGName <- t1$V4

## The 269 TF names agree with transfac.gname from inititialTFs.R
##identical(sort(transfac.gname),unique(sort(tfVecGName)))
##identical(sort(transfac.ncbiID),unique(sort(tfVecGeneID)))
## transfac.ps is the corrected bs list for these 269, ordered alphabetically by TF gene name

## These now differ in Pax2 and Pou4f1, previously associated with 14445000_at and 14296667_at respectively
## drawn into doubt by steves 6/22/06 analysis

tf.not.on.array.inds <- which(tfVecGName %in% c("Pax2") )
##tf.not.on.array.inds <- which(tfVecGName %in% c("Pax2","Pou4f1") )
matVec <- matVec[-tf.not.on.array.inds]
tfVecGeneID <- tfVecGeneID[-tf.not.on.array.inds]
tfVecGName <- tfVecGName[-tf.not.on.array.inds]
### Do not repeat the above three commands without initializing 

## Contruct list to map PWMs to targets 
pwm2tf.ps <- list()
pwm2tf.gname <- list()
pwm2tf.ncbiID <- list()
pwm2tf.tte.ps <- list() ## mapping, for just the expressed TFs
pwm2tf.tte.gname <- list() ## mapping, for just the expressed TFs

for ( mat in allMats ){
  tf.inds <- which(matVec==mat)
  tfs.ncbiID <- tfVecGeneID[tf.inds]
  tfs.gname <- tfVecGName[tf.inds]
  tfs.ps <- as.character(repProbes.cname[tfs.gname])
  ##tfs.ps <-as.character(transfac.ps[tfs.gname])
  pwm2tf.ps[[mat]] <-  tfs.ps
  pwm2tf.gname[[mat]] <-  tfs.gname
  pwm2tf.ncbiID[[mat]] <-  tfs.ncbiID
  pwm2tf.tte.ps[[mat]] <-  tfs.ps[tfs.ps %in% transfac.tfs.expressed]
  pwm2tf.tte.gname[[mat]] <- cname.compare[  pwm2tf.tte.ps[[mat]] ]
}

## Have checked that the following is equal to tranfac.tfs.expressed
##transfac.ps.exp <- unique(sort(as.character(unlist(pwm2tf.tte.ps))))

## Contruct list to map targets to PWMs

tfps2pwm <- list() ## reference by ps
for ( psoi in transfac.ps ){
  collection <- character()
  for ( mat in allMats ){
    if (psoi %in% pwm2tf.ps[[mat]] ){
      collection <- c(collection,mat)
    }
  collection <- unique(sort(collection))
  tfps2pwm[[psoi]] <- collection
  }
}


## PWMs for which the TF is expressed
expressedMats <- vector()
for ( psoi in transfac.tfs.expressed ){
  expressedMats <- c(expressedMats,tfps2pwm[[psoi]])
}
expressedMats <- unique(sort(expressedMats))

constitutiveMats <- setdiff ( unique(sort(as.character(unlist(tfps2pwm[transfac.tfs.constitutive.ps])))), expressedMats )

expressedMatsPrefix <- unique(sort(as.character(sapply(expressedMats,split1))))

## Derive mappings for families from the above

## load familyMap, a mapping from vector of matrix names to family strings
load(paste(annot.dir,"familyMapVT3.RData",sep="/"))
familyMap <- familyMapVT3; rm(familyMapVT3);
familyNames <- unique(sort(as.character(familyMap)))
allFams <- unique(sort(as.character(familyMap)))

## April 2008. This one is now reporting partial match warnings
fam2tf.tte.gname <- list() ## mapping, for just the expressed TFs
for ( mat in expressedMats ){
  family <- as.character(familyMap[mat])
  fam2tf.tte.gname[[family]] <- c(fam2tf.tte.gname[[family]],pwm2tf.tte.gname[[mat]])
}
for ( family in  names(fam2tf.tte.gname) ){
  fam2tf.tte.gname[[family]] <- unique(sort(fam2tf.tte.gname[[family]]))
}

fam2tf.gname <- list() ## mapping, for just the expressed TFs
for ( mat in mouseMatrices ){
  family <- as.character(familyMap[mat])
  fam2tf.gname[[family]] <- c(fam2tf.gname[[family]],pwm2tf.gname[[mat]])
}
for ( family in  names(fam2tf.gname) ){
  fam2tf.gname[[family]] <- unique(sort(fam2tf.gname[[family]]))
}


tte <- transfac.tfs.expressed

##
## Write expressed matrices for each prefix
##
expressedFamilyMats <- list()
ofile <- paste(annot.dir,"expressedFamilyMats.tsv",sep="/")
zz <- file(ofile, "w")  # open an output file connection
for ( emp in expressedMatsPrefix ){
  mems <- intersect(names(which(familyMap==emp)),expressedMats)
  if ( length(mems) > 0 ){
    expressedFamilyMats[[emp]] <- mems
    sumstring <- paste(emp,paste(mems,collapse=","),sep='\t')
    cat(sumstring,"\n",file=zz)
  }
}
close(zz)

ofile <- paste(annot.dir,"tteMaps.RData",sep="/")
save(transfac.tfs.expressed,pwm2tf.ps,pwm2tf.gname,pwm2tf.ncbiID,pwm2tf.tte.ps,pwm2tf.tte.gname,expressedMats,expressedMatsPrefix,tfps2pwm,allFams,fam2tf.tte.gname,fam2tf.gname,expressedFamilyMats,file=ofile)

##
##
## Decide whether to keep:
## which matrix families does a given tf map to?

## names(which(lapply((lapply(fam2tf.tte.gname,'==',tf1.cname)),sum)!=0))


