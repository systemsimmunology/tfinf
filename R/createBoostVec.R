##target.pull <- 60

array.times <- c(0,20,40,60,80,120,240,360,480,1080,1440)

t.index.max <- 11

input.boost.time <- 60

inputs.boosted <- matrix(nrow=length(transfac.tfs.expressed),ncol=t.index.max)
rownames(inputs.boosted) <- transfac.tfs.expressed
for ( psoii in transfac.tfs.expressed ){
  inputs.boosted[psoii,] <- approx(array.times,lps.mat.max1[psoii,],array.times-input.boost.time,rule=2)$y
}

boost.vec <- rep(input.boost.time,length(transfac.tfs.expressed))
names(boost.vec) <- transfac.tfs.expressed


early.tfs <- names(halfChangeAt[halfChangeAt %in% c("min20","min40","min60","min80")])
mid.tfs <- names(halfChangeAt[halfChangeAt %in% c("min120","hr2","hr4")])
late.tfs <- names(halfChangeAt[halfChangeAt %in% c("hr6","hr8","hr18","hr24")])
boost.vec[early.tfs] <- 0
boost.vec[mid.tfs] <- 30
boost.vec[late.tfs] <- 60


##targ.exp.pulled <- approx(array.times,lps.mat.max1[targ,],array.times+target.pull,rule=2)$y

## For some input variables, use the nuclear protein
## lps.mat.max1 will be used for the new input variables
## lps.mat.max1 will remain unboosted
## 
## For now, use no 'boost' time

nucloc.input.set <- c("Rel","Atf3","Egr1","Egr2","Fos","Rela")
for ( coi in nucloc.input.set ){
  boost.vec[repProbes.cname[coi]] <- 0
  lps.mat.max1[repProbes.cname[coi],] <- approx(t.western,protein.nuclear.western[coi,],t.array,rule=2)$y
}


if ( !is.null(boost.vec.refit) ){
  boost.vec[names(boost.vec.refit)] <- boost.vec.refit
}

## June 2007

## For randomization trials,  create a boost.vec for other TFS

addeds <- setdiff(jared.tfs.expressed,transfac.tfs.expressed )
boost.vec.addeds <- rep(60,length=length(addeds))
names(boost.vec.addeds) <- addeds
boost.vec <- c(boost.vec,boost.vec.addeds)


##boost.vec.mat <- sortRows(cbind(cc[names(boost.vec)],ncbiID[names(boost.vec)],names(boost.vec),boost.vec),1)
##colnames(boost.vec.mat) <- c("MGNC","EntrezID","Probeset","Delay")
##write.table(boost.vec.mat,file="TFDelays.tsv",quote=FALSE,row.names=FALSE,sep='\t')
