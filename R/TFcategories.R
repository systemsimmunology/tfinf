
conds <-  c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")

## Rough divisions for prioritization

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")

load(paste(annot.dir,"annotation.objects.RData",sep="/"))
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(annot.dir,"all.ps.list.objects.RData",sep="/"))
load(paste(annot.dir,"tteMaps.RData",sep="/"))

tte <- transfac.tfs.expressed
cc <- cname.compare

load(paste(exp.dir,"all.mus.objects.RData",sep="/"))

## Maximum expression level

maxLev <- apply(lps.mus[transfac.tfs.expressed,],1,max)

mag <- vector(length=length(transfac.tfs.expressed))
names(mag) <- transfac.tfs.expressed

mag[names(which(maxLev>900))] <- "High"
mag[names(which((maxLev>300)&(maxLev<900)))] <- "Middle"
mag[names(which(maxLev<300))] <- "Low"

## Up or down overall ?

ud <- vector(length=length(transfac.tfs.expressed))
names(ud) <- transfac.tfs.expressed

imax <- 11
m <- lps.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120,60*c(4,6,8,18,24))
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
lps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]


ud[names(which(lps.mic[transfac.tfs.expressed]>=0))] <- "Increasing"
ud[names(which(lps.mic[transfac.tfs.expressed]<0))] <- "Decreasing"


## Repair by hand (mic score may look unnatural  up when ups and downs compete )
## Perhaps <= hr4-6 should carry most weight

ud[repProbes.cname["Ahr"]] <- "Increasing"
ud[repProbes.cname["Atf2"]] <- "Increasing"
ud[repProbes.cname["Elk3"]] <- "Decreasing"
ud[repProbes.cname["Mafg"]] <- "Decreasing"
ud[repProbes.cname["Myc"]] <- "Increasing"
ud[repProbes.cname["Nfe2l2"]] <- "Increasing"
ud[repProbes.cname["Nfyb"]] <- "Increasing"


## Where is the extreme value 

ups <- names(which(ud=="Increasing"))
downs <- names(which(ud=="Decreasing"))

extremeAt <- vector(length=length(transfac.tfs.expressed))
names(extremeAt) <- transfac.tfs.expressed

extremeAt[ups] <- colnames(lps.mus)[apply(lps.mus[ups,],1,which.max)]
extremeAt[downs] <- colnames(lps.mus)[apply(lps.mus[downs,],1,which.min)]

## Estimate where signal reaches has maximum for ups

## Compute the time point at which half the maximum is exceeded
halfmax <- (apply(lps.mus,1,max)+apply(lps.mus,1,min))/2
mw <- function( vec ){ min(which(vec)) }
jj <- (lps.mus[ups,]-halfmax[ups])>0
tup <- apply(jj,1,mw)

halfmin <- halfmax ## odd, but it works out to the same expression
jj <- (lps.mus[downs,]-halfmin[downs])<0
tdown <- apply(jj,1,mw)

halfChangeAt <- conds[c(tup,tdown)[tte]] ## concanenate and rearrange 
names(halfChangeAt) <- transfac.tfs.expressed ## (was lost by the conds operation )


## Repair by hand (max min construct really only works for unimodal )
## These will do better with an 'initial baseline' estimate

halfChangeAt[repProbes.cname["Atf2"]] <- "hr4"
halfChangeAt[repProbes.cname["Bhlhb2"]] <- "hr2"
halfChangeAt[repProbes.cname["Crem"]] <- "hr2"
halfChangeAt[repProbes.cname["E2f1"]] <- "hr4"
halfChangeAt[repProbes.cname["Elk3"]] <- "hr4"
halfChangeAt[repProbes.cname["Mafg"]] <- "hr4"
halfChangeAt[repProbes.cname["Myc"]] <- "min60"
halfChangeAt[repProbes.cname["Nfe2l2"]] <- "min40"
halfChangeAt[repProbes.cname["Nfyb"]] <- "hr4"


### Output


write.table(file="ttfEcluster.tsv",
            cbind(transfac.tfs.expressed,mag,ud,halfChangeAt,extremeAt),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

### Mappings

mapvec <- character(length=length(transfac.tfs.expressed))
for ( i in 1:length(transfac.tfs.expressed)){
  ##mapvec[i] <- paste(tfps2pwm[[transfac.tfs.expressed[i]]],collapse=",")
  fams <- unique(sort(sapply(tfps2pwm[[transfac.tfs.expressed[i]]],split1)))
  mapvec[i] <- paste(fams,collapse=",")
}


write.table(file="ttfEpwms.tsv",
            cbind(transfac.tfs.expressed,mapvec),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')



### Plots 
tte.spreadsheet.order <- as.character(read.table("tfsSpreadSheetOrder.tsv",as.is=TRUE)$V1)


### Strong, induced TFs, and matrices

bigUps <- names(which((ud=="Increasing") & (mag=="High")))
midUps <- names(which((ud=="Increasing") & (mag=="Middle")))


## arrange by time order of induction
ind.cond.numeric <- match(halfChangeAt[bigUps],conds)
bu.cname <- names(cc[bigUps])

bigUpFams <- us(sapply(unlist(tfps2pwm[bigUps]),split1))

om <- cbind(ind.cond.numeric,bu.cname)
rownames(om) <- bigUps

#om <- data.frame(om)

### Strong, induced TFs, and matrices

bigDowns <- names(which((ud=="Decreasing") & (mag=="High")))
midDowns <- names(which((ud=="Decreasing") & (mag=="Middle")))

## arrange by time order of induction
ind.cond.numeric <- match(halfChangeAt[bigDowns],conds)
bd.cname <- names(cc[bigDowns])

bigDownFams <- us(sapply(unlist(tfps2pwm[bigDowns]),split1))

bd.om <- cbind(ind.cond.numeric,bd.cname)
rownames(bd.om) <- bigDowns

## Output

ofile <- paste(annot.dir,"TFcategories.RData",sep="/")
save(mag,ud,halfChangeAt,extremeAt,bigUps,bigDowns,midUps,midDowns,file=ofile)

