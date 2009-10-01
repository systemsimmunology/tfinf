
##
## insertNucLoc.R
## 
## Input scaled.mus.objects.RData
## nuclear localization profile for selected genes
## Output:scaled.mus.objects.RData taking into account nuclear localization
##
##
## Run after scaleIntensities.R, but only once !

t1 <- read.table('TFwesternScoredNuclear.tsv',as.is=TRUE,header=TRUE,row.names=1)
colnames(t1) <- c("min0","min30","min60","min120","min240","min360")
protein.nuclear.western <- as.matrix(t1)
rm(t1)
t.western <- c(0,30,60,120,240,360)

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
load(paste(annot.dir,"representativeProbes.RData",sep="/"))
load(paste(exp.dir,"all.mic.RData",sep="/"))
load(paste(exp.dir,"scaled.mus.objects.RData",sep="/"))
load(paste(exp.dir,"all.mus.objects.RData",sep="/"))

## 27 November
## Approximate nucl localizaton values, and import them into max.mat1
## Find changes relative to LPS in the mRNA, and apply that to the LPS nuc loc, to get the approximated profile

array.times.hr2 <- c(0,20,40,60,80,120)
array.times.hr12 <- c(array.times.hr2,240,480,720)

## Attempt to transform nucloc profiles for consistency with array profile

## 
timeShift <- function( profile, timefac ){
  sp <- approx(array.times.hr2,profile,array.times.hr2-timefac,rule=2)$y
  names(sp) <- paste("min",array.times.hr2,sep="")
  sp
}

nucloc.input.set <- c("Rel","Atf3","Egr1","Egr2","Fos","Rela")

stims <- colnames(all.mic)
stimslc <- c("lps","pam2","pam3","polyIC","r848","cpg")
colnames(all.mic) <- stimslc ## needed for proper indexing

## Rel
coi <- "Rel"
Rel.nuc.profile <- protein.nuclear.western[coi,]
psoi <- repProbes.cname[coi]
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
appr.protein.nuclear.western  <- numeric()
for ( i in 2:6 ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[i] *  Rel.nuc.profile )
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Rel.appr.protein.nuclear.western <- appr.protein.nuclear.western 


## Atf3
Atf3.nuc.profile <- protein.nuclear.western["Atf3",]
psoi <- repProbes.cname["Atf3"]
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
scalefacs["polyIC"] <- 1.
scalefacs["cpg"] <- 1.
## The later two mostly appear shifted w.r.t. lps
tb <- numeric(length=5)
names(tb) <- stimslc[2:6]
tb["pam2"] <- 0.
tb["pam3"] <- 0.
tb["polyIC"] <- 60.
tb["r848"] <- 0.
tb["cpg"] <- 60.
appr.protein.nuclear.western  <- numeric()
for ( stim in stimslc[2:6] ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[stim]*timeShift(Atf3.nuc.profile,tb[stim] ))
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Atf3.appr.protein.nuclear.western <- appr.protein.nuclear.western 


## Egr1
Egr1.nuc.profile <- protein.nuclear.western["Egr1",]
psoi <- repProbes.cname["Egr1"]
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
tb <- numeric(length=5)
names(tb) <- stimslc[2:6]
tb["pam2"] <- 0.
tb["pam3"] <- 0.
tb["polyIC"] <- 60.
tb["r848"] <- 0.
tb["cpg"] <- 60.
appr.protein.nuclear.western  <- numeric()
for ( stim in stimslc[2:6] ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[stim]*timeShift(Egr1.nuc.profile,tb[stim] ))
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Egr1.appr.protein.nuclear.western <- appr.protein.nuclear.western 


## Egr2
Egr2.nuc.profile <- protein.nuclear.western["Egr2",]
psoi <- repProbes.cname["Egr2"]
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
tb <- numeric(length=5)
names(tb) <- stimslc[2:6]
tb["pam2"] <- 0.
tb["pam3"] <- 0.
tb["polyIC"] <- 60.
tb["r848"] <- 0.
tb["cpg"] <- 60.
appr.protein.nuclear.western  <- numeric()
for ( stim in stimslc[2:6] ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[stim]*timeShift(Egr2.nuc.profile,tb[stim] ))
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Egr2.appr.protein.nuclear.western <- appr.protein.nuclear.western 

## Fos
Fos.nuc.profile <- protein.nuclear.western["Fos",]
psoi <- repProbes.cname["Fos"]
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
tb <- numeric(length=5)
names(tb) <- stimslc[2:6]
tb["pam2"] <- 0.
tb["pam3"] <- 0.
tb["polyIC"] <- 60.
tb["r848"] <- 0.
tb["cpg"] <- 60.
appr.protein.nuclear.western  <- numeric()
for ( stim in stimslc[2:6] ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[stim]*timeShift(Fos.nuc.profile,tb[stim] ))
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Fos.appr.protein.nuclear.western <- appr.protein.nuclear.western 

## Rela
### Don't know if mic levels are the right scale factor
### e.g. PolyIC goes really low after normalizatioin, tho intensity is on
## the high side. Intensity matters for beta multiplier.
coi <- "Rela"
Rela.nuc.profile <- protein.nuclear.western[coi,]
psoi <- repProbes.cname[coi]
## mic normalization
scalefacs <- all.mic[psoi,]/all.mic[psoi,"lps"]
scalefacs["polyIC"] <- abs(scalefacs["polyIC"] )
## abs val normalization
scalefacs <- c(lps.mus[psoi,"min60"],pam2.mus[psoi,"min60"],pam3.mus[psoi,"min60"],polyIC.mus[psoi,"min60"],r848.mus[psoi,"min60"],cpg.mus[psoi,"min60"])/lps.mus[psoi,"min60"]
## polyIC scalefac is negative. PolyIC mRNA hardly changes over the first two hours.
appr.protein.nuclear.western  <- numeric()
for ( i in 2:6 ){
  appr.protein.nuclear.western <- rbind(appr.protein.nuclear.western, scalefacs[i] *  Rela.nuc.profile )
}
rownames(appr.protein.nuclear.western) <- stimslc[2:6]
Rela.appr.protein.nuclear.western <- appr.protein.nuclear.western  


#### Incorporate into mat max
stim <- "pam2"
for ( coi in nucloc.input.set ){
  psoi <- repProbes.cname[coi]
  wprof <- eval(parse(text=paste(coi,".appr.protein.nuclear.western",sep="")))
  rhs <-  approx(t.western,wprof[stim,],array.times.hr2,rule=2)$y
  pam2.mat.max1[psoi,] <- approx(t.western,wprof[stim,],array.times.hr2,rule=2)$y
}
stim <- "pam3"
for ( coi in nucloc.input.set ){
  psoi <- repProbes.cname[coi]
  wprof <- eval(parse(text=paste(coi,".appr.protein.nuclear.western",sep="")))
  rhs <-  approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
  pam3.mat.max1[psoi,] <- approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
}
stim <- "polyIC"
for ( coi in nucloc.input.set ){
  psoi <- repProbes.cname[coi]
  wprof <- eval(parse(text=paste(coi,".appr.protein.nuclear.western",sep="")))
  rhs <-  approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
  polyIC.mat.max1[psoi,] <- approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
}
stim <- "r848"
for ( coi in nucloc.input.set ){
  psoi <- repProbes.cname[coi]
  wprof <- eval(parse(text=paste(coi,".appr.protein.nuclear.western",sep="")))
  rhs <-  approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
  r848.mat.max1[psoi,] <- approx(t.western,wprof[stim,],array.times.hr12,rule=2)$y
}
stim <- "cpg"
for ( coi in nucloc.input.set ){
  psoi <- repProbes.cname[coi]
  wprof <- eval(parse(text=paste(coi,".appr.protein.nuclear.western",sep="")))
  rhs <-  approx(t.western,wprof[stim,],array.times.hr2,rule=2)$y
  cpg.mat.max1[psoi,] <- approx(t.western,wprof[stim,],array.times.hr2,rule=2)$y
}

ofile <- paste(exp.dir,"scaled.mus.objects.RData",sep="/")
savelist <- c("max.intensity","lps.mat.max1","lps.mat.max1.exp","pam2.mat.max1","pam3.mat.max1","polyIC.mat.max1","r848.mat.max1","cpg.mat.max1" )
save(list=savelist, file=ofile)
