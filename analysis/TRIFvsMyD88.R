
filter <- lps.ps.sig
baddies <- names(which(cname.compare[filter]=="multiple"))
filter <- setdiff(filter,baddies)
 
measure <- all.mic[filter,"LPS"]-(all.mic[filter,"PAM3"]+all.mic[filter,"PolyIC"])
tt <- names(sort(measure,decreasing=TRUE)[1:50])

lps.mic <- all.mic[,"LPS"]

## pam3 polyIC signal integral
imax <- 4
m <- PAM3polyIC.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
pam3polyIC.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

measure.2 <- pam3polyIC.mic[filter]-lps.mic[filter]
tt.2 <- names(sort(measure.2,decreasing=TRUE)[1:50])

cname.compare[names(sort(measure.2,decreasing=TRUE)[1:30])]

## Sets

myd88.naive <- setdiff(lps.ps.sig,polyIC.ps.sig) ## "MyD88 depedent"
trif.naive <- setdiff(lps.ps.sig,pam3.ps.sig) ## "TRIF depedent"

## TRIF genes

## signal loss in the wild-type
trif.ko.effect <- trifKOlps.mic[filter] - all.mic[filter,"LPS"]
tt.3 <- names(sort(trif.ko.effect,decreasing=FALSE)[1:20])

# signaling difference LPS, but not PAM3 means it's TRIF
trif.signaling.diff <- all.mic[filter,"LPS"]-all.mic[filter,"PAM3"]

labs <- cname.compare[filter]
trif.signaling.diff <- trif.signaling.diff[setdiff(names(trif.signaling.diff),baddies)]
trif.ko.effect <- trif.ko.effect[setdiff(names(trif.ko.effect),baddies)]

labs[(trif.signaling.diff<1500) & (abs(trif.ko.effect)<1000)] <- ""
plot(trif.signaling.diff,trif.ko.effect,xlab="Difference between wild-type LPS signal and wild-type PAM3 signal",ylab="Difference between LPS induced signal in KO and wild-type", main="TRIF signaling",ylim=c(-5500,0),xlim=c(0,5500))
text(trif.signaling.diff,trif.ko.effect,labels=labs,pos=1)
abline(0,-1)

## generat this figure with axes
plot(trif.signaling.diff,all.mic[filter,"PolyIC"])


## MyD88 genes

## signal loss in the wild-type
myd88.ko.effect <- myd88KOlps.mic[filter] - all.mic[filter,"LPS"]

# signaling difference LPS, but not PolyIC means it's MyD88
myd88.signaling.diff <- all.mic[filter,"LPS"]- all.mic[filter,"PolyIC"]

x11()
labs <- cname.compare[filter]
labs[(myd88.signaling.diff<2000) & (abs(myd88.ko.effect)<3000)] <- ""
plot(myd88.signaling.diff,myd88.ko.effect,xlab="wt LPS, wt PolyIC signal difference",ylab="MyD88KO signal loss,LPS,wt" ,main="MyD88 signaling")
points(myd88.signaling.diff[myd88.signaling.diff>1000],myd88.ko.effect[myd88.signaling.diff>1000],col='blue')
text(myd88.signaling.diff,myd88.ko.effect,labels=labs,pos=3)
abline(0,-1)

## measures re-visited

myd88.ko.effect <- myd88KOlps.mic[filter] - all.mic[filter,"LPS"]


## generat this figure with axes

## could also consider TLR9=CpG?
plot(myd88.signaling.diff,all.mic[filter,"R848"])

### Compare gene sets

trif.keepers <- filter[(trif.signaling.diff>1500) & (trif.ko.effect < -1000)]


myd88.keepers <- filter[(myd88.signaling.diff>2000) & (myd88.ko.effect < -3000)]

both <- intersect(myd88.keepers,trif.keepers) 

trif.alone <- setdiff(trif.keepers,myd88.keepers)
myd88.alone <- setdiff(myd88.keepers, trif.keepers)



unaffected <- filter[(abs(myd88.ko.effect) < 1000) & (abs(trif.ko.effect) < 1000)]

x11()
labs <- cname.compare[filter]
labs[unaffected] <- ""
plot(myd88.ko.effect,trif.ko.effect,main="MyD88 and TRIF KO effects")
text(myd88.ko.effect,trif.ko.effect,labels=labs,pos=3)



total.mic <- cbind(all.mic[filter,],myd88KOlps.mic[filter],myd88KOpam3.mic[filter],myd88KOpolyIC.mic[filter],trifKOlps.mic[filter],trifKOpam2.mic[filter],myd88.signaling.diff[filter],trif.signaling.diff[filter], myd88.ko.effect[filter],trif.ko.effect[filter])

colnames(total.mic) <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG","MyD88 KO LPS","MyD88 KO PAM3", "MyD88 KO PolyIC","TRIF KO LPS","TRIF KO PAM2","MyD88 signal diff","TRIF signal diff","MyD88 KO Effect","TRIF KO Effect")


xgobi(total.mic,rowlab=cname.compare[filter])



## January 2006

## Find predominant clusters

gg <- names(sort(myd88.ko.effect)[1:100])

sortRows(goSignificance(np[gg]),1)[1:10,]

hh <- decimalRepRA.lps[gg]
sort(tapply(hh,as.factor(hh), length),decreasing=TRUE)


## March 2006
## Plots for 'Can we suggest where a mutation lies?' example

## Pretend the TRIF mutant is an unkown mutant

## Simplify some long names
cname.compare["1427932_s_at"]="1200003I10Rik"
cname.compare["1435137_s_at"]="1200015M12Rik"

unaffected <- filter[(abs(trif.ko.effect) < 2000) & (abs(myd88.signaling.diff) < 2000)]
labs <- cname.compare[filter]
labs[unaffected] <- ""

x11()

par(mfcol=c(1,2))

plot(trif.ko.effect,myd88.signaling.diff,main="\"Unknown\" mutant vs MyD88 signaling strength")
##text(trif.ko.effect[filter],myd88.signaling.diff[filter],labels=labs,pos=3)

plot(trif.ko.effect,trif.signaling.diff,main="\"Unknown\" mutant vs TRIF signaling strength")
##text(trif.ko.effect[filter],trif.signaling.diff[filter],labels=labs,pos=3)


## Pretend the MyD88 mutant is an unkown mutant

unaffected <- filter[(abs(myd88.ko.effect) < 2000) & (abs(myd88.signaling.diff) < 2000)]
labs <- cname.compare[filter]
labs[unaffected] <- ""

x11()

par(mfcol=c(1,2))

plot(myd88.ko.effect,trif.signaling.diff,main="\"Unknown\" effect mutant vs TRIF signaling strength",xlab="\"Unknown\" Mutant Effect",ylab="TLR TRIF signaling component")
##text(myd88.ko.effect[filter],trif.signaling.diff[filter],labels=labs,pos=3)

plot(myd88.ko.effect,myd88.signaling.diff,main="\"Unknown\" mutant effect vs MyD88 signaling strength",xlab="\"Unknown\" Mutant Effect",ylab="TLR MyD88 signaling component")
##text(myd88.ko.effect[filter],myd88.signaling.diff[filter],labels=labs,pos=3)



## Broadcast a matrix to the gaggle?


myd88.ko.effect <- myd88KOlps.mic - all.mic[lps.full.ps.sig,"LPS"]
trif.ko.effect <- trifKOlps.mic - all.mic[lps.full.ps.sig,"LPS"]

m1 <- cbind(myd88.ko.effect,trif.ko.effect)
rownames(m1) <- as.character(np[rownames(m1)])

broadcast (m1, "KO effect")

NFkB <- c("NP_033070","NP_033071","NP_033072","NP_032715","NP_062281")


### Irf 3 targets

source("/users/thorsson/macrophage/ArrayAnalysis/readTFTargets.R")

gg <- as.character(repProbes.np[Irf3.targets])

keepers <- c(gg[c(3,4,6,8)],gname2gid("G1p2"),gname2gid("Ccl5"))

## Works out to 
###  "1419282_at"   "1418930_at"   "1450783_at"   "1418293_at"   "1431591_s_at" "1418126_at"  


### October 2006

## Need to start from 'other' end and look at genes that do *not* have
## evidence for myd88 and/or trif signal

tt <- names(sort(lps.mus[names(which(abs(myd88.ko.effect)<25)),"min120"],decreasing=TRUE)[1:15])


filter1 <- (myd88.ko.effect<0)## KO indicates positive effect
##filter2 <- (myd88.ko.effect> (-200) ) ## KO is small 

filter3 <- lps.mic[lps.ps.sig] > 200  ## Look only at sizable genes
filter4 <- ( ( abs(myd88.ko.effect)/lps.mic[lps.ps.sig] ) < 0.20 )## KO effect is minor
##filter5 <- ( ( abs(myd88.ko.effect)/lps.mic[lps.ps.sig] ) > 0  )## KO effect is minor 
 
filterc <- filter1 & filter3 & filter4 

length(which(filterc))
ss  <- lps.ps.sig[filterc]

### October 2006

### TRIF and MyD88 wt signals revisited


x11()
labs <- cname.compare[filter]
unaffected <- filter[(abs(myd88.signaling.diff) < 1000) & (abs(trif.signaling.diff) < 1000)]
labs[unaffected] <- ""
plot(myd88.signaling.diff,trif.signaling.diff,main="MyD88 and TRIF Signaling Differences")
text(myd88.signaling.diff,trif.signaling.diff,labels=labs,pos=3)


## Doyle et al. Immunity Vol 17, p251, 2002

## TLR3/4-specific primary response genes
doyle.trif.oriented.cname <- c("Cxcl10","Ccl5","Ifnb1","Isg15","Ifit1")
## odd 1431591_s_at was previously Isg15(formerly G1p2) is now LOC677168 'similar to Isg15'
doyle.trif.oriented.cname <- c("Cxcl10","Ccl5","Ifnb1","LOC677168","Ifit1")
doyle.trif.oriented.ps <- as.character(repProbes.cname[doyle.trif.oriented.cname])

## This close to 1
trif.signaling.diff[doyle.trif.oriented.ps]/lps.mic[doyle.trif.oriented.ps]

## This is smallish, but more varialbe 
myd88.signaling.diff[doyle.trif.oriented.ps]/lps.mic[doyle.trif.oriented.ps]




x11()
labs <- cname.compare[filter]
unaffected <- filter[(abs(myd88.signaling.diff) < 1000) & (abs(trif.signaling.diff) < 1000)]
labs[unaffected] <- ""
plot(myd88.signaling.diff,trif.signaling.diff,main="MyD88 and TRIF Signaling Differences")
text(myd88.signaling.diff,trif.signaling.diff,labels=labs,pos=3)


trif.wt.signal.ratio <- trif.signaling.diff[filter]/lps.mic[filter]
myd88.wt.signal.ratio <- myd88.signaling.diff[filter]/lps.mic[filter]





filter1 <- ( abs(1-trif.wt.signal.ratio) < 0.05 )

filter1 <- ( abs(1-trif.wt.signal.ratio) < 0.1 )


filter2 <- ( myd88.wt.signal.ratio < 0.67 ) ## selected to cover doyle.trif.oriented !
filter3 <- lps.mic[filter] > 50

pure.trif <- filter[filter1 & filter2 & filter3 ]
