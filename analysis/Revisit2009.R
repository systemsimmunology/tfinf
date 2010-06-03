
expression.directory <- file.path(Sys.getenv("AA"),"data/20100426.curated.3prime")
annotation.directory <- file.path(Sys.getenv("AA"),"annotation")

load(file.path(expression.directory,"all.mic.RData"))
load(file.path(expression.directory,"all.mus.objects.RData"))
load(file.path(expression.directory,"all.cube.objects.RData"))

load(file.path(annotation.directory,"all.ps.list.objects.RData"))

##mtdir <- "/Users/thorsson/macrophage/MyD88vsTRIF/R"
##tifdir <- "/Users/thorsson/macrophage/TFinfluence/R"
##setwd(tifdir)
##source("/Users/thorsson/macrophage/TFinfluence/R/initialize.R")
##source("/Users/thorsson/macrophage/AffyArray/newAffy/utilities.R")

filter <- lps.full.ps.sig
labs <- cname.compare[filter]
labs[(trif.signaling.diff<1500) & (abs(trif.ko.effect)<1000)] <- ""

 
plot(trif.signaling.diff,trif.ko.effect,xlab="wt LPS, wt PAM3 signal difference",ylab="TRIF KO signal loss,LPS,wt", main="TRIF signaling")
text(trif.signaling.diff,trif.ko.effect,labels=labs,pos=3)
abline(0,-1)




### July 2009. Considering again to include PolyIC signal strength
polyIC.mic <- all.mic[,"PolyIC"]
trif.ko.effect <- - trif.ko.effect ## redefine. All these negatives are hurting my brain.
polyIC.part <- lps.mic[filter]-polyIC.mic[filter] ## contribution of Poly IC to LPS signal

 
polyIC.frac <- polyIC.mic[filter]/lps.mic[filter] ## fraction of LPS signal accounted for by PolyIC
trif.ko.frac <- trif.ko.effect/lps.mic[filter] ## trif ko effect as fraction of LPS
trif.signaling.diff.frac <- trif.signaling.diff/lps.mic[filter] ## LPS, PAM3 signal difference as fraction

x11() 

labs <- cname.compare[filter]

cond1 <- trif.ko.frac < 0.6
cond2 <- polyIC.frac  < 0.3
cond3 <- polyIC.mic[filter] < 30

labs[(cond1|cond2|cond3)] <- ""
 
plot(trif.ko.frac,polyIC.frac[filter],xlab="frac. TRIF signal by KO",ylab="PolyIC fraction of LPS",
     xlim=c(0,1.2),ylim=c(0,1.2),main="TRIF signaling")
text(trif.ko.frac,polyIC.frac[filter],labels=labs,pos=3)
abline(1,0)
lines(c(1,1),c(-10^9,10^9))

trif.ko.frac[psoi];!cond1[psoi];polyIC.frac[psoi];!cond2[psoi];polyIC.mic[eid];!cond3[psoi];labs[psoi]

plot(trif.ko.effect,polyIC.part,xlab="TRIF KO signal loss,LPS,wt",ylab="PolyIC contribu. of LPS", main="TRIF signaling")

## April 2010

## TRIFness
polyIC.frac <- polyIC.frac
trifko.frac <- trif.ko.frac
lpspam3diff.frac <- trif.signaling.diff.frac

w1 <- 1
w2 <- 1
w3 <- 1
## Dot product with 1,1,1 (or weight vector)
dotproduct <- w1*polyIC.frac + w2*trifko.frac + w3*lpspam3diff.frac
lengthw <- sqrt(w1^2 + w2^2 + w3^2)
lengthvec <- sqrt(polyIC.frac^2 + trifko.frac^2 + lpspam3diff.frac^2)
angle <- acos(dotproduct/lengthw/lengthvec)

v <- names(sort(angle,decreasing=FALSE)[1:40])
vv <- names(which(lps.mic[v] > 600 ))
psoi <- rc["Ccl4"]
c(polyIC.frac[psoi],trifko.frac[psoi],lpspam3diff.frac[psoi])

## Distance to the "ideal" point
eucdist.trif <- sqrt((polyIC.frac-w1)^2 + (trifko.frac-w2)^2 + (lpspam3diff.frac-w3)^2)

v <- names(sort(eucdist.trif,decreasing=FALSE)[1:40])
vv <- names(which(lps.mic[v] > 200 ))

psoi <- rc["Ccl4"]
c(polyIC.frac[psoi],trifko.frac[psoi],lpspam3diff.frac[psoi])

##TRIFness <- sqrt(w1*polyIC.frac^2+w2*trifko.frac+w3*lpspam3diff.frac^2)/abs(w1+w2+w3)

vset <- names(which(TRIFness<1.1 & TRIFness>0.9))

## MyD88 ness
pam3.mic <- all.mic[,"PAM3"]
myd88.ko.effect <- - myd88.ko.effect ## redefine. All these negatives are hurting my brain.
 
pam3.frac <- pam3.mic[filter]/lps.mic[filter] ## fraction of LPS signal accounted for by PAM3
myd88.ko.frac <- myd88.ko.effect/lps.mic[filter] ## trif ko effect as fraction of LPS
myd88.signaling.diff.frac <- myd88.signaling.diff/lps.mic[filter] ## LPS, PolyIC signal difference as fraction

## MYD88ness
pam3.frac <- pam3.frac
myd88ko.frac <- myd88.ko.frac
lpspolyICdiff.frac <- myd88.signaling.diff.frac

w1 <- 1
w2 <- 1
w3 <- 1

## Distance to the "ideal" point
eucdist.myd88 <- sqrt((pam3.frac-w1)^2 + (myd88ko.frac-w2)^2 + (lpspolyICdiff.frac-w3)^2)

v <- names(sort(eucdist.myd88,decreasing=FALSE)[1:40])
vv <- names(which(lps.mic[v] > 300 ))

psoi <- rc["Ccl4"]
c(pam3.frac[psoi],myd88ko.frac[psoi],lpspolyICdiff.frac[psoi])

## Dot product with 1,1,1 (or weight vector)
dotproduct <- w1*pam3.frac+w2*myd88ko.frac+w3*lpspolyICdiff.frac
lengthw <- sqrt(w1^2+w2^2+w3^2)
lengthvec <- sqrt(pam3.frac^2+myd88ko.frac^2+lpspolyICdiff.frac^2)
angle <- acos(dotproduct/lengthw/lengthvec)

## Not sure how useful it is to plot one against the other
plot(eucdist.trif[higher.mic],eucdist.myd88[higher.mic],xlim=c(0,2),ylim=c(0,2))
text(eucdist.trif[higher.mic],eucdist.myd88[higher.mic],label=cc[higher.mic])

cc[names(which(MyD88ness<1.05 & MyD88ness>0.95))]

higher.mic <- names(which(lps.mic > 300))

