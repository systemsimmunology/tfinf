## Overall signal integral
imax <- 6
m <- lps.mus[lps.full.ps.sig,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
wt.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

m <- trifKOlps.mus[lps.full.ps.sig,]
imax <- 3
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
trifKOlps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

m <- trifKOpam2.mus[pam2.ps.sig,]
imax <- 3
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
trifKOpam2.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

##
## Signal integral combinations
##

trif.ko.effect <- all.mic[filter,"LPS"] - trifKOlps.mic[filter]
trif.signaling.diff <- all.mic[filter,"LPS"]-all.mic[filter,"PAM3"]
polyIC.mic <- all.mic[,"PolyIC"]

polyIC.part <- lps.mic[filter]-polyIC.mic[filter] ## contribution of Poly IC to LPS signal
polyIC.frac <- polyIC.mic[filter]/lps.mic[filter] ## fraction of LPS signal accounted for by PolyIC
trif.ko.frac <- trif.ko.effect/lps.mic[filter] ## trif ko effect as fraction of LPS
trif.signaling.diff.frac <- trif.signaling.diff/lps.mic[filter] ## LPS, PAM3 signal difference as fraction


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

elps  <- names(which(lps.lambdas[,"min120"] > lambda.cutoff ))
inducts <- names(which(lps.ratios[,"min120"] > 0  ))
ei <- intersect(elps,inducts)

x <- trifko.frac[ei]
y <- polyIC.frac[ei]
z <- lpspam3diff.frac[ei]             

vset <- names(which(TRIFness[ei]> -0.5 & TRIFness[ei]<0.5))
