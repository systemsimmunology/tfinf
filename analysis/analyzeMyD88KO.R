################ Explore differences at time zero and due to KO alone

hh  <- names(sort(myd88KOlps.mus[,"min0"]/lps.mus[,"min0"],decreasing=TRUE))


jjn <- lps.mus[,c(1,2,4,6)]
jjd <- myd88KOlps.mus

jj <- lps.mus[,c(1,2,4,6)]/myd88KOlps.mus


kk <- names(sort(jj[,4]/jj[,1],decreasing=TRUE))


cbind(gid2gname(kk[1:50]),jjn[kk[1:50],])
cbind(cname.compare[kk[50:100]],jjn[kk[50:100],])

## Cumulative change score
sum(jjn["1419607_at",2:4]-jjn["1419607_at",1:3])
sum(jjd["1419607_at",2:4]-jjd["1419607_at",1:3])


nn1 <- apply((jjn[,2:4]-jjn[,1:3]),1,sum)
nn2 <- apply((jjd[,2:4]-jjd[,1:3]),1,sum)
nn <- names(sort(nn1/nn2,decreasing=TRUE))

#Wow this really looks like crud!

nn <- names(sort(nn1-nn2,decreasing=TRUE))

## Better


## Kawai 2001


kw  <- names(sort(myd88KOlps.mus[,"min120"]/myd88KOlps.mus[,"min0"],decreasing=TRUE)) # induced in the knockout. Rank by greatest induction.

kw.sig.frac <- myd88KOlps.mus[kw,"min120"]/lps.mus[kw,"min120"] ## less in the knockout than in wt at min 120

kwsfs <- names(sort(kw.sig.frac[1:50],decreasing=TRUE)) ## which are the "most reduced" in the knockout ( not used again ! )

kwrange <- 1:50

## Plot ko and wt at min120. Include only those genes ranking highly in "ko induction"
## We need to do this filtering - otherwise off-axis points may correspond to constitutively different between ko and wildtype

plot(log10(lps.mus[kw[kwrange],"min120"]),log10(myd88KOlps.mus[kw[kwrange],"min120"]))
text(log10(lps.mus[kw[kwrange],"min120"]),log10(myd88KOlps.mus[kw[kwrange],"min120"]),labels=cname.compare[kw[1:50]],pos=3)
abline(0,1)


kwrange <- 50:100

plot(log10(lps.mus[kw[kwrange],"min120"]),log10(myd88KOlps.mus[kw[kwrange],"min120"]))
text(log10(lps.mus[kw[kwrange],"min120"]),log10(myd88KOlps.mus[kw[kwrange],"min120"]),labels=cname.compare[kw[1:50]],pos=3)
abline(0,1)

#
# LPS PolyIC alone signal

lambda.cutoff <- 57.2

gg <- (lps.lambdas[,"min120"] > lambda.cutoff ) &
     (lps.lambdas[,"min120"] > lambda.cutoff ) &
     (pam2.lambdas[,"min120"] < lambda.cutoff/2 ) &
     (pam3.lambdas[,"min120"] < lambda.cutoff/2 ) &
     (polyIC.lambdas[,"min120"] < lambda.cutoff/2 ) &
     (r848.lambdas[,"min120"] < lambda.cutoff/2 ) &
     (cpg.lambdas[,"min120"] < lambda.cutoff/2 ) 

jj<- probeset[which(gg)]

## Magnitude increase. Works OK. But does turn up high intensity genes
## increasing slightly
mm <- apply(lps.mus[,2:6]-lps.mus[,1:5],1,sum)
mm.s <- names(sort(mm,decreasing=TRUE))

nn <- apply(lps.lambdas[,1:5],1,sum)
nn.s <- names(sort(nn,decreasing=TRUE))



### Ratio at one time point (this miss equivalences at an earlier timepoint)

x11()
ko.induction <- log10(myd88KOlps.mus[lps.full.ps.sig,"min120"]/myd88KOlps.mus[lps.full.ps.sig,"min0"])
wt.induction <- log10(lps.mus[lps.full.ps.sig,"min120"]/lps.mus[lps.full.ps.sig,"min0"])
labels <- cname.compare[lps.full.ps.sig]
plot(wt.induction, ko.induction, main="min120")
text(wt.induction, ko.induction, labels=labels, pos=3)


### Closest so far 
tt <- lps.full.ps.sig[which( abs(wt.induction-ko.induction) > log10(3)  )]
cname.compare[rownames(sortRows(lps.mus[tt,1:6],6,decreasing=TRUE))]

## A case where this fails miserably: Cxcl10="1418930_at"
## Very similar at 120 mins in wt and knockout
## The log difference of the difference was very sensitive to the t=0 difference between KO and wt


## Overall signal integral

imax <- 6
m <- lps.mus[lps.full.ps.sig,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
wt.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]


imax <- 4
m <- myd88KOlps.mus[lps.full.ps.sig,]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
myd88KOlps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

imax <- 4 ##  dropped the restriction to lps.ps.sig - not always relevant
m <- myd88KOpam3.mus
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
myd88KOpam3.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]

imax <- 4
m <- myd88KOpolyIC.mus[lps.full.ps.sig,]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,60,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
myd88KOpolyIC.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]


labels <- cname.compare[lps.full.ps.sig]
plot(wt.mic, myd88KOlps.mic, main="mean increase")
text(wt.mic, myd88KOlps.mic, labels=labels, pos=3)

biggies <- lps.full.ps.sig[which( (wt.mic>100)&(myd88KOlps.mic> (-200.) ))]
## Notes: some show increased profiles with knockout Ifnb1,Cd69






gg.m <- names(sort(myd88KOlps.mic[biggies]/wt.mic[biggies],decreasing=FALSE)[1:150])

gg.m2 <- names(sort(myd88KOlps.mus[gg.m,"min120"]/lps.mus[gg.m,"min120"],decreasing=FALSE))

## When considering MyD88 genes from WT, super high confidence as follows

myd88.ps.sig <- intersect(lps.ps.sig,intersect(pam2.ps.sig,intersect(pam3.ps.sig,intersect(r848.ps.sig,cpg.ps.sig))))

## down to 380, from 1303 in lps (of course some lps can be trif)

x11()
plot(all.mic[myd88.ps.sig,"LPS"], myd88KOlps.mic[myd88.ps.sig])
abline(0,1)

x11()
plot(all.mic[myd88.ps.sig,"PAM3"], myd88KOpam3.mic[myd88.ps.sig])
abline(0,1)

x11()
plot(all.mic[myd88.ps.sig,"PolyIC"], myd88KOpolyIC.mic[myd88.ps.sig])
abline(0,1)

### One more PolyIC induced but affected by MyD88
##> cname.compare["1453238_s_at"] 1453238_s_at "1200016E24Rik"
## Adapt calls this A130040M12Rik  a VL30 retroelement Interesting
## Search EnsMart. Ensembl says:
##1. Affymetrix Probe set: 1453238_s_at
##Affymetric probeset 1453238_s_at hits the genome in 1840 locations.

plot(all.mic[polyIC.ps.sig,"PolyIC"], myd88KOpolyIC.mic[polyIC.ps.sig])

# Doesn't seem to show much different than the myd88.ps.sig filter
## Result: Genes that are MyD88 dependent can also be expressed in PolyIC (presumably under TRIF dependence). However, the expression of those genes is independent of MyD88

### Normalization revisited

filter <- lps.ps.sig ## I think a large number of genes is good here

correction.factor <- quantile(all.mic[filter,"PAM2"],0.75)/quantile(all.mic[filter,"PAM3"],0.75)

plot(correction.factor*all.mic[filter,"PAM3"], all.mic[filter,"PAM2"])
abline(0,1)
abline(0,-1)

## filter dependent
## lps.ps.sig: 0.8096 union(pam2,pam3) 0.75 intersect(pam2,pam3) ).77

## Have to be careful here. These could be totally 'legit' differences, due to kinetics

