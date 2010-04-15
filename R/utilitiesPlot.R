
##  plotter for WT intensities
paiATF3 <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr6","hr8","hr18","hr24")

  ymax <- max(c(lps.mus[probeID,], atf3KOlps.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)

  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  lines( cpg.mus[probeID,], col='yellow',lwd=3)

  lines( c(1,4,6,7), atf3KOlps.mus[probeID,], col='blue',lwd=3,lty=2)
  lines( c(1,4,6), atf3KOpam2.mus[probeID,], col='green',lwd=3,lty=2)
  lines( c(1,4,6), atf3KOpolyIC.mus[probeID,], col='magenta',lwd=3,lty=2)
  lines( c(1,4,6), atf3KOcpg.mus[probeID,], col='yellow',lwd=3,lty=2)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PolyIC","CpG","ATF3KO LPS","ATF3KO PAM2","ATF3KO PolyIC","ATF3KO CpG"),
       col=c('blue','green','magenta','yellow','blue','green','magenta','yellow'),
       lty=c(1,1,1,1,2,2,2,2),
       lwd=3)
  
}

## A plotter for WT intensities
paiTLR7v9 <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  ymax <- max(c(r848.mus[probeID,], cpg.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( r848.mus[probeID,], ylim=c(ymin,ymax), col='purple', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines( cpg.mus[probeID,], col='yellow',lwd=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("R848","CpG"),
       col=c('purple','yellow'),
       lty=c(1,1),
       lwd=3)
  
}

## A plotter for WT intensities
paiTLR2 <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  ymax <- max(c(pam2.mus[probeID,], pam3.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( pam2.mus[probeID,], ylim=c(ymin,ymax), col='green', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines( pam3.mus[probeID,], col='red',lwd=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("PAM2","PAM3"),
       col=c('green','red'),
       lty=c(1,1),
       lwd=3)
  
}


paiWT <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr6","hr8","hr12","hr18","hr24")

  ymax <- max(c(pam2.mus[probeID,],lps.mus[probeID,], pam3.mus[probeID,],polyIC.mus[probeID,], r848.mus[probeID,], cpg.mus[probeID,],polyIC.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

#  plot(c(1:(length(conds)+1)),lps.mus[probeID,], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3, plot=FALSE)

  plot(c(1:(length(conds)+1)),c(1:(length(conds)+1)) , ylim=c(ymin,ymax), col='blue', type='n',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)

  lines( c(1:9,11:12), lps.mus[probeID,], col='blue', lwd=3 )
  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines(  c(1:7,9:10), pam3.mus[probeID,], col='red',lwd=3)
  lines( c(1:7,9:10), polyIC.mus[probeID,], col='magenta',lwd=3)
  lines( c(1:7,9:10), r848.mus[probeID,], col='purple',lwd=3)
  lines( cpg.mus[probeID,], col='yellow',lwd=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG"),
       col=c('blue','green','red','magenta','purple','yellow'),
       lty=c(1,1,1,1,1,1),
       lwd=3)
  
  
}

## A plotter for WT intensities
paiWT.old <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  ymax <- max(c(pam2.mus[probeID,],lps.mus[probeID,1:9], pam3.mus[probeID,],polyIC.mus[probeID,], r848.mus[probeID,], cpg.mus[probeID,],polyIC.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,1:9], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( pam3.mus[probeID,], col='red',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  lines( r848.mus[probeID,], col='purple',lwd=3)
  lines( cpg.mus[probeID,], col='yellow',lwd=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG"),
       col=c('blue','green','red','magenta','purple','yellow'),
       lty=c(1,1,1,1,1,1),
       lwd=3)
  
  
}

## A plotter for Primed  intensities
paiPriming <- function(probeID){

  
  conds <- c("min20","min40","min60","min80","hr2","hr4","hr6","hr8","hr18","hr24")
  all.conds <- c("min0","min20","min40","min60","min80","hr2","hr4","hr6","hr8","hr18","hr24")
  lps.conds <- c("min0","min20","min40","min60","min80","hr2","hr4","hr6","hr8","hr18","hr24")
  priming.conds <- c("min0","min20","min40","min60","min80","hr2","hr4")

  ymax <- max(lps.mus[probeID,1:11],zymosan.lps.mus[probeID,],pma.lps.mus[probeID,],zymosan.mus[probeID,],pma.mus[probeID,])
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,1:11], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)


  lines(c(1,2,3,4,5,6,7),zymosan.lps.mus[probeID,], col='red',lwd=3,lty=2)
  lines(c(1,2,3,4,5,6,7),pma.lps.mus[probeID,], col='red',lwd=3,lty=3)


  lines(c(5,6,7),zymosan.mus[probeID,], col='green',lwd=3,lty=2)
  lines(c(5,6,7),pma.mus[probeID,], col='green',lwd=3,lty=3)

#  points(6,zymosan.mus[probeID],pch=19,col='blue')
#  points(6,pma.mus[probeID],pch=17,col='blue')
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)

  legend (1,0.95*ymax,legend=c("LPS","LPS+Zymosan","LPS+PMA","Zymosan","PMA"),
         col=c('blue','red','red','green','green'),
         lty=c(1,2,3,2,3),
         lwd=3)
          
  ## This was for mixed lines and points!
  ##  legend(1,0.95*ymax,legend=c("LPS","LPS+Zymosan","LPS+PMA","Zymosan","PMA"),
  ##       col=c('blue','blue','blue','blue','blue'),
  ##       lty=c(1,2,3,NA,NA),
  ##       pch=c(NA,NA,NA,19,17),
  ##       lwd=3)
  
  
}



## A plotter illustrating Myd88 dependence
paiMyD88 <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")


  wtrange <- 1:6
  ymax <- max(c(lps.mus[probeID,wtrange],myd88KOlps.mus[probeID,],trifKOlps.mus[probeID,],pam2.mus[probeID,],pam3.mus[probeID,],polyIC.mus[probeID,]),myd88KOpam3.mus[probeID,])
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,wtrange], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines(c(1,2,4,6),myd88KOlps.mus[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,2,4,6),myd88KOpam3.mus[probeID,], col='red',lwd=3,lty=2)
  lines(c(1,4,6),trifKOlps.mus[probeID,], col='blue',lwd=3,lty=3)

  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( pam3.mus[probeID,], col='red',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS",
              "MyD88 KO LPS","TRIF KO LPS","PAM2","PAM3","MyD88 KO PAM3","PolyIC"),
       col=c('blue','blue','blue','green','red','red','magenta'),
       lty=c(1,2,3,1,1,2,1),
       lwd=3)
  
  
}



## A plotter illustrating trif dependence
paiTRIF <- function(probeID,ymax=NULL){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  wtrange <- 1:6
  if ( is.null(ymax) ){
    ymax <- max(c(lps.mus[probeID,wtrange],myd88KOlps.mus[probeID,],trifKOlps.mus[probeID,],pam2.mus[probeID,],pam3.mus[probeID,],polyIC.mus[probeID,]))
  }
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,wtrange], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines(c(1,2,4,6),myd88KOlps.mus[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,4,6),trifKOlps.mus[probeID,], col='blue',lwd=3,lty=3)

  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( pam3.mus[probeID,], col='red',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS",
              "MyD88 KO LPS","TRIF KO LPS","PAM2","PAM3","PolyIC"),
       col=c('blue','blue','blue','green','red','magenta'),
       lty=c(1,2,3,1,1,1),
       lwd=3)
  
  
}

paiLPS <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")


  wtrange <- 1:6
  ymax <- max(c(lps.mus[probeID,wtrange],myd88KOlps.mus[probeID,],trifKOlps.mus[probeID,]))
  ymin <- 0

  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  x11()

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)

  par(op)

  plot( lps.mus[probeID,wtrange], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,
main=main,lwd=3)

  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)

  lines(c(1,2,4,6),myd88KOlps.mus[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,4,6),trifKOlps.mus[probeID,], col='blue',lwd=3,lty=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)

  par(op)
  legend(1,0.95*ymax,legend=c("LPS",
              "MyD88 KO LPS","TRIF KO LPS"),
       col=c('blue','blue','blue'),
       lty=c(1,2,3),
       lwd=3)


}
        

plotNormAffyIntensities <- function(probeID){

    conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")


  ymax <- max(c(pam2.mus.renorm[probeID,],lps.mus.renorm[probeID,2:9], pam3.mus.renorm[probeID,],polyIC.mus.renorm[probeID,], r848.mus.renorm[probeID,], cpg.mus.renorm[probeID,],polyIC.mus.renorm[probeID,]))
  
  ymin <- 0
  
  main <- paste(gid2gname(probeID),", Probe ID:",probeID,sep="")

  x11()
  
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus.renorm[probeID,2:9], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  

  axis(1, 1:length(conds),labels=conds)
  axis(2)
  
  lines( pam2.mus.renorm[probeID,], col='green',lwd=3)
  lines( pam3.mus.renorm[probeID,], col='red',lwd=3)
  lines( polyIC.mus.renorm[probeID,], col='magenta',lwd=3)
  lines( r848.mus.renorm[probeID,], col='purple',lwd=3)
  lines( cpg.mus.renorm[probeID,], col='yellow',lwd=3)

  lines(c(1,3,5),myd88KOlps.mus.renorm[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,3,5),myd88KOpam3.mus.renorm[probeID,], col='red',lwd=3,lty=2)
  lines(c(1,3,5),myd88KOpolyIC.mus.renorm[probeID,], col='magenta',lwd=3,lty=2)
  lines(c(1,4,6),trifKOlps.mus[probeID,], col='blue',lwd=3,lty=3)
  
  lines(c(1,3,5),pam3polyIC.mus.renorm[probeID,], col='skyblue',lwd=3)
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG",
              "MyD88KO LPS","MyD88KO PAM3","MyD88KO PolyIC","TRIF KO LPS","PAM3+PolyIC"),
       col=c('blue','green','red','magenta','purple','yellow','blue','red','magenta','blue','skyblue'),
       lty=c(1,1,1,1,1,1,2,2,2,3,1),
       lwd=3)
  
  
}


## Plot ratios for all condtions for a single probe
pai <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  ymax <- max(c(pam2.mus[probeID,],lps.mus[probeID,1:9], pam3.mus[probeID,],polyIC.mus[probeID,], r848.mus[probeID,], cpg.mus[probeID,],polyIC.mus[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  x11()
  
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,1:9], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( pam3.mus[probeID,], col='red',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  lines( r848.mus[probeID,], col='purple',lwd=3)
  lines( cpg.mus[probeID,], col='yellow',lwd=3)

  lines(c(1,2,4,6),myd88KOlps.mus[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,2,4,6),myd88KOpam3.mus[probeID,], col='red',lwd=3,lty=2)
  lines(c(1,2,4,6),myd88KOpolyIC.mus[probeID,], col='magenta',lwd=3,lty=2)
  lines(c(1,4,6),trifKOlps.mus[probeID,], col='blue',lwd=3,lty=3)
  lines(c(1,4,6),trifKOpam2.mus[probeID,], col='green',lwd=3,lty=3)
  
  lines(c(1,2,4,6),pam3polyIC.mus[probeID,], col='skyblue',lwd=3)
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG",
              "MyD88KO LPS","MyD88KO PAM3","MyD88KO PolyIC","TRIFKO LPS","TRIFKO PAM2","PAM3+PolyIC"),
       col=c('blue','green','red','magenta','purple','yellow','blue','red','magenta','blue','green','skyblue'),
       lty=c(1,1,1,1,1,1,2,2,2,3,3,1),
       lwd=3)
  
  
}


## plot nuclear localizations  (approximate ) 

pNLA <- function(probeID){

  condz <- c("min20","min40","min60","min80","hr2")

  ymax <- max(c(pam2.mat.max1[probeID,],lps.mat.max1[probeID,], pam3.mat.max1[probeID,],polyIC.mat.max1[probeID,], r848.mat.max1[probeID,], cpg.mat.max1[probeID,],polyIC.mat.max1[probeID,]))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Nuc. Loc.",sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mat.max1[probeID,1:6], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(condz)+1),labels=c("min0",condz))
  axis(2)
  
  lines( pam2.mat.max1[probeID,], col='green',lwd=3)
  lines( pam3.mat.max1[probeID,], col='red',lwd=3)
  lines( polyIC.mat.max1[probeID,], col='magenta',lwd=3)
  lines( r848.mat.max1[probeID,], col='purple',lwd=3)
  lines( cpg.mat.max1[probeID,], col='yellow',lwd=3)

  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG"),
       col=c('blue','green','red','magenta','purple','yellow'),
       lty=c(1,1,1,1,1,1),
       lwd=3)
  
  
}



## Plot ratios for all condtions for a single probe
plotAffyRatios <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  ymax <- max(c(pam2.ratios[probeID,],lps.ratios[probeID,], pam3.ratios[probeID,],polyIC.ratios[probeID,], r848.ratios[probeID,], cpg.ratios[probeID,],polyIC.ratios[probeID,]))
  
  ymin <- min(c(pam2.ratios[probeID,],lps.ratios[probeID,], pam3.ratios[probeID,],polyIC.ratios[probeID,], r848.ratios[probeID,],cpg.ratios[probeID,],polyIC.ratios[probeID,]),-1.5)
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")
  
  x11()
  
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.ratios[probeID,], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Log10Ratio',axes=FALSE,main=main,lwd=3)
  

  axis(1, 1:length(conds),labels=conds)
  axis(2)
  
  lines( pam2.ratios[probeID,], col='green',lwd=3)
  lines( pam3.ratios[probeID,], col='red',lwd=3)
  lines( polyIC.ratios[probeID,], col='magenta',lwd=3)
  lines( r848.ratios[probeID,], col='purple',lwd=3)
  lines( cpg.ratios[probeID,], col='yellow',lwd=3)

  lines(c(1,3,5),myd88KOlps.ratios[probeID,], col='blue',lwd=3,lty=2)
  lines(c(1,3,5),myd88KOpam3.ratios[probeID,], col='red',lwd=3,lty=2)
  lines(c(1,3,5),myd88KOpolyIC.ratios[probeID,], col='magenta',lwd=3,lty=2)
  
  lines(c(1,3,5),pam3polyIC.ratios[probeID,], col='skyblue',lwd=3)
  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(5.5,1.0,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG",
              "MyD88KO LPS","MyD88KO PAM3","MyD88KO PolyIC","PAM3+PolyIC"),
       col=c('blue','green','red','magenta','purple','yellow','blue','red','magenta','skyblue'),
       lty=c(1,1,1,1,1,1,2,2,2,1),
       lwd=3)
  
  
}



#---------------------------------------------------------------##

plotLPSProfile <- function(probeID){

  conds <- c("min20","min40","min60","min80","hr2","hr4","hr8","hr24")

  main <- paste(as.character(commonName[probeID]),", Probe ID:",probeID,sep="")
  
  ##main <- paste(gid2gname(probeID),", Probe ID:",probeID,sep="")
  
  mews <- lps.mus[probeID,]
  conds <- names(mews)
  ymax <- max(mews)

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  plot( mews, ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  axis(1, 1:length(conds),labels=conds)
  axis(2)

  ##legend(length(conds)-2,ymax/4., legend=c("LPS"), col=c('blue'), lwd=3)
  
  par(op)

  
}

## plot predicted matrix for a probeID
## predMat must be 6 rows, one for eac stim
## 5 cols for times 20 to 120 mins
paiPredStim <- function(probeID,predMat,t.index.min=1){

  conds <- c("min20","min40","min60","min80","min120")

  ymax <- max(c(pam2.mus[probeID,],lps.mus[probeID,1:6], pam3.mus[probeID,],polyIC.mus[probeID,], r848.mus[probeID,], cpg.mus[probeID,],polyIC.mus[probeID,],pred))
  
  ymin <- 0
  
  main <- paste(cname.compare[probeID],", Probe ID:",probeID,sep="")

  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  par(op)

  plot( lps.mus[probeID,1:6], ylim=c(ymin,ymax), col='blue', type='l',xlab='Time', ylab='Intensity',axes=FALSE,main=main,lwd=3)
  
  axis(1, 1:(length(conds)+1),labels=c("min0",conds))
  axis(2)
  
  lines( pam2.mus[probeID,], col='green',lwd=3)
  lines( pam3.mus[probeID,], col='red',lwd=3)
  lines( polyIC.mus[probeID,], col='magenta',lwd=3)
  lines( r848.mus[probeID,], col='purple',lwd=3)
  lines( cpg.mus[probeID,], col='yellow',lwd=3)


  t.index.max <- 6
  xvals <- 1:t.index.max

  ##lines(xvals[(t.index.min+1):t.index.max],predMat["lps",], col='blue',lty=2,lwd=3)
  lines(xvals[(t.index.min+1):t.index.max],predMat["pam2",], col='green',lty=2,lwd=3)
  lines(xvals[(t.index.min+1):t.index.max],predMat["pam3",], col='red',lty=2,lwd=3)
  lines(xvals[(t.index.min+1):t.index.max],predMat["polyIC",], col='magenta',lty=2,lwd=3)
  lines(xvals[(t.index.min+1):t.index.max],predMat["r848",], col='purple',lty=2,lwd=3)
  lines(xvals[(t.index.min+1):t.index.max],predMat["cpg",], col='yellow',lty=2,lwd=3)

  
  op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=1.5)
  
  par(op)
  legend(1,0.95*ymax,legend=c("LPS","PAM2","PAM3","PolyIC","R848","CpG"),
       col=c('blue','green','red','magenta','purple','yellow'),
       lty=c(1,1,1,1,1,1),
       lwd=3)
  
  
}


quickGroupPlot <- function ( psois, main="" ){
  profileplotRatio(lps.ratios[psois,],cname.compare[psois],main=main)
##  profileplotIntensity(lps.ratios[psois,],cname.compare[psois],main=main)
}
quickGP <- quickGroupPlot
qgp <- quickGroupPlot


plotFit <- function(observed,predicted,main.title="",time.constant=999999,deltaT=0){

  times <- c(0,20,40,60,80,120) ## fixed for now, could use as an argument

  # ymax should bound the data above by integer * power of ten
  yy <- max(c(data.vals,pred))
  dec <- ceiling(log10(yy))-1
  ymax <- 10^dec*ceiling(yy/10^dec)
  
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  
  plot(times,observed, ylim=c(0,ymax), col='blue', type='l',xlab='Time', ylab='Intensity Gain',axes=TRUE,main=main.title,lwd=3)
  lines(times,predicted, col='green',lwd=3)

  lines(c(deltaT,deltaT),c(0,ymax),col='magenta',lwd=3) 
  
  if ( time.constant != 999999 ){
    lines(c(time.constant+deltaT,time.constant+deltaT),c(0,ymax),col='red',lwd=3)
    legend(80,ymax/4., legend=c("data","prediction","start time","start+induction times"), col=c('blue','green','magenta','red'), lwd=3)
  } else {
    legend(80,ymax/4., legend=c("data","prediction","start time"), col=c('blue','green','magenta'), lwd=3)
  }

  
  par(op)

}


## Quick version of profileplot.R
## For ratio data, the legend
## is assumed to be commonNames
## no title on plot
ratioProfile <- function( gids ){

  profileplot(ratios[gids,],gid2gname(gids),"")

}
