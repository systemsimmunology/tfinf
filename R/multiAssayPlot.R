
## Feb 25, 2006
## Plot multiple assays for TFs

## Assume we have (normalized to 1)
## protein.cytoplasmic.western ; t.western
## protein.nuclear.western ; t.western
## protein.il6.bound ; t.chip
## mRNA.realtime ; t.realtime
## mRNA.array ; t.array

"multiAssayPlot" <- function(tf,t.max){

  colorvec = c('magenta','green','blue','turquoise','red')
  op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5,mar=(c(5, 4, 6, 5)+ 0.1))
  ymax <- 1
  ymin <- 0
  main <- tf
  ##op <- par(font=2,lwd=2,font.axis=2,font.main=2,font.lab=2,font.sub=2,cex=1.5)
  par(op)

  t.all.range <- (t.all <= t.max)
  t.western.range <- (t.western <= t.max)
  t.chip.range <- (t.chip <= t.max)
  t.realtime.range <- (t.realtime <= t.max)
  t.array.range <- (t.array <= t.max)
  
  plot(t.all[t.all.range],xlim=c(0,t.max),ylim=c(ymin,ymax),xlab='Time (Minutes)', main=tf,ylab='normalized units',axes=FALSE,type="n")
  axis(1,at=t.all[t.all.range])
  axis(2)

  have.data <- rep(FALSE,5)
  
  ## western, cytoplasmic
  if ( tf %in% rownames(protein.cytoplasmic.western) ){
    have.data[1] <- TRUE
    lines( t.western[t.western.range],protein.cytoplasmic.western[tf,t.western.range], col=colorvec[1],lty=1,type='l',lwd=3)
  }

  ## western, nuclear
  if ( tf %in% rownames(protein.nuclear.western) ){
    have.data[2] <- TRUE
    lines( t.western[t.western.range],protein.nuclear.western[tf,t.western.range], col=colorvec[2],lty=1,type='l',lwd=3)
  }

  ## chIP
  if ( tf %in% rownames(protein.il6.bound) ){
    have.data[3] <- TRUE
    lines( t.chip[t.chip.range],protein.il6.bound[tf,t.chip.range], col=colorvec[3],lty=1,type='l',lwd=3)
  }

  ## mRNA realtime
  if ( tf %in% rownames(mRNA.realtime) ){
    have.data[4] <- TRUE
    lines( t.realtime[t.realtime.range],mRNA.realtime[tf,t.realtime.range], col=colorvec[4],lty=1,type='l',lwd=3)
  }

  ## mRNA array
  if ( tf %in% rownames(mRNA.array) ){
    have.data[5] <- TRUE
    lines( t.array[t.array.range],mRNA.array[tf,t.array.range], col=colorvec[5],lty=1,type='l',lwd=3)
  }

  ##op <- par(font=1,lwd=2,font.axis=2,font.main=1,font.lab=1,font.sub=1,cex=0.6)
  ##par(op)
  
  legend(0.3*t.max,0.95*ymax,legend=c("Cyt. Prot.","Nuc. Prot","IL6 chr. IP","mRNA RT","mRNA array")[have.data],col=colorvec[have.data],lwd=3)
  
}
