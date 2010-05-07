
exp.dir <- file.path(Sys.getenv("AA"),"data/20100426.curated.3prime")
exp.dir <- file.path(Sys.getenv("TFINF"),"expression_data")

max.intensity <- apply(lps.mus,1,max) 
lps.mat.max1 <- lps.mus/max.intensity ## this one will be modified for nucloc
lps.mat.max1.exp <- lps.mus/max.intensity ## this one will retain expression vars 

##max.intensity <- apply(pam2.mus,1,max) 
pam2.mat.max1 <- pam2.mus/max.intensity

##max.intensity <- apply(pam3.mus,1,max) 
pam3.mat.max1 <- pam3.mus/max.intensity

##max.intensity <- apply(polyIC.mus,1,max) 
polyIC.mat.max1 <- polyIC.mus/max.intensity

##max.intensity <- apply(r848.mus,1,max) 
r848.mat.max1 <- r848.mus/max.intensity

##max.intensity <- apply(cpg.mus,1,max) 
cpg.mat.max1 <- cpg.mus/max.intensity

ofile <- paste(exp.dir,"scaled.mus.objects.RData",sep="/")
savelist <- c("max.intensity","lps.mat.max1","lps.mat.max1.exp","pam2.mat.max1","pam3.mat.max1","polyIC.mat.max1","r848.mat.max1","cpg.mat.max1" )

save(list=savelist, file=ofile)
