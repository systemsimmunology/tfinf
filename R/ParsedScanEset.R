
seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")

load(paste(annot.dir,"expressed.scanned.ensembl.RData",sep="/"))
load(paste(seq.dir,"Parsed.Scan.Results.allMouseAug2010.dat",sep="/"))

TRE.OUT.eset <- TRE.OUT[expressed.scanned.ensembl]

save(TRE.OUT.eset,file=paste(seq.dir,"Parsed.Scan.eset.RData",sep="/"))
