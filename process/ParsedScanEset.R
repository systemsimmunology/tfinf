out.dir <- file.path(Sys.getenv("TFINF"),"data/process_output")
seq.dir <- file.path(Sys.getenv("TFINF"),"data/motifscan_input")

load(paste(out.dir,"expressed.scanned.ensembl.RData",sep="/"))
load(paste(seq.dir,"Parsed.Scan.Results.dat",sep="/"))

## 10-11-12 Don't know why, but a single gene :ensembl ID 73406 is not in the scanset:
## ENSMUSG00000073406	14963
## 14963	H2-Bl	histocompatibility 2, blastocyst

TRE.OUT.eset <- TRE.OUT[intersect(expressed.scanned.ensembl,names(TRE.OUT))]
save(TRE.OUT.eset,file=paste(out.dir,"Parsed.Scan.eset.RData",sep="/"))

# hate to do this, but here goes
expressed.scanned.ensembl <- names(TRE.OUT.eset)
save(expressed.scanned.ensembl,file=paste(out.dir,"expressed.scanned.ensembl.RData",sep="/"))

