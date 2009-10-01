
ddata.dir <- file.path(Sys.getenv("TFINF"),"derived_data")

load(file=paste(ddata.dir,"collection.1reg.RData",sep="/"))

long.vec <- as.vector(sort(collection.1reg))
n.samps <- length(long.vec)

rands.1reg.0.05 <- long.vec[round(n.samps/20)]
rands.1reg.0.10 <- long.vec[round(n.samps/10)]
rands.1reg.0.20 <- long.vec[round(n.samps/5)]
rands.1reg.0.50 <- long.vec[round(n.samps/2)]

load(file=paste(ddata.dir,"collection.2regs.RData",sep="/"))

long.vec <- as.vector(sort(collection.2regs))
n.samps <- length(long.vec)

rands.2regs.0.05 <- long.vec[round(n.samps/20)]
rands.2regs.0.10 <- long.vec[round(n.samps/10)]
rands.2regs.0.20 <- long.vec[round(n.samps/5)]
rands.2regs.0.50 <- long.vec[round(n.samps/2)]
