## Feb 2009

load("models.rmsf.16Jan2009.RData")

## Contains c("mods.curated.ode.rmsf","mods.hs.et.05.ode.rmsf","mods.hs.et.01.ode.rmsf","mods.hs.et.001.ode.rmsf","mods.enrs.001.ode.rmsf","mods.singles.fromCorTFPairs.01.ode.rmsf","mods.enrp.sings.01.ode.rmsf","mods.enrp.dubs.01.ode.rmsf")


###
### June 2008
### 

singles.tags <- c("mods.curated","mods.hs.et.05","mods.enrs.001","mods.singles.fromCorTFPairs.01","mods.enrp.sings.01")
doubles.tags <- "mods.enrp.dubs.01"
model.tags <- c(singles.tags,doubles.tags)
modnames <- paste(model.tags,".ode.rmsf",sep="")

sto  <- function(inString){
 return (eval(parse(text=inString)))
}
modvec.1 <- createModVecFromMultMod ( modnames[1], sto(modnames[1]), n.cands=1)
modvec.2 <- createModVecFromMultMod ( modnames[2], sto(modnames[2]), n.cands=1)
modvec.3 <- createModVecFromMultMod ( modnames[3], sto(modnames[3]), n.cands=1)
modvec.4 <- createModVecFromMultMod ( modnames[4], sto(modnames[4]), n.cands=1)
modvec.5 <- createModVecFromMultMod ( modnames[5], sto(modnames[5]), n.cands=1)
modvec.6 <- createModVecFromMultMod ( modnames[6], sto(modnames[6]), n.cands=2)

modvec.now <- mergeModVecs(modvec.1,modvec.2)
modvec.now <- mergeModVecs(modvec.now,modvec.3)
modvec.now <- mergeModVecs(modvec.now,modvec.4)
modvec.now <- mergeModVecs(modvec.now,modvec.5)
modvec.now <- mergeModVecs(modvec.now,modvec.6)
modvec.now <- mergeModVecs(modvec.now,modvec.7)
modvec.e  <- modvec.now
rm(modvec.now)

save(modvec.e,file="modvec.e.RData")
