
## To be turned into a command line script
## inputs: ensids, Parsed.Scan filename, output RData filename

filename.subset <- "ensids"

ensids <- as.character(read.table(filename.subset)$V1)

load("Parsed.Scan.Results.allMouseMarch2008.dat")
## gives TRE.OUT

not.there <- setdiff(ensids,names(TRE.OUT))
if ( length(not.there) > 0 ){
  for ( i in 1:length(not.there)) {
    cat("Not found:", not.there[i],"\n")
  }
}

there <- intersect(ensids,names(TRE.OUT))

TRE.OUT.ensids <- TRE.OUT[there]

save(TRE.OUT.ensids,file="Parsed.Scan.ensids.RData")
