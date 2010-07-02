


##
## Requirements
##
##load("/Users/thorsson/tfinf/annotations/annotation.objects.RData")
load("/Users/thorsson/tfinf/derived_data/models.rmsf.10May2010.RData")
load("/Users/thorsson/tfinf/annotations/allMouseMappings.August2009.RData")
load("/Users/thorsson/tfinf/sequence_data/sigPairedSites.28Sep2009.RData")
##load("/Users/thorsson/tfinf/interaction_data/pdnaModels.10May2010.RData")
load("/Users/thorsson/tfinf/annotations/tteMaps.RData")
source("/Users/thorsson/tfinf/analysis/evidenceUtilities.R")

mods <- mods.enrp.dubs.01.ode.rmsf[1:3]
res <- evidenceForMods(mods,n.cands=2)

## Set of all models for multiple genes
mts <- mods.s44
psois <- names(mts)


## Gene index
i <- 2
psoi <- psois[i]

## Neighbor set
friends <- nbrs[[psoi]]
set <- c(psoi,friends)
profileplotRatio(lps.ratios[set,],cc[set],"")

## Write filenames for display with Alistairs script
writeFastaNbrs(psoi,"testfile")

## collect model evidence. All TF models for that gene
evidenceForMods( mts[i], n.cands=2 )

## Do we output pairs on disk for each model?

## Folders for each gene model ?



