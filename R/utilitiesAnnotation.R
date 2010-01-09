

# more efficient version of the old function
gid2gname <- function ( gids ){
  as.character(cname.compare[gids])
  ##  as.character(commonName[gids])
}


## Since this can be one to many, allow only single as input
gname2gid <- function ( gname ){
  names(which(cname.compare==gname))
##  names(which(commonName==gname))
}

# could make the usage more consitent here, and above
cnameFromProbeSet <- function ( psvec ){
  as.character(commonName[psvec])
}
 
npFromProbeSet <- function ( psvec ){
  as.character(np[psvec])
}

probesetFromNP <- function( singleNP ){
  names(which(np==singleNP))
}
 
gidFromProbeSet <- function ( psvec ){
  as.character(gid[psvec])
}
 
cnameFromNP <- function ( npvec ) {
  returnVec <- character()

  for ( npc in npvec ){
    gg <- gid2gname(names(which(np==npc)))
    gg <- gg[1] ## can have multiple hits. They *should* all be the same
    returnVec <- c(returnVec,gg)
  }

  returnVec
}


npFromCname <- function( cvec ){
  returnVec <- character()
  for ( cname in cvec ){
    gg <- np[names(which(cname.compare==cname))]
    gg <- as.character(gg[1]) ## can have multiple hits. They *should* all be the same
    returnVec <- c(returnVec,gg)
  }
  returnVec
}
