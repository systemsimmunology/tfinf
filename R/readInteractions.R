## August 30, 2005
##
## readInteractions.R
##
##
## No special treatment of redundant interactions
## 18 Jan 2006. Simplify to get rid of time-consuming loop

"readInteractions" <- function (filename){
  
  all.ints <- read.table(filename, sep=' ', header=FALSE)
  
  mem1s <- as.character(all.ints$V1)
  intModes <- as.character(all.ints$V2)
  mem2s <- as.character(all.ints$V3)


  keepers <- (mem1s %in% np) & ( mem2s %in% np )
  not.same <- !(mem1s==mem2s)
  keepers <- keepers & not.same

  returnList <- list(partner1=mem1s[keepers],partner2=mem2s[keepers],interactType=intModes[keepers])
  returnList

  ##counter <- 0
  ##partner1 <- character()
  ##partner2 <- character()
  ##interactType <- character()  

  ## Can't quite recall why this filter is in place
  ## It does seem that if the interactor is not in the np vector, there is no representative probe
  ## Human NPs e.g. ?

  ## This is really slow. Replace with vector evaluations 

  ##for (i in 1:length(mem1s) ) {
  ##  intMode <- intModes[i]
  ##  mem1 = mem1s[i]
  ##  if(mem1 %in% np ){
  ##    mem2 = mem2s[i]
  ##    if(mem2 %in% np ){
  ##      if ( mem1!=mem2 ){# ditch the self-loops
  ##        counter <- counter+1
  ##        partner1[counter] <- mem1
  ##        partner2[counter] <- mem2
  ##        interactType[counter] <- intMode
  ##      }
  ##    }
  ##  }
  ##}

  ##returnList <- list(partner1=partner1,partner2=partner2,interactType=interactType)
  ##returnList
  
}
  
