
## Return named row, cols, if applicable
## No check for 'collapse to vector'
which.max.mat <- function ( inMat ) {
  ind1 <- which.max(inMat) ## R returns a 'vector-type' index
  n.row <- dim(inMat)[1]
  col.ind <- ind1 %/% n.row + 1
  row.ind <- ind1 %% n.row
  if ( row.ind == 0 ){
    col.ind <- col.ind - 1 
    row.ind <- n.row
  }
  if ( inMat[row.ind,col.ind] != max(inMat) ){
    stop("Wrong index found")
  }
  ret1 <- row.ind
  ret2 <- col.ind
  if ( !is.null(dimnames(inMat)[[1]]) ){
    ret1 <- dimnames(inMat)[[1]][row.ind]
  }
  if ( !is.null(dimnames(inMat)[[2]]) ){
    ret2 <- dimnames(inMat)[[2]][col.ind]
  }
  return ( c(ret1,ret2) )
}

## Return named row, cols, if applicable
## No check for 'collapse to vector'
which.min.mat <- function ( inMat ) {
  ind1 <- which.min(inMat) ## R returns a 'vector-type' index
  n.row <- dim(inMat)[1]
  col.ind <- ind1 %/% n.row + 1
  row.ind <- ind1 %% n.row
  if ( row.ind == 0 ){
    col.ind <- col.ind - 1 
    row.ind <- n.row
  }
  if ( inMat[row.ind,col.ind] != min(inMat) ){
    stop("Wrong index found")
  }
  ret1 <- row.ind
  ret2 <- col.ind
  if ( !is.null(dimnames(inMat)[[1]]) ){
    ret1 <- dimnames(inMat)[[1]][row.ind]
  }
  if ( !is.null(dimnames(inMat)[[2]]) ){
    ret2 <- dimnames(inMat)[[2]][col.ind]
  }
  return ( c(ret1,ret2) )
}



## Minimum for three dimensional matrix
## Return named row, cols, if applicable
## No check for 'collapse to vector'
which.min.3mat <- function ( in3Mat ) {
  ind <- which.min(in3Mat) ## R returns a 'vector-type' index

  dims <- dim(in3Mat)
  n1 <- dims[1]
  n2 <- dims[2]
  n3 <- dims[3]

  n1.ind <- ind %% n1
  remainder <-  ind %/% n1 + 1
  n2.ind <- remainder %% n2
  n3.ind <- remainder %/% n2 + 1


  ## Correction needed for 'mod=0' cases
  ## Remain to be fully tested
  ## 2 matrix code below, for comparison
  if ( n1.ind == 0 ) {
    n1.ind <- n1
    n2.ind <- n2.ind - 1
  }
  if ( n2.ind == 0 ) {
    n2.ind <- n2
    n3.ind <- n3.ind - 1
  }
  ##  n.row <- dim(inMat)[1]
  ##col.ind <- ind1 %/% n.row + 1
  ##row.ind <- ind1 %% n.row
  ##if ( row.ind == 0 ){
  ##  col.ind <- col.ind - 1 
  ##  row.ind <- n.row
  ##}

  if ( in3Mat[n1.ind,n2.ind,n3.ind] != min(in3Mat) ){
    stop("Wrong index found")
  }
  ret1 <- n1.ind
  ret2 <- n2.ind
  ret3 <- n3.ind
  if ( !is.null(dimnames(in3Mat)[[1]]) ){
    ret1 <- dimnames(in3Mat)[[1]][n1.ind]
  }
  if ( !is.null(dimnames(in3Mat)[[2]]) ){
    ret2 <- dimnames(in3Mat)[[2]][n2.ind]
  }
  if ( !is.null(dimnames(in3Mat)[[3]]) ){
    ret3 <- dimnames(in3Mat)[[3]][n3.ind]
  }
  return ( c(ret1,ret2,ret3) )
}


##
## Testers 
## 

## Square matrix 
## m1 <- matrix(sample(1:100,9),nrow=3,ncol=3)

## Rectangular 
## m2 <- matrix(sample(1:100,12),nrow=3,ncol=4)
