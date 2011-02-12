## minimalist partial least squares
pls <- function(F, y, K=1, scale=TRUE, verb=TRUE){

  if(scale){
      if(!inherits(F, "simple_triplet_matrix"))
        { F <- t(t(F)/apply(F,2,sd)) } else{ F$v <-  F$v/sdev(F)[F$j] }
    }
 
  phi <- matrix(ncol=K, nrow=ncol(F))
  z <- matrix(ncol=K, nrow=nrow(F))
  yhat <- NULL

  if(ncol(as.matrix(y)) > 1){stop( "PLS only works for univariate y (vector or single column matrix).")}
  v <- normalize(as.numeric(y))
  
  if(verb){ cat("Directions ") }
  for(k in 1:K){
    if(verb){ cat(paste(k,", ",sep="")) }

    ## inverse regression, equiv: t(lm(F~v[,k])$coef)[,2]
    phi[,k] <- corr(F,v[,k])
    
    ## project the fitted direction
    if(inherits(F, "simple_triplet_matrix")){
      z[,k] <- tcrossprod_simple_triplet_matrix(F, t(phi[,k]))
    } else { z[,k] <- F%*%phi[,k] }

    ## orthogonalize
    if(k<K){ v <- cbind(v, lm(v[,k] ~ z[,k])$resid) }
    yhat <- cbind(yhat, lm(as.numeric(y)~z[,1:k])$fitted)
  }  
  if(verb){ cat("done.\n")}
  
  out <- list(y=y, F=F, z=z, phi=phi, v=v, yhat = yhat, fwdmod = lm(as.numeric(y)~z), scale=scale)
  class(out) <- "pls"
  return(out) }

## converting count to frequency matrix
freq <- function(x, byrow=TRUE){
  if(inherits(x, "simple_triplet_matrix")){
    f <- x
    if(byrow){ f$v <- x$v/row_sums(x)[x$i] }
    else{ f$v <- x$v/col_sums(x)[x$i] }
    return(f)
  } else{
    if(byrow){ return( x/rowSums(x) ) }
    else{ return(t(t(x)/colSums(x))) }
  } }

## converting a count/freq matrix to tfidf
tfidf <- function(x, freq=FALSE){

  if(!freq){tf <- freq(x)}
  else{tf <- x}
  
  if(inherits(x, "simple_triplet_matrix")){ x$v <- 1 }
  else{ x[x>0] <- 1 }

  idf <- log( nrow(x) ) - log(1+col_sums(x>0))

  if(inherits(x, "simple_triplet_matrix")){ x$v <- tf$v*idf[x$j] }
  else{ x <- t(t(tf)*idf) }

  return( x ) }

## correlation for slam simple_triplet_matrix and regular matrix
corr <- function(x, y){
  if(!inherits(x, "simple_triplet_matrix")){ return(cor(x,y) ) }

  s <- sdev(x)
  z <-  simple_triplet_matrix(i=x$i, j=x$j, v=x$v/s[x$j])
  v <- t(normalize(y))
  
  r <- tcrossprod_simple_triplet_matrix(t(z), v)/(nrow(x)-1)
  dimnames(r) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
  return( r ) }
  
## column standard deviation for a simple_triplet_matrix 
sdev <- function(x){
  if(!inherits(x, "simple_triplet_matrix")){ return(apply(x,2,sd)) }
  n <- nrow(x)
  x2 <-  simple_triplet_matrix(i=x$i, j=x$j,  v=x$v^2)
  return( sqrt(col_sums(x2)/(n-1) - col_sums(x)^2/(n^2 - n)) ) }

##  normalizing design matrices
normalize <- function(x, m=NULL, s=NULL){
  x <- as.matrix(x)
  if(is.null(m)){ m <- apply(x,2,mean) }
  if(is.null(s)){ s <- apply(x,2,sd) }
  return( t((t(x) - m)/s) ) }
