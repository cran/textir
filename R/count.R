## Tools for manipulation of text count matrices ##

## converting count to frequency matrix
freq <- function(x, byrow=TRUE){
    if(byrow){ return( x/row_sums(x) ) }
    else{ return(t(t(x)/col_sums(x))) }
}

## converting a count/freq matrix to tfidf
tfidf <- function(x, freq=FALSE){

  if(!freq){tf <- freq(x)}
  else{tf <- x}
  
  idf <- log( nrow(x) ) - log(col_sums(x>0) + 10^{-100}) # e-100 for non-inf with never occurring terms
  x <- t( t(tf) * idf )
  
  return( x ) }

## correlation for slam simple_triplet_matrix and regular matrix
corr <- function(x, y){
  if(!inherits(x, "simple_triplet_matrix")){ return(cor(x,y) ) }

  n <- nrow(x)
  v <- t(normalize(y))
  
  r <- tcrossprod_simple_triplet_matrix(t(x)/sdev(x), v)/(nrow(x)-1)
  dimnames(r) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
  return( r ) }
  
## column standard deviation for a simple_triplet_matrix 
sdev <- function(x){
  if(!inherits(x, "simple_triplet_matrix")){ return(apply(x,2,sd)) }
  n <- nrow(x)
  sqrt(col_sums(x^2)/(n-1) - col_sums(x)^2/(n^2 - n))
  return( sqrt(col_sums(x^2)/(n-1) - col_sums(x)^2/(n^2 - n)) ) }

##  normalizing design matrices
normalize <- function(x, m=NULL, s=NULL, undo=FALSE){
  x <- as.matrix(x)
  if(!undo){
    if(is.null(m)){ m <- apply(x,2,mean) }
    if(is.null(s)){ s <- apply(x,2,sd) }
    return( t((t(x) - m)/s) ) }
  else{
    if(is.null(m) || is.null(s)){ stop("can't un-normalize without mean and sd.") }
    return( t( t(x)*s + m ) )}
}
    
## Dirichlet RNG
rdir <- function(n, alpha)
{
    x <- matrix(rgamma(length(alpha)*n,alpha),nrow=n,byrow=TRUE)
    return(t(x/rowSums(x))) }
