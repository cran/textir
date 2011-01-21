##  normalizing design matrices
normalize <- function(x, m=NULL, s=NULL){
  x <- as.matrix(x)
  if(is.null(m)){ m <- apply(x,2,mean) }
  if(is.null(s)){ s <- apply(x,2,sd) }
  return( t((t(x) - m)/s) ) }

## converting count to frequency matrix
freq <- function(x){
  if(class(x)[1] == "TermDocumentMatrix"){ x <- t(x) }
  x <- as.matrix(x)
  return( x/rowSums(x) ) }

## partial least squares
pls <- function(F, y, K=1){
  F <- as.matrix(F)
  phi <- matrix(ncol=K, nrow=ncol(F))
  z <- matrix(ncol=K, nrow=nrow(F))
  yhat <- NULL
  v <- matrix(y)
  cat("Directions ")
  for(k in 1:K)
    { cat(paste(k,", ",sep=""))
      phi[,k] <- cor(F,v[,k]) #inverse regression, equiv: t(lm(F~v[,k])$coef)[,2] 
      z[,k] <- F%*%phi[,k] # fitted direction
      v <- cbind(v, lm(v[,k] ~ z[,k])$resid) # orthogonalize
      yhat <- cbind(yhat, lm(y~z[,1:k])$fitted)
    }
  
  cat("done.\n")
  out <- list(y=y, F=F, z=z, phi=phi, v=v[,-K], yhat = yhat, fwdmod = lm(y~z) )
  class(out) <- "pls"
  return(out) }

