## minimalist partial least squares
pls <- function(X, y, K=1, scale=TRUE, verb=TRUE){

  if(scale){
    sdf = sdev(X)
    if(!inherits(X, "simple_triplet_matrix"))
      { X <- t(t(X)/sdf) } else{ X$v <-  X$v/sdf[X$j] }
    scale = sdf
  }
 
  phi <- matrix(ncol=K, nrow=ncol(X))
  z <- matrix(ncol=K, nrow=nrow(X))
  yhat <- NULL

  if(ncol(as.matrix(y)) > 1){stop( "PLS only works for univariate y (vector or single column matrix).")}
  v <- normalize(as.numeric(y))
  
  if(verb){ cat("Directions ") }
  for(k in 1:K){
    if(verb){ cat(paste(k,", ",sep="")) }

    ## inverse regression, equiv: t(lm(X~v[,k])$coef)[,2]
    phi[,k] <- corr(X,v[,k])
    
    ## project the fitted direction
    if(inherits(X, "simple_triplet_matrix")){
      z[,k] <- tcrossprod_simple_triplet_matrix(X, t(phi[,k]))
    } else { z[,k] <- X%*%phi[,k] }

    ## orthogonalize
    if(k<K){ v <- cbind(v, lm(v[,k] ~ z[,k])$resid) }
    yhat <- cbind(yhat, lm(as.numeric(y)~z[,1:k])$fitted)
  }  
  if(verb){ cat("done.\n")}
  
  fwdmod = lm(as.numeric(y)~z)

  dimnames(z) <- list(dimnames(X)[[1]], direction=paste("z",1:ncol(z),sep=""))
  dimnames(yhat) <- list(dimnames(X)[[1]], model=paste("pls",1:ncol(z),sep=""))
  dimnames(phi) <- list(dimnames(X)[[2]], factor=paste("v",1:ncol(z) -1 ,sep=""))
  dimnames(v) <- list(dimnames(X)[[1]], factor=paste("v",1:ncol(z) -1 ,sep=""))
  
  out <- list(y=y, X=X, directions=z, loadings=phi, factors=v, fitted = yhat, fwdmod = fwdmod, scale=scale)
  class(out) <- "pls"
  return(out) }


## S3 plot method 
plot.pls <- function(x, K=NULL, xlab="response", ylab=NULL, ...){
  if(is.null(K)){ K <- 1:ncol(x$fitted) }
  y <- x$y 
  K <- as.vector(K)
  par(mfrow=c(1,length(K)))
  ylb <- ylab
  for(k in K){
    if(is.null(ylab)){ ylb=paste("pls(",k,") fitted values") }
    plot(x$fitted[,k] ~y, xlab=xlab, ylab=ylb, ...)
    legend("topleft", legend=paste("corr =",
                        round(cor(as.numeric(y),x$fitted[,k]),2)), bty="n", cex=1.4)
  }
}
  

## S3 method summary function
summary.pls <- function(object, ...){
  print(object)
  cat("Forward regression summary:\n")
  print(summary(object$fwdmod))
}

## S3 method summary function
print.pls <- function(x, ...){
  cat(paste("\nA pls(", ncol(x$directions), ") object, reduced from ", ncol(x$X), " input variables. \n\n", sep="")) }

 ## S3 method predict function
predict.pls <- function(object, newdata, response=TRUE, ...)
{
  if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
  if(object$scale != 0){
    if(!inherits(newdata, "simple_triplet_matrix"))
      { newdata <- t(t(newdata)/object$scale) } else{ newdata$v <-  newdata$v/object$scale[newdata$j] }
  }
  if(inherits(newdata, "simple_triplet_matrix")){
    z <- tcrossprod_simple_triplet_matrix(newdata, t(object$loadings))
    } else { z <- newdata%*%object$loadings }
  if(response){
    fitted <- cbind(1,z)%*%object$fwdmod$coef
    return(fitted)
  } else{  return(z) }
}

  

