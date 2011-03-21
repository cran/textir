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


## s3 plot method for pls
plot.pls <- function(x, K=NULL, xlab="response", ylab=NULL, ...){
  if(is.null(K)){ K <- 1:ncol(x$yhat) }
  y <- x$y 
  K <- as.vector(K)
  par(mfrow=c(1,length(K)))
  ylb <- ylab
  for(k in K){
    if(is.null(ylab)){ ylb=paste("pls(",k,") fitted values") }
    plot(x$yhat[,k] ~y, xlab=xlab, ylab=ylb, ...)
    legend("topleft", legend=paste("R2 =",
                        round(cor(as.numeric(y),x$yhat[,k])^2,2)), bty="n", cex=1.4)
  }
}
  
