#### plot functions for textir ####

plot.mnlm<- function(x, covar=NULL, y=NULL, cat=NULL, delta=0.05,...)
{
  if(is.null(covar)){ covar <- 1 }
  if(is.null(cat)){
    X <- x$X
    frq <- freq(X)
    if(inherits(frq,"simple_triplet_matrix")){ z <- tcrossprod_simple_triplet_matrix(frq, t(x$coef[,covar+1])) }
    else{ z <- frq%*%x$coef[,covar+1] }
    
    if(is.null(y) || length(y)!=nrow(X)){  y <- x$V[,covar] }
    
    plot(z ~ y, 
         xlab=paste("response factor",covar),
         ylab=paste("fitted dimension",covar), ...)
    legend("topleft", legend=paste("R2 =", round(cor(x$V[,covar], z)^2,2)), bty="n", cex=1.2)
  }
  else{  coefplot(x, cat, covar+1, delta, ...) }
}
  

plot.pls <- function(x, K=NULL, ...){
  if(is.null(K)){ K <- 1:ncol(x$yhat) }
  y <- x$y 
  K <- as.vector(K)
  par(mfrow=c(1,length(K)))
  for(k in K){
    plot(x$yhat[,k] ~y, xlab="response", ylab=paste("pls(",k,") fitted values"), ...)
    legend("topleft", legend=paste("R2 =",
                        round(cor(as.numeric(y),x$yhat[,k])^2,2)), bty="n")
  }
}
  

## utility function for coefficient plots
coefplot <- function(fit, cat, covar,
                     delta=NULL, s=NULL, r=NULL,
                     from=NULL, to=NULL, length=200, ...)
{
  B <- fit$coef
  V <- cbind(1, fit$V)
  X <- as.simple_triplet_matrix(fit$X)
  p <- nrow(B)-1
  K <- ncol(B)-1
  m <- row_sums(X)
  if(fit$maplam){ if(is.null(s)){ s <- fit$lampar[1] }
                  if(is.null(r)){ r <- fit$lampar[2] } }
  if(is.null(delta)){ delta <- fit$delta[cat,covar] }
  
  # find the converged lhd
  eta <- V%*%t(B)
  denom <- rowSums(exp(eta))
  Bsum <- sum(abs(B[,-1]))
  l <- -sum(X$v*eta[(X$j-1)*nrow(X) + X$i]) + sum(m*log(denom))
  if(fit$maplam){
    c <- (s+p*K-1)*Bsum/(r+Bsum)
  } else { c <- Bsum*fit$lambda }
  
  g <- -sum(V[,covar]*(as.matrix(X[,cat]) - m*exp(eta[,cat])/denom))
 
  E <- denom - exp(eta[,cat])
  e <- E
  e[E<exp(eta[,cat] - delta)] <- exp(eta[,cat] - delta)[E<exp(eta[,cat] - delta)]
  e[E>exp(eta[,cat] + delta)] <- exp(eta[,cat] + delta)[E>exp(eta[,cat] + delta)]

  F <- e/E + E/e + 2
  H <- sum(V[,covar]^2*m/F)

  # grid along our chosen phi_catcovar
  if(is.null(from)){ from <- B[cat,covar] - 1.2*delta }
  if(is.null(to)){ to <- B[cat,covar] + 1.2*delta }

  b <- seq(from, to, length=length)
  Bcat <- matrix(rep(B[cat,],length),nrow=K+1)
  Bcat[covar,] <- b
  
  l1dif <- -sum(as.matrix(X[,cat])*V[,covar])*(b - B[cat,covar])
  denomdif <- exp(V%*%Bcat) - rep(exp(V%*%B[cat,]), length)
  l2dif <- colSums(m*(log(denom+denomdif)-log(denom)))

  Bsumvec <- Bsum + (covar>1)*(abs(b)-abs(B[cat,covar]))
  if(fit$maplam){
    cgrid <- (s+p*K-1)*Bsumvec/(r+Bsumvec)
  } else{ cgrid <- fit$lambda*Bsumvec }

  L <- l + l1dif + l2dif + cgrid
  Bnd <- l + g*(b-B[cat,covar]) + 0.5*H*(b-B[cat,covar])^2  + cgrid
  
  plot(b, L-(l+c), type="l",
       xlab=paste("phi [", dimnames(X)[[2]][cat], "]"),
       ylab="L( phi ) - solved objective", ...)
  points(B[cat,covar], 0, pch=20, cex=1.5, ...)
  lines(b, Bnd-(l+c), lty=2, ...)

}
    
    
