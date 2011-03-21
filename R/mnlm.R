##### Estimation for Regularized Logit Multinomial Regression  ######

## Main function; most happens in .C 
mnlm <- function(counts, covars, normalize=FALSE, lambda=NULL, start=NULL,
                 tol=0.1, tmax=1000, delta=1, dmin=0, bins=0, verb=TRUE)
{

  on.exit(.C("mnlm_cleanup", PACKAGE = "textir"))
  
  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
  if(inherits(counts, "TermDocumentMatrix")){ counts <- t(counts) }
  counts <- as.simple_triplet_matrix(counts)
  ## covariates input and check
  covars <- as.matrix(covars)
  if(nrow(counts) != nrow(covars)){ stop("different number of predictor and response rows") }
  
  ## standardize for mean 0 variance 1.
  covarSD <- apply(covars,2,sd)
  if(prod(covarSD!=0)!=1){ stop("You have a constant covariate; do not include an intercept term in 'covars'.") }
  covarMean <- NULL
  if(normalize)
    { covarMean <- colMeans(covars)
      covars <- normalize(covars, m=covarMean, s=covarSD)
    }

  ## check and possibly bin observations for fast inference
  binned <- mncheck(counts, covars, bins)
  X <- binned$X # adds a null column onto X
  V <- cbind(1,binned$V)
  
  p <- ncol(X)-1  
  n <- nrow(X)
  m <- row_sums(X)
  d <- ncol(V)
  
  ## initialization for the coefficient estimates
  if(is.null(start)){
    if(verb){ cat("Finding initial coefficient estimates... ") }
    F <- freq(X)
    phi <- NULL
    if(p<=10000){ phi <- suppressWarnings(try(corr(F[,-1],V[,-1]), silent=TRUE)) }
    if(p > 10000 || class(phi)=="try-error"){
      phi <- tcrossprod_simple_triplet_matrix(t(F[,-1]),t(matrix(V[,-1], ncol=d-1))) }
    q <- col_means(F)
    start <- cbind( log(q[-1])-log(q[1]), phi)
    if(verb){ cat("done.\n") }
  }
  else{ if(nrow(start) != p || ncol(start) != d) stop("bad starting coefficient matrix") }
  coef <- rbind(rep(0,d), start)

  ## the delta step size matrix
  Dmat <- matrix(rep(delta, d*(p+1)), ncol=d)
  Dmat[1,] <- 0.0;
  
  if(length(lambda)==1){ maplam <- 0 }
  else
    { maplam = 1
      if(is.null(lambda)){ lambda <- c(p/4, p/4) }
      else if(length(lambda) != 2 || prod(lambda>0) == 0){ stop("bad lambda argument") } }
  lampar <- lambda

  map <- .C("Rmnlogit",
            n = as.integer(n),
            p = as.integer(p),
            d = as.integer(d),
            m = as.integer(m),
            tol = as.double(tol),
            niter = as.integer(tmax),
            N = as.integer(length(X$v)),
            X = as.double(X$v),
            xi = as.integer(cbind(X$i,X$j)-1),
            V = as.double(V),
            coef = as.double(coef),
            L = double(tmax + 1),
            maplam = as.integer(maplam),
            lampar = as.double(lampar),
            lambda = double(maplam*(tmax+1)),
            dmin = as.double(dmin),
            delta = as.double(Dmat),
            verb = as.integer(verb),
            PACKAGE="textir")
            

  coef <- matrix(map$coef, ncol=d)
  delta <- matrix(map$delta, ncol=d)
  niter <- map$niter
  L <- map$L[1:niter]
  if(maplam){
    lambda <- map$lambda[1:niter]
    lampar <- map$lampar
  }
  else{ lampar <- NULL }

  ## construct the return object
  out <- list(intercept=coef[-1,1],
              loadings=matrix(coef[-1,-1], nrow=p, dimnames=list(phrase=dimnames(counts)[[2]], direction=dimnames(covars)[[2]])),
              X=X, counts=counts,
              covars=covars, V=V, covarMean=covarMean, covarSD=covarSD,
              maplam=maplam, lampar=lampar, lambda=lambda, 
              delta=delta, L=L, niter=niter, tol=tol, tmax=tmax, start=start)
  
  ## class and return
  class(out) <- "mnlm"
  invisible(out)
}

## s3 plot method for mnlm
plot.mnlm<- function(x, covar=NULL, v=NULL, cat=NULL, delta=0.05, xlab=NULL, ylab=NULL, ...)
{
  if(is.null(covar)){ covar <- 1 }
  if(is.null(cat)){
    
    X <- as.simple_triplet_matrix(x$counts)
    frq <- freq(X)
    if(inherits(frq,"simple_triplet_matrix")){ z <- tcrossprod_simple_triplet_matrix(frq, t(x$loadings[,covar])) }
    else{ z <- frq%*%x$loadings[,covar] }
    
    if(is.null(v) || length(v)!=nrow(X)){
      if(is.null(xlab)){ xlab <- paste("response ",covar) }
      v <- x$covars[,covar] }
    else if(is.null(xlab)){ xlab <- "V" }

    if(is.null(ylab)){ ylab <- paste("fitted direction",covar) }
    
    plot(z ~ v, xlab=xlab, ylab=ylab, ...)
    legend("topleft", legend=paste("corr =", round(cor(as.numeric(v), z),2)), bty="n", cex=1.2)
  }
  else{  coefplot(fit=x, cat=cat, covar=covar+1, xlab=xlab, ylab=ylab, delta=delta, ...) }
}

 ## S3 method predict function
predict.mnlm <- function(object, newcounts, ...)
{
  if(is.vector(newcounts)){ newcounts <- matrix(newcounts, nrow=1) }
  F <- freq(as.simple_triplet_matrix(newcounts))
  
  if(class(object)!="mnlm"){ stop("object class must be `mnlm'.") }

  phi <- object$loadings
  if(nrow(phi) != ncol(F)){ stop("Dimension mismatch: nrow(phi) != ncol(X)") }

  return(tcrossprod_simple_triplet_matrix(F, t(phi))) }


## cubic function solver: y = x^3 + ax^2 + bx + c
cubic <- function(a, b, c, quiet=FALSE, plot=FALSE)
  {
    soln <- .C("Rcubic", coef = as.double(c(a,b,c)), num = double(1), PACKAGE="textir")
    roots <- soln$coef
    if(soln$num == 1)
      { if(!quiet){ cat( paste( "\n Single root at x =", round(roots[1],5),"\n\n")) }
        type="one real" }
    if(soln$num == 3)
      { if(!quiet){ cat( paste( "\n Three real roots at x1 =",
                   round(roots[1],5), "x2 =", round(roots[2],5), "x3 =", round(roots[3],5), "\n\n") ) }
        type="three real" }
    if(soln$num == 2)
      { if(!quiet){ cat( paste( "\n One real root at x =", round(roots[1],5),
                   " and complex roots x =",  round(roots[2],5), "+-", round(roots[3],5), "i\n\n") ) }
        type= "one real, two complex" }

    if(plot)
      { bnd <- mean(abs(c(a,b,c)))/10
        if(type=="three real"){ x <- seq(min(roots)-bnd,max(roots)+bnd,length=1000) }
        else{ x <- seq(min(roots[1:2])-bnd,max(roots[1:2])+bnd,length=1000) }
        plot(x, x^3 + a*x^2 + b*x + c, type="l", main=type,ylab=paste("x^3 + ",a,"x^2 + ",b,"x + ",c,sep=""))
        abline(h=0, col=8)
        if(type=="three real"){ points(roots,rep(0,3), pch=21, bg=grey(.5)) }
        else{ points(roots[1],0, pch=21, bg=grey(.5)) } }

    return(list(type=type, coef=c(a,b,c), roots=roots)) 

  }

#######  Undocumented "mnlm" utility functions #########

## coefficient plots
coefplot <- function(fit, cat, covar, xlab=NULL, ylab=NULL,
                     delta=0.05, s=NULL, r=NULL,
                     from=NULL, to=NULL, length=100, ...)
{
  cat <- cat+1 # shift to account for null
  B <- cbind(fit$intercept, fit$loadings)
  B <- rbind(rep(0,ncol(B)),B)
  V <- fit$V
  X <- fit$X
  p <- nrow(B)-1
  K <- ncol(B)-1
  m <- row_sums(X)
  if(fit$maplam){ if(is.null(s)){ s <- fit$lampar[1] }
                  if(is.null(r)){ r <- fit$lampar[2] } }
  
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

  if(is.null(xlab)){ xlab <- paste("phi [", dimnames(X)[[2]][cat], "]") }
  if(is.null(ylab)){ ylab <- "L( phi ) - solved objective" }
  plot(b, L-(l+c), type="l",
       xlab=xlab, ylab=ylab, ...)
  points(B[cat,covar], 0, pch=20, cex=1.5, ...)
  lines(b, Bnd-(l+c), lty=2, ...)

}
   
mncheck <- function(X, V, bins){

  if(bins<=1){ return( list(X=cbind(rep(0.01,nrow(X)),X), V=V, I=NULL) ) }             

  V <- as.matrix(V)
  R <- apply(V, 2, range)
  B <- mapply(seq, from=R[1,], to=R[2,], MoreArgs=list(length=bins))
  O <- apply(V, 2, order)-1

  out  <- .C("Rbin",
             n = as.integer(nrow(V)),
             d = as.integer(ncol(V)),
             b = as.integer(bins),
             B = as.double(B),
             V = as.double(V),
             O = as.integer(O),
             PACKAGE="textir")
  
  F <- apply(matrix(out$O, nrow=nrow(V), ncol=ncol(V)), 2, as.factor)
  I <- interaction(as.data.frame(F), drop=TRUE)

  Im <- function(v){ return(tapply(v, I, mean)) }
  Vm <- apply(V,2,Im)

  xij <- interaction(as.numeric(I)[X$i], X$j, drop=TRUE)
  vals <- tapply(X$v, xij, sum) 
  ij <- matrix(as.numeric(unlist(strsplit(names(vals), split="\\."))), ncol=2, byrow=TRUE)
  Xs <- simple_triplet_matrix(i=ij[,1], j=ij[,2], v=vals, dimnames=list(I=dimnames(Vm)[[1]], cat=dimnames(X)[[2]]))
  Xs <- cbind(rep(0.01,nrow(Xs)),Xs)
  
  return( list(X=Xs, V=Vm, I=as.numeric(I)) )

}
  
