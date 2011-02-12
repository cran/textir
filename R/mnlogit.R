
## laplace regularized logit multinomial regression 

mnlm <- function(counts, covars, normalize=FALSE, lambda=NULL, start=NULL,
                 tol=0.1, tmax=1000, delta=1, dmin=0, bins=0, verb=TRUE)
{  
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

  ## possibly bin observations for fast inference
  binned <- txtbin(counts, covars, bins)
  X <- binned$X
  V <- cbind(1,binned$V)
  
  p <- ncol(X)-1  
  n <- nrow(X)
  m <- row_sums(X)
  d <- ncol(V)
  
  ## initialization for the coefficient estimates
  if(is.null(start)){
    if(verb){ cat("Finding initial coefficient estimates... ") }
    F <- freq(X)
    phi <- corr(F[,-1],V[,-1])
    q <- col_means(F)
    start <- cbind( log(q)-log(q[1]), rbind(rep(0,d-1),phi))
    if(verb){ cat("done.\n") }
  }
  else{ if(nrow(start) != p+1 || ncol(start) != d) stop("bad starting coefficient matrix") }
  coef <- start 

  ## the delta step size matrix
  Dmat <- matrix(rep(delta, d*(p+1)), ncol=d)
  Dmat[1,] <- 0.0;
  
  if(length(lambda)==1){ maplam <- 0 }
  else
    { maplam = 1
      if(is.null(lambda)){ lambda <- c(sqrt(2), p/4, p/4) }
      else if(length(lambda) != 3 || prod(lambda>0) == 0){ stop("bad lambda argument") } }
  lampar <- lambda

  map <- .C("cgd_mnlogit",
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
    lampar <- map$lampar[2:3]
  }
  else{ lampar <- NULL }

  ## construct the return object
  out <- list(coef=coef, L=L, niter=niter, X=counts, V=covars,
              covarMean=covarMean, covarSD=covarSD,
              maplam=maplam, lampar=lampar, lambda=lambda, delta=delta,
              tol=tol, tmax=tmax, start=start)
  
  ## class and return
  class(out) <- "mnlm"
  invisible(out)
}


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

txtbin <- function(X, V, bins){

  if(bins<=1){ return( list(X=X, V=V, I=NULL) ) }             

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
  Xs <- simple_triplet_matrix(i=ij[,1], j=ij[,2], v=vals, dimnames=list(I=dimnames(V)[[1]], cat=dimnames(X)[[2]]))

  return( list(X=Xs, V=Vm, I=as.numeric(I)) )

}
  
