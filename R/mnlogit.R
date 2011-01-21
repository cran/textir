
## laplace regularized logit multinomial regression 

mnlm <- function(counts, covars, normalize=FALSE, lambda=NULL, start=NULL,
                 tol=0.1, tmax=1000, delta=1, dmin=0, verb=TRUE)
{  
  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  X <- as.simple_triplet_matrix(counts)
  p <- ncol(X)-1  
  n <- nrow(X)
  m <- row_sums(X)
  
  ## covariates input, check, and possibly standardize for mean 0 variance 1.
  covars <- as.matrix(covars)
  if(n != nrow(covars)){ stop("different number of predictor and response rows") }
  covarMean <- covarSD <- NULL
  if(normalize)
    { covarMean <- colMeans(covars)
      covarSD <- apply(covar,2,sd)
      covar <- normalize(covar, m=covarMean, s=covarSD)
    }
  V <- cbind(rep(1,n), covars)
  d <- ncol(V)

  ## initialization for the coefficient estimates
  if(is.null(start)){
    cat("Finding initial coefficient estimates... ")
    F <- freq(counts)
    phi <- as.matrix(cor(F[,-1],covars))
    q <- colMeans(F)
    start <- cbind( log(q)-log(q[1]), rbind(rep(0,d-1),phi))
    cat("done.\n")
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
  out <- list(coef=coef, L=L, niter=niter, X=X, V=V, covarMean=covarMean, covarSD=covarSD,
              maplam=maplam, lampar=lampar, lambda=lambda, delta=delta,
              tol=tol, tmax=tmax, start=start)
  
  ## class and return
  class(out) <- "mnlm"
  invisible(out)
}
