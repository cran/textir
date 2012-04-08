##### Estimation for Regularized Logit Multinomial Regression  ######

## Main function; most happens in .C 
mnlm <- function(counts, covars, normalize=TRUE, penalty=c(shape=1,rate=1/2), start=NULL,
                 tol=1e-2, bins=0, verb=FALSE, ...)
{
  
  on.exit(.C("mnlm_cleanup", PACKAGE = "textir"))
  
  ## check and possibly bin observations for fast inference
  checked <- mncheck(counts=counts, covars=covars, normalize=normalize, bins=bins, verb=verb, ...)
  X <- checked$X # stm; adds a null column onto X if ncol(X) > 2
  V <- checked$V # stm; adds an intercept onto the covariates

  counts = checked$counts # simple triplet form
  covars = checked$covars # possibly normalized
  covarMean = checked$covarMean
  covarSD = checked$covarSD
  
  p <- ncol(X)-1  
  n <- nrow(X)
  m <- row_sums(X)
  d <- ncol(V)
  
  ## initialization for the coefficient estimates
  if(is.null(start)){
    if(verb){ cat("Finding initial coefficient estimates... ") }
    F <- freq(X)
    phi <- NULL
    if(d >= 50){ phi <- matrix(0, nrow=p, ncol=d-1) }
    else{
      if(p<=10000){ phi <- suppressWarnings(try(corr(F[,-1],as.matrix(V[,-1])), silent=TRUE)) }
      if(p > 10000 || class(phi)=="try-error"){
        phi <- tcrossprod_simple_triplet_matrix(t(F[,-1]),t(as.matrix(V)[,-1,drop=FALSE]))
      }
    }
    q <- col_means(F)
    start <- cbind( log(q[-1])-log(max(q[1],0.0001)), phi)
    if(verb){ cat("done.\n") }
  }
  else{ if(nrow(start) != p || ncol(start) != d) stop("bad starting coefficient matrix") }
  coef <- rbind(rep(0,d), start)

  if(d>p) G <- t(-tcrossprod_simple_triplet_matrix(t(V), as.matrix(t(X))))
  else G <- -tcrossprod_simple_triplet_matrix(t(X), as.matrix(t(V)))

  ## RE create remean and revar
  if(checked$RE){
    rev <- log((checked$N+1)/checked$N)
    rem <- log(checked$N) - 0.5*rev
    xi <- c(0,cumsum(row_sums(X!=0)))
    revec = c(rem,rev)
  }
  else{ revec <- rev <- rem <- xi <- NULL }

  prior = buildprior(penalty, d, verb)
  
  map <- .C("Rmnlogit",
            n = as.integer(n),
            p = as.integer(p),
            d = as.integer(d),
            m = as.double(m),
            tol = as.double(tol),
            Nx = as.integer(length(X$v)),
            X = as.double(X$v),
            xind = as.integer(cbind(X$i,X$j)-1),
            xi = as.integer(xi),
            Nv = as.integer(length(V$v)),
            V = as.double(V$v),
            vind = as.integer(cbind(V$i,V$j)-1),
            vk = as.integer(c(0,cumsum(col_sums(V!=0)))),
            coef = as.double(coef),
            fitted = double(length(X$v)),
            maplam = as.integer(prior$maplam),
            lampar = as.double(prior$lampar),
            dmin = as.double(checked$dmin),
            delta = as.double(checked$delta),
            G = as.double(G),
            RE = as.integer(checked$RE),
            revec = as.double(revec),
            verb = as.integer(verb),
            PACKAGE="textir")
            

  ## drop the null category and format output
  coef <- matrix(map$coef, ncol=d)
  map$fitted[is.nan(map$fitted)|!is.finite(map$fitted)|is.na(map$fitted)] <- 1/ncol(X)
  
  if(ncol(X)>2){
    if(sum(X[,1])>0){ m <- m-X[,1]$v }
    xhat <- simple_triplet_matrix(i=X$i, j=X$j, v=map$fitted*m[X$i], dimnames=dimnames(X))[,-1]
    X <- X[,-1]
    intercept <- matrix(coef[-1,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]],"intercept"))
    loadings <- matrix(coef[-1,-1], nrow=p, dimnames=list(category=dimnames(counts)[[2]], covariate=dimnames(covars)[[2]]))
  } else{
    xhat <- as.matrix(simple_triplet_matrix(i=X$i, j=X$j, v=map$fitted, dimnames=dimnames(X)))
    xhat[xhat[,2]==0,2] <- 1 - xhat[xhat[,2]==0,1]
    xhat <- xhat[,2] 
    intercept <- matrix(coef[-1,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]][-1],"intercept"))
    loadings <- matrix(coef[-1,-1], nrow=p, dimnames=list(category=dimnames(counts)[[2]][-1], covariate=dimnames(covars)[[2]]))
    p <- p+1 }

  ## clean un-needed summaries
  if(bins==0){ X=NULL; V=NULL; }

  ## calculate the estimated penalties 
  rateincrease = 2.0^map$lampar[1]
  if(rateincrease > 1){ cat(sprintf("WARNING: Prior shape and rate were multiplied by %g in order to get convergence.\n", rateincrease)) }
  lambda = abs(loadings)
  if(is.list(penalty)){
    for(k in 2:d){
      if(length(penalty[[k]])>1){ penalty[[k]] <- penalty[[k]]*rateincrease
                                  lambda[,k-1] <- penalty[[k]][1]/(penalty[[k]][2] + lambda[,k-1]) }
      else{ lambda[,k-1] <- penalty[[k]] } } }
  else{
    if(length(penalty)>1){ penalty = penalty*rateincrease
                           lambda <- penalty[1]/(penalty[2] + lambda) }
    else{ lambda <- penalty } }

  ## construct the return object
  out <- list(intercept=intercept, loadings=loadings,
              X=X, counts=counts, V=V, covars=covars, 
              normalized=normalize, binned=bins>0,
              covarMean=covarMean, covarSD=covarSD, 
              prior=penalty, lambda=lambda,
              fitted=xhat)
  
  ## class and return
  class(out) <- "mnlm"
  invisible(out)
}

## s3 plot method for mnlm
plot.mnlm<- function(x, type=c("response","reduction"), covar=NULL, v=NULL, xlab=NULL, ylab=NULL, col=NULL, ...)
{
  if(type[1]=="reduction"){
    if(ncol(x$counts)<=2){
      cat("No useful sufficient reductions for a binary response.  Use type=`repsonse' instead.\n") 
      return(invisible()); }

    if(is.null(covar)){ covar <- 1 }
    if(is.null(col)){ col=1 }
    
    X <- as.simple_triplet_matrix(x$counts)
    frq <- freq(X)
    if(inherits(frq,"simple_triplet_matrix")){ z <- tcrossprod_simple_triplet_matrix(frq, t(x$loadings[,covar])) }
    else{ z <- frq%*%x$loadings[,covar] }
    
    if(is.null(v) || length(v)!=nrow(X)){
      if(is.null(xlab)){ xlab <- paste("covariate ",covar) }
      v <- as.matrix(x$covars[,covar]) }
    else if(is.null(xlab)){ xlab <- "V" }
    
    if(is.null(ylab)){ ylab <- paste("fitted direction",covar) }
    
    plot(z ~ v, xlab=xlab, ylab=ylab, col=col, ...)
    legend("topleft", legend=paste("corr =", round(cor(as.numeric(v), z),2)), bty="n", cex=1.2)
  }
  else{
    if(is.vector(x$fitted)){
        if(is.null(xlab)){ xlab <- "response" }
        if(is.null(ylab)){ ylab <- "fitted probability" }
        if(is.null(col)){ col=c(2,4) }
        if(!is.null(x$X)){
          resp = as.matrix(x$X[,2])/row_sums(x$X)
          plot(x$fitted ~ resp, bg=8, pch=21, cex=2, xlab=xlab, ylab=ylab, ...) }
        else{ plot(x$fitted ~ factor(as.matrix(x$counts[,2])), xlab=xlab, ylab=ylab, col=col, varwidth=TRUE, ... ) }
      }
    else{
      if(max(x$counts)==1 && is.null(x$X)){
        if(is.null(xlab)){ xlab <- "response" }
        if(is.null(ylab)){ ylab <- "fitted probability" }
        if(is.null(col)){ col=rainbow(ncol(x$counts)) }
        resp <- factor(dimnames(x$counts)[[2]][x$counts$j])
        plot(x$fitted$v ~ resp, xlab=xlab, ylab=ylab, col=col, varwidth=TRUE, ... )
      }
      else{
        if(is.null(xlab)){ xlab <- "observed count" }
        if(is.null(ylab)){ ylab <- "fitted count" }
        if(x$binned){ counts <- factor(x$X$v) }
        else{ counts <- factor(x$counts$v) }
        isfull <- counts%in%levels(counts)[table(counts)>=5]
        if(is.null(col)){ col=rainbow(nlevels(counts)) }
        numcnt <- as.numeric(levels(counts))
        plot(x$fitted$v ~ counts, at=numcnt, xlab=xlab, ylab=ylab, col=col, varwidth=TRUE, xlim=range(numcnt)+c(-.25,.25), ... )    
        points(as.numeric(as.character(counts[!isfull])), x$fitted$v[!isfull])
      }
    }
  }
}

 ## S3 method predict function
predict.mnlm <- function(object, newdata, type=c("response","reduction"), ...)
{
  if(type[1]=="reduction"){
    if(ncol(object$counts)<=2){
      cat("No useful sufficient reductions for a binary response.  Use type=`response' instead.\n") 
      return(invisible()); }
    if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
    F <- freq(as.simple_triplet_matrix(newdata))
    
    if(class(object)!="mnlm"){ stop("object class must be `mnlm'.") }
    
    phi <- object$loadings
    if(nrow(phi) != ncol(F)){ stop("Dimension mismatch: nrow(phi) != ncol(X)") }
    
    return(tcrossprod_simple_triplet_matrix(F, t(phi)))
  }
  else{
    if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
    if(object$normalized){ newdata <- normalize(newdata, m=object$covarMean, s=object$covarSD) }
    if(ncol(newdata)!=ncol(object$loadings)){ stop("newdata must be a matrix with the same columns as object$covars") }
    newdata <- as.simple_triplet_matrix(newdata)

    expeta <- exp(tcrossprod_simple_triplet_matrix(cbind(rep(1,nrow(newdata)),newdata),cbind(object$intercept,object$loadings)))
    expeta[!is.finite(expeta)] <- max(expeta[is.finite(expeta)])
    if(ncol(expeta)==1){ P <- expeta/(1+expeta) }
    else{ P <- expeta/row_sums(expeta) }
    dimnames(P) <- list(dimnames(newdata)[[1]], probability=dimnames(object$loadings)[[1]])
    return(P)
  }
}


 ## S3 method summary function
summary.mnlm <- function(object, y=NULL, ...){

  print(object)

  ps <- round(100*sum(object$loadings==0)/length(object$loadings),1) 
  cat(paste("   Loadings matrix is ",
            ps, "% sparse.\n\n", sep=""))
  
  if(ncol(object$counts)>2){
    z <- predict(object, newdata=object$counts, type="reduction") 
    if(!is.null(y)){
      reg <- lm(y~z)
      cat(paste("   Sufficent reduction R2 for y: ", round(cor(reg$fitted,y)^2,3),
                " (residual scale of ", round(sd(reg$resid),3), ")\n", sep="")) }

    if(ncol(object$covars)<=10){
      cat(paste("   Correlations in each sufficient reduction direction:\n     "))
      vars <- dimnames(object$covars)[[2]]
      if(is.null(vars)){ vars <- paste("var",1:ncol(z), sep="") }
      for( i in 1:ncol(z)){
        cat(paste(vars[i], ": ", round(cor(z[,i], as.matrix(object$covars[,i])),3), ". ", sep="")) }
      cat("\n")
    }
    cat("\n")
  }
  else if(is.null(object$X)){
    if(is.vector(object$fitted)){
      cat(paste("   False positives and negatives:\n "))
      cuts <- c(.1,.25,.5,.75,.9)
      F <- matrix(nrow=2,ncol=5, dimnames=list(c("    %fp","    %fn"),
                                   " classification cut-off" = cuts)) 
      for(i in 1:5){
        p <- object$fitted >= cuts[i];
        q <- object$fitted < cuts[i];
        F[1,i] <- sum( as.matrix(object$counts)[p,1])/sum(p)
        F[2,i] <- sum( as.matrix(object$counts)[q,2])/sum(q) }
      F[is.na(F)] <- 0
      print(round(F*100,1))
      cat("\n")
    } else{
      cat("   cor(counts,fitted) =",paste(round(cor(object$fitted$v, object$counts$v),2),
                "for nonzero entries.\n\n")) }
  }
}

print.mnlm <- function(x, ...){
  if(ncol(x$counts)>2){ nr = paste(nrow(x$loadings),"response categories") }
  else{ nr = "binary response" }
  
  cat(paste("\n   mnlm object for", nr, "and", ncol(x$loadings), "covariates.\n\n")) }

## quadratic function solver: y = x^2 + bx + c
quadratic <- function(b, c, quiet=FALSE, plot=FALSE)
  {
    soln <- .C("Rquadratic", coef = as.double(c(b,c)), num = double(1), PACKAGE="textir")
    roots <- soln$coef
    if(soln$num == 1)
      { if(!quiet){ cat( paste( "\n Single root at x =", round(roots[1],5),"\n\n")) }
        type="one real" }
    if(soln$num == 2)
      { if(!quiet){ cat( paste( "\n Two real roots at x1 =",
                   round(roots[1],5), "x2 =", round(roots[2],5), "\n\n") ) }
        type="two real" }
    if(soln$num == 0)
      { if(!quiet){ cat( paste( "\n Complex roots x =",  round(roots[1],5), "+-", round(roots[2],5), "i\n\n") ) }
        type= "two complex" }

    if(plot)
      { bnd <- mean(abs(c(b,c)))/10
        if(type=="two complex"){  x <- seq( roots[1] - abs(roots[2]) - bnd, roots[1] + abs(roots[2]) + bnd, length=100 ) }
        else{ x <- seq(min(roots)-bnd,max(roots)+bnd,length=1000) }
        plot(x, x^2 + b*x + c, type="l", main=type,ylab=paste("x^2 + ",b,"x + ",c,sep=""))
        abline(h=0, col=8)
        if(type=="two real"){ points(roots,rep(0,2), pch=21, bg=grey(.5)) }
        else{ points(roots[1],0, pch=21, bg=grey(.5)) } }

    return(list(type=type, coef=c(b,c), roots=roots)) 

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

#######  Undocumented "mnlm"-related utility functions #########

## check the inputs and bin
mncheck <- function(counts, covars, normalize, bins, verb, delta=1, dmin=0.01, nullfactor=0, randeffects=FALSE){

  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix or a factor)
  if(is.null(dim(counts))){ 
    counts <- factor(counts)
    counts <- simple_triplet_matrix(i=1:length(counts), j=as.numeric(counts), v=rep(1,length(counts)),
                                    dimnames=list(NULL, response=levels(counts))) }
  if(inherits(counts, "TermDocumentMatrix")){ counts <- t(counts) }
  counts <- as.simple_triplet_matrix(counts)

  ## covariates input and check
  if(nrow(counts) != nrow(as.matrix(covars))){ stop("different number of predictor and response rows") }
  if(is.null(dim(covars))){
    covars <- matrix(covars, dimnames=list(dimnames(counts)[[1]], deparse(substitute(covars)))) }
  
  ## standardize for mean 0 variance 1.
  covarMean <- covarSD <- NULL
  if(normalize)
    {
      if(is.simple_triplet_matrix(covars)){ covarMean <- 0 }
      else{
        covars <- as.matrix(covars)
        covarMean <- colMeans(covars) }
      
      covarSD <- sdev(covars)
      if(any(covarSD==0)){ warning("You have a column of 'covars' with sd=0; do not include an intercept term.")
                           covarSD[covarSD==0] <- 1  }
      covars <- normalize(covars, m=covarMean, s=covarSD)
    }

  ## bin, add the null category, and organize V
  X <- counts
  V <- covars


  if(nullfactor < 0){ stop("You need a positive nullfactor") }

  if(bins<=1){  
    if(ncol(X)>2){ X <- cbind(as.vector(row_sums(X)*nullfactor), X) }
    X <- X
    V <- cbind(rep(1,nrow(V)),V)
    N <- rep(1,nrow(V))
    I <- NULL
  }
  else{
    V <- as.matrix(V)  #  no point in sparse format when binning
    R <- apply(V, 2, range)
    B <- mapply(seq, from=R[1,], to=R[2,], MoreArgs=list(length=bins+1))[-1,]
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

    Vm <- apply(V, 2, function(v) tapply(v, I, mean) )

    xij <- interaction(as.numeric(I)[X$i], X$j, drop=TRUE)
    vals <- tapply(X$v, xij, sum) 
    ij <- matrix(as.numeric(unlist(strsplit(names(vals), split="\\."))), ncol=2, byrow=TRUE)
    Xs <- simple_triplet_matrix(i=ij[,1], j=ij[,2], v=vals, dimnames=list(I=dimnames(Vm)[[1]], cat=dimnames(X)[[2]]))
    if(ncol(Xs)>2){ Xs <- cbind(as.vector(row_sums(Xs)*nullfactor),Xs) } 

    if(verb){
      cat("\nResponse Bins: \n")
      print(Vm)
      cat("\n") }


    X = Xs
    V = cbind(1,Vm)
    N = as.vector(table(I))
    I = as.numeric(I)
  }

  ## stm for V
  V <- as.simple_triplet_matrix(V)

  ## format V with col major index ordering
  o <- order(V$j)
  V$i <- V$i[o]
  V$j <- V$j[o]
  V$v <- V$v[o]
  
  ## format X row major index ordering
  if(randeffects){
    o <- order(X$i)
    X$i <- X$i[o]
    X$j <- X$j[o]
    X$v <- X$v[o] }

  return( list(counts=counts, covars=covars, covarMean=covarMean, covarSD=covarSD, X=X, V=V, I=I, N=N, delta=delta, dmin=dmin, RE=randeffects) )

}

## build the penalties
buildprior <- function(penalty, d, verb=0){
  if(is.list(penalty)){
    if(length(penalty) != d){ stop("bad penalty argument") }
    maplam = c()
    lampar = c()
    for(i in 1:d){
      if(length(penalty[[i]]) == 1){
        if(penalty[[i]] < 0){
          maplam = c(maplam,-1)
          lampar = c(lampar, c(0,0))
          if(verb && d<10) cat(sprintf("Coefficients taken as given for column %d\n", i-1))
        }
        else{
          maplam = c(maplam,0)
          lampar = c(lampar, c(penalty[[i]], 0))
          if(verb && d<10) cat(sprintf("Fixed L1 penalty at %g for column %d\n", penalty[[i]], i-1))
        }
      }
      else if(length(penalty[[i]]) == 2){
        if(any(penalty[[i]]<=0)) stop("Gamma hyperprior requires parameters > 0")
        maplam = c(maplam,1)
        lampar = c(lampar, c(penalty[[i]][1], penalty[[i]][2]))
        if(verb && d<10) cat(sprintf("Gamma(%g,%g) prior on penalties for column %d\n", penalty[[i]][1], penalty[[i]][2], i-1))
      }
      else{ stop("bad lambda argument within your penalty list") }
    }
  }
  else{
    if(length(penalty)==1){
      maplam <- rep(0,d)
      lampar <- c(c(0,0),rep(c(penalty,0),d-1)) 
      if(verb) cat(sprintf("Fixed L1 penalty of %g.\n", penalty))
    }
    else if(length(penalty)==2){
      if(any(penalty<=0)) stop("Gamma hyperprior requires parameters > 0")
      maplam = c(0,rep(1,d-1))
      lampar = c(c(0,0), rep(c(penalty[1], penalty[2]),d-1))
      if(verb) cat(sprintf("Gamma(%g,%g) prior on penalties.\n", penalty[1], penalty[2]))
    }
    else{ stop("bad lambda argument: length(penalty) > 2") }
  }
  return(list(maplam=maplam, lampar=lampar)) }
  
