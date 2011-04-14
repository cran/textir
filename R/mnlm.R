##### Estimation for Regularized Logit Multinomial Regression  ######

## Main function; most happens in .C 
mnlm <- function(counts, covars, normalize=FALSE, penalty=c(1,0.2), start=NULL,
                 tol=0.1, tmax=1000, delta=1, dmin=0, bins=0, verb=FALSE)
{
  
  on.exit(.C("mnlm_cleanup", PACKAGE = "textir"))
  
  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix or a factor)
  if(is.null(dim(counts))){ 
    counts <- as.factor(counts)
    counts <- simple_triplet_matrix(i=1:length(counts), j=as.numeric(counts), v=rep(1,length(counts)),
                                    dimnames=list(obsv=1:length(counts), response=levels(counts)),
                                    nrow=length(counts), ncol=nlevels(counts)) }
  if(inherits(counts, "TermDocumentMatrix")){ counts <- t(counts) }
  
  counts <- as.simple_triplet_matrix(counts)
  ## covariates input and check
  if(nrow(counts) != nrow(as.matrix(covars))){ stop("different number of predictor and response rows") }
  if(is.null(dim(covars))){
    covars <- matrix(covars, dimnames=list(dimnames(counts)[[1]], deparse(substitute(covars)))) }
  
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
  X <- binned$X # adds a null column onto X if ncol(X) > 2
  V <- cbind(1, binned$V)
  
  p <- ncol(X)-1  
  n <- nrow(X)
  m <- row_sums(X)
  d <- ncol(V)
  
  vmax <- apply(abs(V),2,max)
  bmax <- log(10^6)/vmax

  ## initialization for the coefficient estimates
  if(is.null(start)){
    if(verb){ cat("Finding initial coefficient estimates... ") }
    F <- freq(X)
    phi <- NULL
    if(d >= 100){ phi <- matrix(0, nrow=p, ncol=d-1) }
    else if(p<=10000){ phi <- suppressWarnings(try(corr(F[,-1],V[,-1]), silent=TRUE)) }
    if(d < 100 && p > 10000 || class(phi)=="try-error"){
      phi <- tcrossprod_simple_triplet_matrix(t(F[,-1]),t(matrix(V[,-1], ncol=d-1)))
    }
    q <- col_means(F)
    start <- cbind( log(q[-1])-log(q[1]), phi)
    if(verb){ cat("done.\n") }
  }
  else{ if(nrow(start) != p || ncol(start) != d) stop("bad starting coefficient matrix") }
  coef <- rbind(rep(0,d), start)
  
  if(length(penalty)==1){ maplam <- 0 }
  else
    { maplam = 1
      if(length(penalty) != 2 || prod(penalty>0) == 0){ stop("bad lambda argument") } }
  lampar <- penalty

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
            bmax = as.double(bmax),
            L = double(tmax + 1),
            resids = double(length(X$v)),
            xhat = double(length(X$v)),
            maplam = as.integer(maplam),
            lampar = as.double(lampar),
            dmin = as.double(dmin),
            delta = as.double(delta),
            verb = as.integer(verb),
            PACKAGE="textir")
            

  coef <- matrix(map$coef, ncol=d)
  L <- map$L[1:map$niter]
  resids <- simple_triplet_matrix(i=X$i, j=X$j, v=map$resids, dimnames=dimnames(X))
  xhat <- simple_triplet_matrix(i=X$i, j=X$j, v=map$xhat, dimnames=dimnames(X))

  if(ncol(X)>2){
    resids <- resids[,-1]
    xhat <- xhat[,-1]*1.001
    X <- X[,-1]
    intercept <- matrix(coef[-1,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]],"intercept"))
    loadings <- matrix(coef[-1,-1], nrow=p, dimnames=list(category=dimnames(counts)[[2]], covariate=dimnames(covars)[[2]]))
  } else{
    if(max(counts)==1){ xhat <- as.matrix(xhat)
                        xhat[xhat[,2]==0,2] <- 1 - xhat[xhat[,2]==0,1]
                        xhat <- xhat[,2] }
    intercept <- matrix(coef[-1,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]][-1],"intercept"))
    loadings <- matrix(coef[-1,-1], nrow=p, dimnames=list(category=dimnames(counts)[[2]][-1], covariate=dimnames(covars)[[2]]))
    p <- p+1 }
  
  if(bins==0){ X=NULL; V=NULL; }
  if(!normalize){ covarSD = NULL }

  ## construct the return object
  out <- list(intercept=intercept, loadings=loadings,
              penalty=penalty, normalized=normalize, binned=bins>0,
              X=X, counts=counts, V=V, covars=covars, covarMean=covarMean, covarSD=covarSD, 
              L=L, residuals=resids, fitted=xhat)
  
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
      v <- x$covars[,covar] }
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
        plot(x$fitted ~ factor(as.matrix(x$counts[,2])), xlab=xlab, ylab=ylab, col=col, varwidth=TRUE, ... )
      }
    else{
      if(max(x$counts)==1){
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
        plot(x$fitted$v ~ counts, xlab=xlab, ylab=ylab, col=col, xaxt="n", varwidth=TRUE,  ... )    
        ax <- axTicks(1)
        ax[1] = 1
        axis(1, at=ax)
        points(as.numeric(counts[!isfull]), x$fitted$v[!isfull])
      }
    }
  }
}

 ## S3 method predict function
predict.mnlm <- function(object, newdata, type=c("response","reduction"), ...)
{
  if(type[1]=="reduction"){
    if(ncol(object$counts)<=2){
      cat("No useful sufficient reductions for a binary response.  Use type=`repsonse' instead.\n") 
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
    newdata <- as.matrix(newdata)
    if(object$normalized){ newdata <- normalize(newdata, m=object$covarMean, s=object$covarSD) }
    if(ncol(newdata)!=ncol(object$loadings)){ stop("newdata must be a matrix with the same columns as object$covars") }
    
    eta <- cbind(rep(1,nrow(newdata)),newdata)%*%t(cbind(object$intercept,object$loadings))
    if(ncol(eta)==1){ P <- exp(eta)/(1+exp(eta)) }
    else{ P <- exp(eta)/row_sums(exp(eta)) }
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
    
    cat(paste("   Correlations in each sufficient reduction direction:\n     "))
    vars <- dimnames(object$covars)[[2]]
    if(is.null(vars)){ vars <- paste("var",1:ncol(z), sep="") }
    for( i in 1:ncol(z)){
      cat(paste(vars[i], ": ", round(cor(z[,i], object$covars[,i]),3), ". ", sep="")) }
    cat("\n\n")
  }
  else{
    if(is.vector(object$fitted)){
      cat(paste("   False positives and negatives:\n "))
      cuts <- c(.1,.25,.5,.75,.9)
      F <- matrix(nrow=2,ncol=5, dimnames=list(c("    %fp","    %fn"),
                                   " classification cut-off" = cuts)) 
      for(i in 1:5){
        p <- object$fitted >= cuts[i];
        q <- object$fitted < cuts[i];
        F[1,i] <- sum(as.matrix(object$counts)[p,1])/sum(p)
        F[2,i] <- sum(as.matrix(object$counts)[q,2])/sum(q) }
      F[is.na(F)] <- 0
      print(round(F*100,1))
      cat("\n")
    }
    else{
      cat("   cor(counts,fitted) =",paste(round(cor(object$fitted$v, object$counts$v),2),
                "for nonzero entries.\n\n"))
    }
  }
}

print.mnlm <- function(x, ...){
  if(ncol(x$counts)>2){ nr = paste(nrow(x$loadings),"response categories") }
  else{ nr = "binary response" }
  
  cat(paste("\n   mnlm object for", nr, "and", ncol(x$loadings), "covariates.\n\n")) }

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

mncheck <- function(X, V, bins){

  if(bins<=1){
    if(ncol(X)>2){ X <- cbind(as.vector(0.001*row_sums(X)/0.999), X) }
    return( list(X=X, V=V, I=NULL) ) }             

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
  if(ncol(Xs)>2){ Xs <- cbind(as.vector(0.001*row_sums(Xs)/0.999),Xs) }
  
  return( list(X=Xs, V=Vm, I=as.numeric(I)) )

}
