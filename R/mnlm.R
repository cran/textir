##### Estimation for Regularized Logit Multinomial Regression  ######

## Main function; most happens in .C 
mnlm <- function(counts, covars, normalize=TRUE, penalty=c(shape=1,rate=1/2), start=NULL,
                 tol=1e-2, bins=0, verb=FALSE, ...)
{
  
  on.exit(.C("mnlm_cleanup", PACKAGE = "textir"))
  
  ## check and possibly bin observations for fast inference
  checked <- mncheck(counts=counts, covars=covars, normalize=normalize, bins=bins, penalty=penalty, verb=verb, ...)
  X <- checked$X # stm; adds a null column onto X if ncol(X) > 2
  V <- checked$V # stm; adds an intercept onto the covariates

  counts = checked$counts # simple triplet form
  never = checked$never
  fullnames = checked$fullnames
  covars = checked$covars # possibly normalized
  covarMean = checked$covarMean
  covarSD = checked$covarSD
  prior = checked$prior
  
  p <- ncol(X)  
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
      if(p<=10000){ phi <- suppressWarnings(try(corr(F,as.matrix(V[,-1])), silent=TRUE))
                    phi[is.nan(phi)] <- 0 }
      if(p > 10000 || class(phi)=="try-error"){
        phi <- tcrossprod_simple_triplet_matrix(t(F),t(as.matrix(V)[,-1,drop=FALSE]))
      }
    }
    q <- col_means(F)
    start <- cbind( log(q)-log(mean(q)), phi )
    if(verb){ cat("done.\n") }
  }
  else{
    if(ncol(start) != length(fullnames) || nrow(start) != d) stop("bad starting coefficient matrix")
    start <- t(start) ## transpose to match coef
    start[is.nan(start)] <- 0
    if(any(never)) start <- start[!never,] 
  }
  coef <- start
  if(p==2) coef[1,] <- rep(0,ncol(coef))

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
            qn = as.integer(checked$QN),
            verb = as.integer(verb),
            PACKAGE="textir")
  
  ## format output
  coef <- matrix(map$coef, ncol=d)
  map$fitted[is.nan(map$fitted)|!is.finite(map$fitted)|is.na(map$fitted)] <- 1/ncol(X)
  if(p==2){
    xhat <- as.matrix(simple_triplet_matrix(i=X$i, j=X$j, v=map$fitted, dimnames=dimnames(X)))
    xhat[xhat[,2]==0,2] <- 1 - xhat[xhat[,2]==0,1]
    xhat <- xhat[,2] 
    intercept <- matrix(coef[-1,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]][-1],"intercept"))
    loadings <- matrix(coef[-1,-1], nrow=1, dimnames=list(category=dimnames(counts)[[2]][-1], covariate=dimnames(covars)[[2]])) }
  else{
    xhat <- simple_triplet_matrix(i=X$i, j=X$j, v=map$fitted*m[X$i], dimnames=dimnames(X))
    intercept <- matrix(coef[,1], ncol=1, dimnames=list(category=dimnames(counts)[[2]],"intercept"))
    loadings <- matrix(coef[,-1], nrow=p, dimnames=list(category=dimnames(counts)[[2]], covariate=dimnames(covars)[[2]])) }
  loadings <- as.simple_triplet_matrix(loadings)

  ## clean un-needed summaries
  if(bins==0){ X=NULL; V=NULL; }

  ## calculate the estimated penalties 
  rateincrease = 2^map$lampar[1]
  if(rateincrease > 1){ cat(sprintf("WARNING: Hyperprior parameters were multiplied by %g to ensure global convergence.\n", rateincrease)) }

  ## deal with never observed categories
  if(any(never)){
    loadings <- simple_triplet_matrix(i=match(dimnames(loadings)[[1]][loadings$i],fullnames),
                                      j=loadings$j, v = loadings$v, nrow=length(fullnames),
                                      dimnames=list(category=fullnames, covariate=dimnames(loadings)[[2]]))
    fullalpha <- rep(log(1/(10+sum(counts))),length(fullnames))
    fullalpha[!never] <- intercept
    intercept <- matrix(fullalpha, ncol=1, dimnames=list(category=fullnames,"intercept"))
    counts <- simple_triplet_matrix(i=counts$i, j=match(dimnames(counts)[[2]][counts$j],fullnames), v=counts$v,
                                    ncol=length(fullnames), dimnames=list(dimnames(counts)[[1]],fullnames))
    xhat <- simple_triplet_matrix(i=xhat$i, j=match(dimnames(xhat)[[2]][xhat$j],fullnames), v=xhat$v,
                                    ncol=length(fullnames), dimnames=list(I=dimnames(xhat)[[1]],cat=fullnames))
    if(!is.null(X)){
      X <- simple_triplet_matrix(i=X$i, j=match(dimnames(X)[[2]][X$j],fullnames), v=X$v,
                                 ncol=length(fullnames), dimnames=list(I=dimnames(X)[[1]],cat=fullnames)) }
  }
    

  ## construct the return object
  out <- list(intercept=intercept, loadings=loadings,
              X=X, counts=counts, V=V, covars=covars, 
              normalized=normalize, binned=bins>0,
              covarMean=covarMean, covarSD=covarSD, 
              prior=penalty, fitted=xhat)
  
  ## class and return
  class(out) <- "mnlm"
  invisible(out)
}

## s3 plot method for mnlm
plot.mnlm<- function(x, type=c("response","reduction","roc"), covar=NULL, v=NULL, xlab=NULL, ylab=NULL, col=NULL, ...)
{
  if(type[1]=="roc"){
    if(max(x$counts$v) > 1) stop("ROC curves are provided only for max(count)=1 classification data")
    
    pred <- predict(x, newdata=x$covars)
    if(is.null(col) || length(col) != ncol(pred)) col = rainbow(ncol(pred))
    if(is.null(xlab)) xlab = "1-Specificity"
    if(is.null(ylab)) ylab = "Sensitivity"
    
    roc <- function(y,p,add,lc, ...){
      Q <- p > matrix(rep(sort(p),length(p)),ncol=length(p),byrow=TRUE)
      fp <- colSums((y==levels(y)[1])*Q)/sum(y==levels(y)[1])
      tp <- colSums((y==levels(y)[2])*Q)/sum(y==levels(y)[2])
      if(!add) plot(fp, tp, type="l", xlab=xlab, ylab=ylab, col=lc, ...)
      else lines(fp, tp, col=lc) }
    
    if(ncol(pred)==1) y <- factor(as.matrix(x$counts[,2])[,1]>0) 
    else  y <- factor(as.matrix(x$counts[,1])[,1]>0) 
    roc(y, pred[,1], FALSE, col[1], ...)
    abline(a=0,b=1,lty=2,col=8)
    if(ncol(pred)>1){
      for(j in 2:ncol(pred)) roc( factor(as.matrix(x$counts[,j])[,1]>0), pred[,j], TRUE, col[j], ...)         
    }
  }
  else if(type[1]=="reduction"){
    if(ncol(x$counts)<=2){
      cat("No useful sufficient reductions for a binary response.  Use type=`repsonse' instead.\n") 
      return(invisible()); }

    if(is.null(covar)){ covar <- 1 }
    if(is.null(col)){ col=1 }
    
    X <- as.simple_triplet_matrix(x$counts)
    frq <- freq(X)
    if(inherits(frq,"simple_triplet_matrix")){ z <- tcrossprod_simple_triplet_matrix(frq, t(as.matrix(x$loadings[,covar]))) }
    else{ z <- frq%*%as.matrix(x$loadings[,covar]) }
    
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
        if(is.null(col)){ col=rainbow(nlevels(counts)) }
        numcnt <- as.numeric(levels(counts))
        plot(x$fitted$v ~ counts, at=numcnt, xlab=xlab, ylab=ylab, col=col, varwidth=TRUE, xlim=range(numcnt)+c(-.25,.25), ... )    
      }
    }
  }
}

## S3 method coef function
coef.mnlm <- function(object, origscale=TRUE, ...){
  if(origscale && !is.null(object$covarSD)){
    object$loadings$v = object$loadings$v/object$covarSD[object$loadings$j]
    object$intercept <- object$intercept - col_sums(object$covarMean*t(object$loadings))
  } 
  return( as.matrix(t(cbind(object$intercept, object$loadings))) )
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
    
    phi <- as.matrix(object$loadings)
    if(nrow(phi) != ncol(F)){ stop("Dimension mismatch: nrow(phi) != ncol(X)") }
    
    return(tcrossprod_simple_triplet_matrix(F, t(phi)))
  }
  else{
    if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
    if(object$normalized){ newdata <- normalize(newdata, m=object$covarMean, s=object$covarSD) }
    if(ncol(newdata)!=ncol(object$loadings)){ stop("newdata must be a matrix with the same columns as object$covars") }
    newdata <- as.simple_triplet_matrix(newdata)

    expeta <- exp(tcrossprod_simple_triplet_matrix(cbind(rep(1,nrow(newdata)),newdata),
                                                   cbind(object$intercept,as.matrix(object$loadings))))
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

  ps <- round(100*sum(object$loadings==0)/(nrow(object$loadings)*ncol(object$loadings)),1) 
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
      cuts <- c(.1,.25,.5,.75,.9)
      F <- matrix(nrow=2,ncol=5, dimnames=list(c("sensitivity %","specificity %"),
                                   " classification cut-off" = cuts)) 
      for(i in 1:5){
        q <- object$fitted >= cuts[i];
        F[1,i] <- sum(as.matrix(object$counts)[q,2])/sum(as.matrix(object$counts)[,2])
        F[2,i] <- sum(as.matrix(object$counts)[!q,1])/sum(as.matrix(object$counts)[,1]) }
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
mncheck <- function(counts, covars, normalize, bins, penalty, verb, delta=1, dmin=0.01, randeffects=FALSE, quasinewton=TRUE){

  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix or a factor)
  if(is.null(dim(counts))){ 
    counts <- factor(counts)
    counts <- simple_triplet_matrix(i=1:length(counts), j=as.numeric(counts), v=rep(1,length(counts)),
                                    dimnames=list(NULL, response=levels(counts))) }
  if(inherits(counts, "TermDocumentMatrix")){ counts <- t(counts) }
  counts <- as.simple_triplet_matrix(counts)
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("cat",1:ncol(counts),sep="") }

  ## deal with never observed categories
  fullnames <- dimnames(counts)[[2]]
  never <- sdev(counts)==0
  counts <- counts[,!never]
  
  ## covariates input and check
  if(nrow(counts) != nrow(as.matrix(covars))){ stop("different number of predictor and response rows") }
  if(is.null(dim(covars))){
    covars <- matrix(covars, dimnames=list(dimnames(counts)[[1]], deparse(substitute(covars)))) }
  if(is.data.frame(covars)){ covars <- sapply(covars,as.numeric) }
  
  ## standardize for mean 0 variance 1.
  covarMean <- NULL
  covarSD <- sdev(covars)
  if(normalize)
    {
      if(is.simple_triplet_matrix(covars)){ covarMean <- 0 }
      else{
        covars <- as.matrix(covars)
        covarMean <- colMeans(covars) }
      covars <- normalize(covars, m=covarMean, s=covarSD)
    }

  ## bin and organize V
  X <- counts
  V <- covars

  if(bins<=1){  
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
    Xs <- rollup(X,1,I)

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

  ## build prior
  prior <- buildprior(penalty, ncol(V), ncol(X), verb=verb, fix=c(FALSE, covarSD==0))
  if(!normalize) covarSD <- NULL

  return( list(counts=counts, covars=covars, covarMean=covarMean, covarSD=covarSD, fullnames=fullnames, never=never,
               prior=prior, X=X, V=V, I=I, N=N, delta=delta, dmin=dmin, RE=randeffects, QN=quasinewton) )

}

## build the penalties
buildprior <- function(penalty, d, p, verb=0, fix){
  if(!is.list(penalty)){
    penlist <- vector(mode='list', length=d)
    if(p==2){ penlist[[1]] <- c(0,0) } # unpenalized
    else{ penlist[[1]] <- c(0,1) }# Norm(0,1) on the intercept
    penlist[2:d] <- rep(data.frame(penalty),d-1) 
  } else{
    if(length(penalty) != d){ stop("bad penalty argument") }
    penlist <- penalty
  }
  maplam = c()
  lampar = c()
  printem <- verb && d<5

  for(i in 1:d){
    if(length(penlist[[i]]) == 1) penlist[[i]] <- c(penlist[[i]],0)
    if(fix[i]) penlist[[i]] <- c(-1,0)
    if(length(penlist[[i]]) == 2){
      if(any(penlist[[i]]<0)){ ## fixed coef
        maplam = c(maplam,-1)
        lampar = c(lampar, c(0,0))
        if(printem) cat(sprintf("fixed coefficients for column %d\n", i-1)) }
      else if(penlist[[i]][2]==0){ ## L1 penalty
        maplam = c(maplam,0)
        lampar = c(lampar, penlist[[i]])
        if(printem){
          if(penlist[[i]][1]==0) cat(sprintf("No prior penalty for column %d\n", i-1))
          else cat(sprintf("laplace(%g) prior for column %d\n", penlist[[i]][1], i-1)) } }
      else if(penlist[[i]][1]==0){ ## L2 penalty
        maplam = c(maplam,2)
        lampar = c(lampar, penlist[[i]])
        if(printem) cat(sprintf("normal(0,%g) prior for column %d\n", 1/penlist[[i]][2], i-1)) }
      else{ ## gamma lasso penalty
        maplam = c(maplam,1)
        lampar = c(lampar, penlist[[i]])
        if(printem) cat(sprintf("gamma(%g,%g)-lasso prior for column %d\n",
                                penlist[[i]][1], penlist[[i]][2], i-1)) }
    }
    else{ stop("bad argument length within your penalty list") }
  }
  return(list(maplam=maplam, lampar=lampar)) }
  
