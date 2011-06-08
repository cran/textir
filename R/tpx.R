#######  Undocumented "tpx" utility functions #########

## ** Only referenced from topics.R

## theta initialization
tpxinit <- function(X, initheta, K1, alpha, verb, ...){

  if(is.matrix(initheta)){
    if(ncol(initheta)!=K1){ stop("mis-match between initheta and K.") } 
    return(freq(initheta + 0.01/ncol(X), byrow=FALSE)) }

  if(is.matrix(alpha)){
    if(nrow(alpha)!=ncol(X) || ncol(alpha)!=K1){ stop("bad matrix alpha dimensions; check your K") }
    return(freq(alpha, byrow=FALSE)) }

  ilength <- initheta[1]
  if(ilength < 1){ ilength <- 1 }

  if(verb){ cat("Initial topics for K = ") }
  nK <- length( Kseq <-  unique(ceiling(seq(2,K1,length=ilength))) )
  initheta <- tpxThetaStart(X, matrix(col_sums(X)/sum(X), ncol=1), matrix(rep(1,nrow(X))), 2)
  
  ## loop over topic numbers
  for(i in 1:nK){

    ## Solve for map omega in NEF space
    fit <- tpxfit(X=X, theta=initheta, alpha=alpha, tmax=2, verb=0, admix=TRUE, grp=NULL, ...)
    if(verb){ cat(paste(Kseq[i],",", sep="")) }

    if(i<nK){ initheta <- tpxThetaStart(X, fit$theta, fit$omega, Kseq[i+1]) }
    else{ initheta <- fit$theta }
  }
  if(verb){ cat("done.\n") }
  return(initheta)
}

## Topic estimation and selection for a list of K values
tpxSelect <- function(X, K, bf, initheta, alpha, tol, kill, verb, admix, grp,  ...){

  if(length(K)==1 && bf==FALSE){
    if(verb){ cat(paste("Fitting the",K,"topic model.\n")) }
    return( tpxfit(X=X, theta=initheta, alpha=alpha,
                   tol=tol, verb=verb, admix=admix, grp=grp, ...) ) }

  if(is.matrix(alpha)){ stop("Matrix alpha only works for fixed K") }
  
  if(verb){ cat(paste("Fit and Bayes Factor Estimation for K =",K[1]))
            if(length(K)>1){ cat(paste(" ...", max(K))) }
            cat("\n") }

  ## dimensions
  n <- nrow(X)
  p <- ncol(X)
  nK <- length(K)
    
  BF <- NULL
  D <- NULL
  iter <- 0
  
  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))
  
  ## allocate and initialize
  best <- -Inf
  bestfit <- NULL  
  
  ## loop over topic numbers
  for(i in 1:nK){
    
    ## Solve for map omega in NEF space
    fit <- tpxfit(X=X, theta=initheta, alpha=alpha, tol=tol, verb=verb, admix=admix, grp=grp,...)
    
    BF <- c(BF, tpxML(X=X, theta=fit$theta, omega=fit$omega, alpha=fit$alpha, L=fit$L, admix=admix, grp=grp) - null)
    R <- tpxresids(X=X, theta=fit$theta, omega=fit$omega, grp=grp)
    D <- c(D, R$d)
    
    if(verb>0) cat(paste("log BF(", K[i], ") =", round(BF[i],2)))
    if(verb>1) cat(paste(" [ ", fit$iter,"steps, disp =",round(D[i],2)," ]\n")) else if(verb >0) cat("\n")
    
    if(is.nan(BF[i])){ 
      cat("NAN for Bayes factor.\n")
      return(bestfit)
      break} 
    
    if(BF[i] > best){ # check for a new "best" topic
      best <- BF[i]
      bestfit <- fit
    } else if(kill>0 && i>kill){ # break after kill consecutive drops
      if(prod(BF[i-0:(kill-1)] < BF[i-1:kill])==1) break }
    
    if(i<nK){
      if(!admix){ initheta <- tpxinit(X,2,K[i+1], alpha, 0) }
      else{ initheta <- tpxThetaStart(X, fit$theta, fit$omega, K[i+1], R=R) }
    }
  }

  names(BF) <- names(D) <- paste(K[1:length(BF)]) 
 
  return(list(theta=bestfit$theta, omega=bestfit$omega, alpha=bestfit$alpha, BF=BF, D=D, K=K[which.max(BF)])) }

                    
## ** main workhorse function.  Only Called by the above wrappers.
## topic estimation for a given number of topics (taken as ncol(theta))
tpxfit <- function(X, theta, alpha, tol=0.001, verb, admix, grp,
                   nef=TRUE, tmax=1000, wtol=10^{-4}, qnewt=1, sqp=TRUE, check=TRUE)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
  K <- ncol(theta)
  n <- nrow(X)
  p <- ncol(X)
  m <- row_sums(X)
  if(!nef) qnewt <- FALSE # can't guarentee finite map to NEF space
  if(is.null(alpha)){ alpha <- 1/(K*p) }
  if(is.matrix(alpha)){ if(nrow(alpha)!=p || ncol(alpha)!=K){ stop("bad matrix alpha dimensions") }}

  ## recycle these in tpcweights to save time
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  
  ## Initialize
  omega <- tpxOmegaStart(X=X, theta=theta)
  if(!admix){ omega <- matrix(apply(omega,2, function(w) tapply(w,grp,mean)), ncol=K) }
  iter <- 0
  dif <- tol+1
  if(verb>0){
    cat("abs(theta.diff): " )
    digits <- max(1, -floor(log(tol, base=10))) }
  
  Y <- NULL # only used for qnewt > 0 
  L <- NULL
  
  ## Iterate towards MAP
  while(abs(dif) > tol && iter < tmax){ 
    oldtheta <- theta
    
    ## sequential quadratic programming for conditional Y solution
    if(sqp && admix){ Wfit <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                start=omega, theta=theta,  verb=0, nef=nef, wtol=wtol, tmax=20) }
    else{ Wfit <- omega }

    ## joint parameter EM update
    move <- tpxEM(X=X, m=m, theta=theta, omega=Wfit, alpha=alpha, nef=nef, admix=admix, grp=grp)
    
    ## possible quasi-newton acceleration
    if(qnewt > dif){
      qn <- tpxQN(move=move, Y=Y, X=X, L=L, alpha=alpha, verb=verb, check=check, admix=admix, grp=grp)
      move <- qn$move
      Y <- qn$Y
      L <- qn$L }
    
    ## Iterate and check
    iter <- iter+1
    theta <- move$theta
    omega <- move$omega
    dif <- sum(abs(oldtheta-theta))
        
     ## LHD diff can be negative near convergence, since the Ws are not exact solutions
    if(verb>0 && (iter-1)%%ceiling(10/verb)==0){
      cat( paste( round(dif,digits), ", ", sep="") )
    }
  }

  ## final log posterior
  L <- tpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, nef=nef, admix=admix, grp=grp) 

  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }
  
  out <- list(theta=theta, omega=omega, K=K, alpha=alpha, L=L, iter=iter)
  invisible(out) }

 
## ** called from topics.R (predict) and tpx.R
## Conditional solution for topic weights given theta
tpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start) 
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="textir")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }

## ** Called only in tpx.R

## single EM update. two versions: admix and mix
tpxEM <- function(X, m, theta, omega, alpha, nef, admix, grp)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(theta)

  if(admix){ Xhat <- (X$v/tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j))*(omega[X$i,]*theta[X$j,])
             Zhat <- .C("Rzhat", n=as.integer(n), p=as.integer(p), K=as.integer(K), N=as.integer(nrow(Xhat)),
                         Xhat=as.double(Xhat), doc=as.integer(X$i-1), wrd=as.integer(X$j-1),
                        zj = as.double(rep(0,K*p)), zi = as.double(rep(0,K*n)), PACKAGE="textir")
             theta <- freq(matrix(Zhat$zj+alpha, ncol=K), byrow=FALSE)
             omega <- freq(matrix(Zhat$zi+nef/K, ncol=K)) }
  else{
    qhat <- tpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ## EM update
    theta <- freq(tcrossprod_simple_triplet_matrix( t(X), t(qhat) ) + alpha, byrow=FALSE)
    omega <- freq(matrix(apply(qhat*m,2, function(x) tapply(x,grp,sum)), ncol=K)+nef/K )  }
    
  return(list(theta=theta, omega=omega)) }

## Quasi Newton update for q>0 
tpxQN <- function(move, Y, X, L, alpha, verb, check, admix, grp)
{
  if(is.null(Y) || ncol(Y) < 3){
    return(list(Y=cbind(Y, tpxToNEF(theta=move$theta, omega=move$omega)), move=move, L=L)) }

  if(check && is.null(L)){ L <- tpxlpost(X=X, theta=move$theta, omega=move$omega,
                                         alpha=alpha, nef=TRUE, admix=admix, grp=grp) }
  
  ## update Y accounting
  Y[,1:2] <- Y[,-1]
  Y[,3] <- tpxToNEF(theta=move$theta, omega=move$omega)

  ## Check qnewt secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sVU <- sum(V*U)
  Ymove <- V*(sVU/(sum(U^2)-sVU)) 
  qnup <- tpxFromNEF(Y[,3] + Ymove,n=nrow(move$omega),
                     p=nrow(move$theta), K=ncol(move$theta))

  ## check for a likelihood improvement
  if(check){
    Lqnup <- try(tpxlpost(X=X, theta=qnup$theta, omega=qnup$omega,
                          alpha=alpha, nef=TRUE, admix=admix, grp=grp), silent=TRUE)
    if(inherits(Lqnup, "try-error")){
      if(verb>10){ cat("(QN: try error) ") }
       L <- tpxlpost(X=X, theta=move$theta, omega=move$omega,
                     alpha=alpha, nef=TRUE, admix=admix, grp=grp)
      return(list(Y=Y, move=move, L=L)) }
    if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }      
    if(Lqnup < L){
       L <- tpxlpost(X=X, theta=move$theta, omega=move$omega,
                     alpha=alpha, nef=TRUE, admix=admix, grp=grp)
       return(list(Y=Y, move=move, L=L)) }
    else{
      L <- Lqnup
      Y[,3] <- tpxToNEF(theta=qnup$theta, omega=qnup$omega)
      return( list(Y=Y, move=qnup, L=L) )
    }
  }
}
  
## unnormalized log posterior (objective function)
tpxlpost <- function(X, theta, omega, alpha, nef, admix=TRUE, grp=NULL)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  K <- ncol(theta)

  if(admix){ L <- sum( X$v*log(tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)) ) }
  else{ L <- sum(tpxMixQ(X, omega, theta, grp)$lqlhd) }
  if(is.null(nrow(alpha))){ if(alpha != 0){ L <- L + sum(alpha*log(theta))  } } # unnormalized prior
  if(nef){ L <- L + sum(log(omega))/K }
  
  return(L) }

## log marginal likelihood
tpxML <- function(X, theta, omega, alpha, L, admix, grp){
  ## get the indices
  K <- ncol(theta)
  p <- nrow(theta)
  n <- nrow(X)
  nw <- nrow(omega)

  ## return BIC for simple finite mixture model
  if(!admix){
    qhat <- tpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ML <- sum(X$v*log(row_sums(qhat[X$i,]*theta[X$j,])))
    return( ML - 0.5*( K*p + (K-1)*nw )*log(sum(X)) ) } 

  ML <- L  # assumes NEF parameterization has been used.
  ML <- ML + (K*p +(K-1)*n)*log(2*pi)/2  # (d/2)log(2pi)  
  if(is.null(nrow(alpha))){ # theta prior normalizing constant
    ML <- ML + K*( lgamma(p*(alpha+1)) - p*lgamma(alpha+1) )  }
  else{ ML <- ML + sum(lgamma(col_sums(alpha+1)) - col_sums(lgamma(alpha+1))) } # matrix version
  ML <- ML - (nw+1)*lfactorial(K) # omega(phi) prior normalizing constant + multiplier for # of modes
  
  ## block-diagonal approx to determinant of the negative log hessian matrix
  q <- tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
  D <- tpxHnegDet(X=X, q=q, theta=theta, omega=omega, alpha=alpha)
  ML <- ML - 0.5*sum( D  )   # -1/2 log |-H|

  return(ML) }

## find residuals for X$v
tpxresids <- function(X, theta, omega, grp=NULL)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }

  m <- row_sums(X)

  if(is.null(grp)){ xhat <- tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)*m[X$i] }
  else{
    qhat <- tpxMixQ(X=X, omega=omega, theta=theta, grp=grp, qhat=TRUE)$qhat
    xhat <- row_sums(qhat[X$i,]*theta[X$j,])*m[X$i] }
  
  s <- sqrt(xhat*(m[X$i]-xhat)/m[X$i])
  e <- (X$v-xhat) 

  r <- e/s
  d <- mean(r^2)
  
  return( list(s=s, e=e, r=r, d=d) ) }

## fast initialization functions for theta (after increasing K) and omega (given theta)
tpxThetaStart <- function(X, theta, omega, K, R=NULL)
  {
    if(is.null(R)){ R <- tpxresids(X, theta=theta, omega=omega) }
    X$v <- R$e*(R$r>3) + 1/ncol(X)
    Kpast <- ncol(theta)
    Kdiff <- K-Kpast
    if(Kpast != ncol(omega) || Kpast >= K){ stop("bad K in tpxThetaStart") }
    initheta <- freq(Kpast*theta+rowMeans(theta), byrow=FALSE)
    n <- nrow(X)
    ki <- matrix(1:(n-n%%Kdiff), ncol=Kdiff)
    for(i in 1:Kdiff){ initheta <- cbind(initheta, (col_sums(X[ki[,i],])+1/ncol(X))/(sum(X[ki[,i],])+1)) }
    return( initheta )
  }

tpxOmegaStart <- function(X, theta)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
    if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
    omega[omega <= 0] <- .5
    return( freq(omega, byrow=TRUE) )
  }


## fast computation of sparse P(X) for X>0
tpxQ <- function(theta, omega, doc, wrd){

  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in tpxQ") }
  
  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="textir" )

  return( out$q ) }

## model and component likelihoods for mixture model
tpxMixQ <- function(X, omega, theta, grp, qhat=FALSE){
  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="textir")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }

## negative log hessian block diagonal matrix for theta & omega
tpxHnegDet <- function(X, q, theta, omega, alpha){
  K <- ncol(theta)
  n <- nrow(omega)

  ## sparse Xij/Qij^2 
  Xq <- X
  Xq$v <- Xq$v/q^2
  
  ## negative 2nd derivitive matrices for theta
  HT <- tcrossprod_simple_triplet_matrix(t(Xq), apply(omega, 1, function(v) v%o%v ) )
  HT[,K*(0:(K-1))+1:K] <- HT[,K*(0:(K-1))+1:K] + alpha/theta^2 # will break for alpha<=1
  DT <- apply(HT, 1, tpxlogdet)

  ## ditto for omega
  HW <- matrix(.C("RnegHW",
                  n = as.integer(nrow(omega)),
                  p = as.integer(nrow(theta)),
                  K = as.integer(K-1),
                  omeg = as.double(omega[,-1]),
                  thet = as.double(theta[,-1]),
                  doc = as.integer(X$i-1),
                  wrd = as.integer(X$j-1),
                  cnt = as.double(X$v),
                  q = as.double(q),
                  N = as.integer(length(q)),
                  H = double(n*(K-1)^2),
                  PACKAGE="textir")$H,
               nrow=(K-1)^2, ncol=n)
  DW <- apply(HW, 2, tpxlogdet)
           
  return( c(DT,DW) ) }

## functions to move theta/omega to and from NEF.  
tpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="textir")$Y)
}

## 'From' NEF representation back to probabilities
tpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="textir")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}

## utility log determinant function for speed/stabilty
tpxlogdet <- function(v){
    v <- matrix(v, ncol=sqrt(length(v)))
    
    if( sum(zeros <- colSums(v)==0)!=0 ){
      cat("warning: boundary values in laplace approx\n")
      v <- v[-zeros,-zeros] }
   
    return(determinant(v, logarithm=TRUE)$modulus) }

