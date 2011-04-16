#######  Undocumented "tpx" utility functions #########

## ** Only referenced from topics.R

## Bayes Factor estimation for a list of K values
tpxBF <- function(X, K, initheta, alpha, tol, kill, verb, init=FALSE, ...){

  ## dimensions
  n <- nrow(X)
  p <- ncol(X)
  nK <- length(K)
  
  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))

  ## allocate and initialize
  BF <- c()
  best <- -Inf
  bestfit <- NULL
  rX <- X
  if(init) verb <- 0
  
  ## loop over topic numbers
  for(i in 1:nK){

    ## Solve for map omega in NEF space
    fit <- tpxfit(X=X, theta=initheta, alpha=alpha, tol=tol, verb=verb, ml=TRUE, ...)
    if(!init){ BF <- c(BF, tpxML(X=X, theta=fit$theta, omega=fit$omega, alpha=fit$alpha) - null)
               if(verb>0) cat(paste("log BF(", K[i], ") =", round(BF[i],2)))
               if(verb>1) cat(paste(" [", fit$iter,"steps ]\n")) else if(verb >0) cat("\n")
    
               if(is.nan(BF[i])){ 
                 cat("NAN for Bayes factor.\n")
                 return(bestfit)
                 break} 
    
               if(BF[i] > best){ # check for a new "best" topic
                 best <- BF[i]
                 bestfit <- fit
               } else if(kill>0 && i>kill){ # break after kill consecutive drops
                 if(prod(BF[i-0:(kill-1)] < BF[i-1:kill])==1) break; }
             }
    else{ cat(paste(K[i],",", sep="")) }
    
    if(i<nK){ # new topic from residuals
      R <- tpxresids(X, theta=fit$theta, omega=fit$omega)
      rX$v <- R$e*(R$r>3) + 1/p
      initheta <- cbind(freq(K[i]*fit$theta+rowMeans(fit$theta), byrow=FALSE), tpxThetaInit(rX, K[i+1]-K[i])) } 
  }
  if(!init){
    names(BF) <- paste(K[1:length(BF)]) 
    k<-which.max(BF) }
  else{
    k<-nK
    bestfit <- fit
  }
  return(list(theta=bestfit$theta, omega=bestfit$omega, alpha=bestfit$alpha, BF=BF, K=K[k])) }

## ** called from topics.R and tpx.R

## topic estimation for a given number of topics (taken as ncol(theta))
tpxfit <- function(X, theta, alpha, tol, verb, nef=TRUE, tmax=1000, wtol=10^{-5}, qnewt=1, sqp=TRUE, ml=FALSE)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
  K <- ncol(theta)
  n <- nrow(X)
  p <- ncol(X)
  nef = nef||ml
  if(!nef) qnewt <- 0 # can't guarentee finite map to NEF space
  if(is.null(alpha)){ alpha <- 1+1/(K*p) }

  ## recycle these in tpcweights to save time
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  
  ## Initialize
  omega <- tpxOmegaInit(X=X, theta=theta)
  L <-  tpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, nef=nef)
  if(qnewt>0)
    { yvec <- tpxToNEF(theta=theta, omega=omega)
      Y <- matrix(nrow=length(yvec), ncol=qnewt+2)
      Y[,qnewt+2] <- yvec  }
  iter <- 0
  dif <- tol+1
  if(verb==1) cat("LHD diff: " )
  
  ## Iterate towards MAP
  while(abs(dif) > tol && iter < tmax){ 

    ## sequential quadratic programming for conditional Y solution
    if(sqp){ WSQP <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                start=omega, theta=theta,  verb=verb-3, nef=nef, wtol=wtol, tmax=20) }
    else{ WSQP <- omega }

    ## joint parameter EM update
    move <- tpxEM(X=X, theta=theta, omega=WSQP, alpha=alpha, nef=nef)
    
    ## possible quasi-newton acceleration
    if(qnewt>0){
      Y[,1:(qnewt+1)] <- Y[,-1]
      Y[,qnewt+2] <- tpxToNEF(theta=move$theta, omega=move$omega)
      ## Check qnewt secant conditions and solve F(x) - x = 0.
      if(iter>qnewt){
        if(qnewt > 1){ Ymove <- tpxQN(qnewt, move, Y) }
        else{
          V <- Y[,3]-Y[,2]
          U <- Y[,2]-Y[,1]
          sVU <- sum(V*U)
          Ymove <- V*(sVU/(sum(U^2)-sVU))
        }
        qnup <- tpxFromNEF(Y[,qnewt+2] + Ymove,n=n, p=p, K=K) 
        qnup$L <- try(tpxlpost(X=X, theta=qnup$theta, omega=qnup$omega,
                               alpha=alpha, nef=TRUE), silent=TRUE)
        if(!inherits(qnup$L, "try-error"))
          if(qnup$L > move$L)
            { if(verb-2 > 0) cat(paste(round(qnup$L - move$L,2),"QN diff at "))
              move <- qnup 
              Y[,qnewt+2] <- tpxToNEF(theta=move$theta, omega=move$omega)  } }
    }#QN
    
    ## Iterate and check
    theta <- move$theta
    omega <- move$omega
    dif <- move$L-L
    L <- move$L
    iter <- iter+1

     ## diff can be negative near convergence, since the Ws are not exact solutions
    if(verb>1){ cat( paste("step", iter, ": L =", L, "diff =", dif, "\n") ) }
    else if(verb==1 && (iter-1)%%10==0){
      cat( paste( round(dif,1), ", ", sep="") )  }
  }

      
  if(verb==1){  cat("done\n")  }        
  
  out <- list(theta=theta, omega=omega, K=K, alpha=alpha, iter=iter, tol=tol, wtol=wtol)
  invisible(out) }

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

## single EM update 
tpxEM <- function(X, theta, omega, alpha, nef)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(theta)
  
  Xhat <- (X$v/tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j))*(omega[X$i,]*theta[X$j,])

  Zhat <- .C("Rzhat", n=as.integer(n), p=as.integer(p), K=as.integer(K), N=as.integer(nrow(Xhat)),
                         Xhat=as.double(Xhat), doc=as.integer(X$i-1), wrd=as.integer(X$j-1),
                         zj = as.double(rep(0,K*p)), zi = as.double(rep(0,K*n)), PACKAGE="textir") 

  theta <- freq(matrix(Zhat$zj+alpha-1, ncol=K), byrow=FALSE)
  omega <- freq(matrix(Zhat$zi+nef/K, ncol=K))
  
  L <- tpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, nef=nef)
  
  return(list(theta=theta, omega=omega, L=L)) }

## Quasi Newton update for q>1 
tpxQN <- function(qnewt, move, Y)
{
  U <- as.matrix(Y[,1:qnewt + 1]-Y[,1:qnewt])
  V <- as.matrix(Y[,1:qnewt + 2]-Y[,1:qnewt + 1])
  M <- try(solve(t(U)%*%U - t(U)%*%V), silent=TRUE)
  if(inherits(M, "try-error")){ return(0) }
  VM <- V%*%M
  UF <- t(U)%*%V[,qnewt]
  return(VM%*%UF)
}
  
## unnormalized log posterior (objective function)
tpxlpost <- function(X, theta, omega, alpha, nef)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  L <- sum( X$v*log(tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)) ) # unnormalized likelihood
  if(alpha != 1){ L <- L + (alpha-1)*sum(log(theta))  } # unnormalized prior
  if(nef){ L <- L + sum(log(omega))/ncol(theta) }
  
  return(L) }

## log marginal likelihood
tpxML <- function(X, theta, omega, alpha){
  ## get the indices
  K <- ncol(theta)
  p <- nrow(theta)
  n <- nrow(X)

  ## probabilities for X>0
  q <- tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)

  ## Build the marginal likelihood, starting with the un-normalized posterior
  ML <- sum( X$v*log(q) ) # unnormalized likelihood

  if(alpha != 1){ ML <- ML + (alpha-1)*sum(log(theta))  } # unnormalized prior
  ML <- ML + K*( lgamma(p*alpha) - p*lgamma(alpha) ) # theta prior normalizing constant  

  ML <- ML + sum(log(omega))/K # necessary here in NEF parameterization
  ML <- ML + (n+1)*lfactorial(K) # omega(phi) prior normalizing constant + multiplier for # of modes

  ML <- ML + (K*p +(K-1)*n)*log(2*pi)/2  # (d/2)log(2pi)  

  # block-diagonal approx to determinant of the negative log hessian matrix
  D <- tpxHDet(X=X, q=q, theta=theta, omega=omega, alpha=alpha)

  ## finish and return
  ML <- ML - 0.5*sum( D  )  # -1/2 log |-H|

  return(ML) }

## find residuals for X$v
tpxresids <- function(X, theta, omega)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }

  m <- row_sums(X)
  xhat <- tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)*m[X$i]
  s <- sqrt(xhat*(m[X$i]-xhat)/m[X$i])
  e <- (X$v-xhat)
  r <- e/s
  
  return( list(s=s, e=e, r=r) ) }

## random start for theta by allocating entire docs to topics
tpxThetaInit <- function(X, K)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    theta <- c()
    n <- nrow(X)
    ki <- matrix(1:(n-n%%K), ncol=K)
    for(i in 1:K){ theta <- cbind(theta, (col_sums(X[ki[,i],])+1/ncol(X))/(sum(X[ki[,i],])+1)) }
    return( theta )
  }

tpxOmegaInit <- function(X, theta)
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

  
## negative log hessian block diagonal matrix for theta & omega
tpxHDet <- function(X, q, theta, omega, alpha){
  K <- ncol(theta)
  n <- nrow(omega)

  ## sparse Xij/Qij^2 
  Xq <- X
  Xq$v <- Xq$v/q^2
  
  ## negative 2nd derivitive matrices for theta
  HT <- tcrossprod_simple_triplet_matrix(t(Xq), apply(omega, 1, function(v) v%o%v ) )
  if(alpha > 1){ HT[,K*(0:(K-1))+1:K] <- HT[,K*(0:(K-1))+1:K] + (alpha-1)/theta^2 }
  DT <- apply(HT, 1, tpxlogdet)

  ## ditto for omega
  HW <- matrix(.C("RcalcHW",
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

## 'From' back to probabilities
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

