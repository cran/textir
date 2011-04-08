##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
topics <- function(counts, K, alpha=NULL, initheta=NULL, tol=0.1, 
                   bf=FALSE, kill=2, ord=TRUE, verb=1, ...)
  ## tpxfit optional arguments and defauls are nef=TRUE, tmax=1000, wtol=10^(-5), qnewt=1, sqp=TRUE,             
{
  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }

  if(prod(row_sums(counts)>0) !=1){ stop("You've included empty documents.") }
  
  X <- as.simple_triplet_matrix(counts)
  p <- ncol(X) 

  ## check the prior parameters for theta
  if(prod(alpha>1) != 1){ stop("use alpha > 1\n") }
                            
  ## check the list of candidate K values
  if(prod(K>1)!=1){ stop(cat("use K values > 1\n")) }
  K <- sort(K)

  ## inital theta
  if(!is.null(initheta)){
    if(ncol(initheta)!=K[1]){ stop("mis-match between initheta and K.") } 
    initheta <- freq(initheta + 0.01/p, byrow=FALSE)  }
  else{
    if(verb){ cat("Finding initial topics for K = ") }
    initheta <- tpxBF(X=X,
                      K=unique(ceiling(seq(2,K[1],length=10))),
                      initheta=tpxThetaInit(X, 2),
                      alpha=alpha, tol=1, tmax=2, kill=0, verb=verb, init=TRUE)$theta
  if(verb){ cat("done.\n") }
  }

  ## either search for marginal MAP K and return bayes factors, or just fit
  BF <- NULL
  tpx <- NULL
  iter <- 0
  if(bf || length(K)>1){
    if(verb){ cat(paste("Fit and Bayes Factor Estimation for K =",K[1]))
              if(length(K)>1){ cat(paste(" ...", max(K))) }
              cat("\n") }
    tpx <- tpxBF(X=X, K=K, initheta=initheta, alpha=alpha, tol=tol, kill=kill, verb=verb, ...)
    initheta <- tpx$theta
    BF <- tpx$BF
    K <- tpx$K
  }
  else{
    if(verb){ cat(paste("Fitting the",K,"topic model.\n")) }
    tpx <- tpxfit(X=X, theta=initheta, alpha=alpha, tol=tol, verb=verb, ...)
  }
  
  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(counts)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=dimnames(counts)[[1]], topic=paste(1:K)) )
  ## residuals (only for positive count entries)
  resids <- X
  resids$v <- tpxresids(X, theta=theta, omega=omega)$r
  ## topic object
  out <- list(K=K, theta=theta, omega=omega, BF=BF, residuals=resids, X=X)
  class(out) <- "topics"
  invisible(out) }

## S3 method predict function
predict.topics <- function(object, newcounts, ...)
  ## tpxweights optional arguments and defauls are verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000
{
  if(is.vector(newcounts)){ newcounts <- matrix(newcounts, nrow=1) }
  if(class(newcounts)[1] == "TermDocumentMatrix"){ newcounts <- t(newcounts) }
  X <- as.simple_triplet_matrix(newcounts)
  
  if(class(object)!="topics"){ stop("object class must be `topics'.") }

  theta <- object$theta
  if(nrow(theta) != ncol(X)){ stop("Dimension mismatch: nrow(theta) != ncol(X)") }
  start <- tpxOmegaInit(X=X, theta=theta)
  
  ## re-order xvec in doc-blocks, and build indices
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  xv <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1

  return(tpxweights(n=nrow(X), p=ncol(X), xv=xv, wrd=wrd, doc=doc, start=start, theta=theta, ...))  }

## S3 method summary function
summary.topics <- function(object, nwrd=5, tpk=NULL, ...){

  K <- object$K
  if(is.null(tpk)){ tpk <- 1:K }
  else if(prod(tpk %in% 1:K)!=1){ stop("requested tpk's are not in 1:K") }
  
  cat(paste("\nTopic % Usage: \n\n"))
  print(round( (col_means(object$omega)*100)[tpk], 1) )

  cat(paste("\nTop", nwrd, "phrases by topic-over-null log odds:\n\n"))
  Q0 <- col_sums(object$X)/sum(object$X)
  topwords <- c()
  for(k in tpk){
    odds <- (log(object$theta[,k]) - log(Q0))[Q0!=0]
    ko <- order(odds, decreasing=TRUE)
    topk <- dimnames(object$theta)[[1]][Q0!=0][ko[1:nwrd]]
    topwords <- cbind(topwords, topk)
    cat(paste("[",k, "] '"))
    cat(topk, sep="', '")
    cat("'\n\n") }

  topwords <- as.matrix(topwords)
  dimnames(topwords)[[2]] <- paste(tpk)
  invisible(topwords)
}

## Colors for topic plotting
TOPICOLS <- matrix(nrow=6,
                   c(grey(c(.9,.8,.7,.55,.35,0)), #GREY
                     "#FFCCCC", "#FA8072", "#FF5333", "#EE0000", "#CC1100", "#800000", #RED
                     "#BDFCC9", "#98FB98", "#49E20E", "#009900", "#006400", "#004F00", #GREEN	
                     "#BBFFFF", "#67E6EC", "#00C5CD", "#0198E1", "#0147FA",  "#000080"))  #BLUE

              	

plot.topics <- function(x, type=c("weight","resid"), group=NULL, labels=NULL, 
                        col=NULL, xlab=NULL, ylab=NULL, main=NULL, tpk=NULL,
                        lgd.K=NULL, cex.lgdc = 1, cex.lgdt = 1, cex.rmar= 1, ...){

  if(type[1]=="resid"){

    if(is.null(col[1])){col <- 8}
    if(is.null(xlab)){ xlab="Standardized residuals for non-zero counts" }
    if(is.null(main)){ main="" }
    
    hist(x$resid$v, col=8, border=grey(.9),
         xlab=xlab,
         main="", cex.lab=cex.lgdt, font.lab=3)
    
    return(invisible())
  }

  if(is.null(tpk)){ tpk <- 1:x$K }
  
  n <- nrow(x$omega)
  if(is.null(lgd.K)){ lgd.K <- max(.1*length(tpk),.75) }
  
  if(is.null(group)){
    ltop = .65*n
    lbot = .35*n
    
    if(is.null(col)){ col<-1 }
    tpcl <- c(0, TOPICOLS[1:6,col%%5])
    brks <- seq(0,1,length=8)
    W <- x$omega[,tpk]
    tplg <- c("w=0","w=1")
  }
  else{ # bunch of extra commands to get shading for two groups
    group <- as.factor(group)
    ltop = .85*n
    lbot = .15*n
    
    if(length(group)!=n){ stop("Your group membership length doesn't match omega.") }
    if(nlevels(group)!=2){ stop("Sorry, grouping must be a two-level factor") }
    if(is.null(col) || length(col)<2){ col <- 1:2 }
    
    tpcl <- c(TOPICOLS[6:1,col[1]%%5],"#FFFFFF",TOPICOLS[,col[2]%%5])
    brks <- c(seq(-1, -0.1,length=7),seq(0.1,1,length=7))
    W <- x$omega[,tpk]*(c(-1,1)[as.logical(group)+1])
    if(is.null(labels)){labels <- c("F","T")}
    tplg=rep("w=1",2)
  }

  ## plot parameters
  xlg <- length(tpk)+lgd.K
  old.mar <- par()$mar
  par(xpd=TRUE, mar=c(5.1,4.1,2.1,5*cex.rmar))
  if(is.null(ylab)){ ylab="Document" }
  if(is.null(xlab)){ xlab="Topic" }
  if(is.null(main)){ main="Topic-Loading Weights" }

  ## finally plot
  image(y=1:n, x=1:length(tpk), t(W), ylab=ylab, xlab=xlab, main=main, col=tpcl, font.lab=3, xaxt="n",
        breaks=brks, ...)
  axis(side=1, at=1:length(tpk), labels=tpk, tick=FALSE, line=-.5)
  points(rep(xlg,length(tpcl)), seq(lbot,ltop,length=length(tpcl)), col=tpcl, pch=15, cex=3*cex.lgdc)
  text(rep(xlg,2), y=c(lbot-.08*n, ltop+.08*n), tplg, cex=cex.lgdt)
  if(!is.null(labels)){ text(rep(xlg,2), y=c(lbot-.14*n, ltop+.14*n), labels, font=3, cex=cex.lgdt) }
  par(mar=old.mar)
}

  
