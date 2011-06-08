##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
topics <- function(counts, K, shape=NULL, initopics=5, tol=0.01, admix=TRUE, grp=NULL,
                   bf=FALSE, kill=2, ord=TRUE, verb=1, ...)
  ## tpxfit optional arguments are nef=TRUE, tmax=1000, wtol=10^(-4), qnewt=1, sqp=1, check=1
{
  ## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }

  if(prod(row_sums(counts)>0) !=1){ stop("You've included empty documents.") }
  
  X <- as.simple_triplet_matrix(counts)
  p <- ncol(X) 

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }
                
  ## check the list of candidate K values
  if(prod(K>1)!=1){ stop(cat("use K values > 1\n")) }
  K <- sort(K)

  ## check grp if mixture
  if(!admix){
    if(is.null(grp) || length(grp)!=nrow(X)){  grp <- rep(1,nrow(X)) }
    else{ grp <- factor(grp) }
  }
  
  ## initialize
  initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1], shape, verb)
  
  ## either search for marginal MAP K and return bayes factors, or just fit
  tpx <- tpxSelect(X, K, bf, initopics, alpha=shape, tol, kill, verb, admix, grp, ...)
  K <- tpx$K
  
  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(counts)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(admix){ dimnames(omega)[[1]] <- dimnames(counts)[[1]] }
  ## residuals (only for positive count entries)
  rpos <- tpxresids(X, theta=tpx$theta, omega=tpx$omega, grp=grp)$r
  ## topic object
  out <- list(K=K, theta=theta, omega=omega,
              BF=tpx$BF, dispersion=tpx$D, residuals=rpos,
              X=X, admix=admix, grp=grp)
  class(out) <- "topics"
  invisible(out) }

## S3 method predict function
predict.topics <- function(object, newcounts, grp=NULL, ...)
  ## tpxweights optional arguments and defauls are verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000
{
  if(is.vector(newcounts)){ newcounts <- matrix(newcounts, nrow=1) }
  if(class(newcounts)[1] == "TermDocumentMatrix"){ newcounts <- t(newcounts) }
  X <- as.simple_triplet_matrix(newcounts)
  
  if(class(object)!="topics"){ stop("object class must be `topics'.") }

  theta <- object$theta
  if(nrow(theta) != ncol(X)){ stop("Dimension mismatch: nrow(theta) != ncol(X)") }

  if(!object$admix){
    if(is.null(grp)){ grp <- rep(1, nrow(newcounts)) }
    Q <- matrix(tpxMixQ(X, omega=object$omega, theta=theta, grp=grp)$lQ, ncol=ncol(theta))
    return( (1:ncol(theta))[apply(Q,1,which.max)] ) }
  
  start <- tpxOmegaStart(X=X, theta=theta)
  
  ## re-order xvec in doc-blocks, and build indices
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  xv <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1

  return(tpxweights(n=nrow(X), p=ncol(X), xv=xv, wrd=wrd, doc=doc, start=start, theta=theta, ...))  }

## S3 method summary function
summary.topics <- function(object, nwrd=5, tpk=NULL, verb=TRUE, ...){

    
  K <- object$K
  if(is.null(tpk)){ tpk <- 1:K }
  else if(prod(tpk %in% 1:K)!=1){ stop("requested tpk's are not in 1:K") }
  
  if(verb){ cat(paste("\nTopic % Usage: \n\n")) 
            print(round( (col_means(object$omega)*100)[tpk], 1) ) }

  if(nwrd>0){ if(verb){ cat(paste("\nTop", nwrd, "phrases by topic-over-null term lift:\n\n")) }
              Q0 <- col_sums(object$X)/sum(object$X)
              topwords <- c()
              for(k in tpk){
                odds <- (log(object$theta[,k]) - log(Q0))[Q0!=0]
                ko <- order(odds, decreasing=TRUE)
                topk <- dimnames(object$theta)[[1]][Q0!=0][ko[1:nwrd]]
                topwords <- cbind(topwords, topk)
                if(verb){ cat(paste("[",k, "] '"))
                          cat(topk, sep="', '")
                          cat("'\n\n") }
              }
              topwords <- as.matrix(topwords)
              dimnames(topwords)[[2]] <- paste(tpk)
            }
  else{ topwords <- NULL
        if(verb){ cat("\n\n") } }
              
  if(!is.null(object$BF) && !is.null(object$dispersion) && verb)
    {
      cat("Log Bayes factor and estimated dispersion, by number of topics:\n\n")
      print(round(rbind(lBF=object$BF,r2=object$dispersion),2))
      cat(paste("\nSelected the K =",object$K,"topic model\n\n"))
    }

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
    
    hist(x$resid, col=col, border=grey(.9),
         xlab=xlab,
         main="", cex.lab=cex.lgdt, font.lab=3)
    
    return(invisible())
  }
  
  n <- nrow(x$omega)
  if(n==1){
    if(is.null(main)){ main="" }
    if(is.null(xlab)){ xlab="topic" }
    if(is.null(ylab)){ ylab="weight" }
    if(is.null(col)){ col = 8 }
    return(plot(c(x$omega), type="h", lwd=10, col=col, main=main, ylab=ylab, xlab=xlab)) }

  if(is.null(tpk)){ tpk <- 1:x$K }
  if(is.null(lgd.K)){ lgd.K <- max(.1*length(tpk),.75) }
  
  if(is.null(group) || !x$admix){
    ltop = .65*n
    lbot = .35*n
    
    if(is.null(col)){ col<-1 }
    tpcl <- c(0, TOPICOLS[1:6,col[1]%%5])
    W <- x$omega[,tpk]
    if(x$admix){
      brks <- seq(0,1,length=8)
      tplg <- c("w=0","w=1") }
    else{
      brks <- seq(min(W),max(W),length=8)
      tplg <- c(round(min(W),3),round(max(W),3))
    }
  } else{ # bunch of extra commands to get shading for two groups
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
  if(is.null(ylab)){
    if(x$admix){ ylab="Document" }
    else{ ylab="group" } }
  if(is.null(xlab)){ xlab="Topic" }
  if(is.null(main)){ main="Topic-Loading Weights" }

  ## finally plot
  image(y=1:n, x=1:length(tpk), t(W), ylab=ylab, xlab=xlab,
        main=main, col=tpcl, font.lab=3, xaxt="n", yaxt="n", breaks=brks, ...)
  axis(side=1, at=1:length(tpk), labels=tpk, tick=FALSE, line=-.5)
  if(x$admix){ axis(2) }
  else{ axis(side=2, at=1:n, labels=levels(as.factor(x$grp))) }
  points(rep(xlg,length(tpcl)), seq(lbot,ltop,length=length(tpcl)), col=tpcl, pch=15, cex=3*cex.lgdc)
  text(rep(xlg,2), y=c(lbot-.08*n, ltop+.08*n), tplg, cex=cex.lgdt)
  if(!is.null(labels)){ text(rep(xlg,2), y=c(lbot-.14*n, ltop+.14*n), labels, font=3, cex=cex.lgdt) }
  par(mar=old.mar)
}

  
