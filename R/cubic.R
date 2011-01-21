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
