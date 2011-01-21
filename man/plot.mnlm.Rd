\name{plot.mnlm}
\alias{plot.mnlm}
\alias{coefplot}
\title{
Multinomial logistic regression Plots
}
\description{ Plot function for \code{mnlm} objects, the output of multinomial logistic regression. }
\usage{
\method{plot}{mnlm}(x, covar=NULL, y=NULL, cat=NULL, delta=NULL, ...)
}
\arguments{
  \item{x}{An output object from the \code{mnlm} function.}
\item{covar}{The covariate direction to plot.  Defaults to 2 (first past
		 the intercept).}
\item{y}{Optional argument for the response to be plotted against (if, e.g., you 
wish to plot against unnormalized response).}
\item{cat}{An optional category in \code{1:(ncol(fit$X))}.  If specified, plot the 
	     conditional bound and likelihood functions for that category around
	     the converged solution.  }
\item{delta}{Optional argument for the delta used by the bounding function (and plot \code{xlim}) when \code{cat} is specified. 
		      Defaults to \code{fit$delta[cat,covar]}. }
 \item{...}{Additional \code{plot} arguments}
}

\details{ Can be used to reproduce all of the plots (except wordles;
see \url{http://wordle.net}) from Taddy (2011). If \code{cat=NULL}, just
plots the fitted reductions \code{freq(counts)\%*\%coef[,covar]} for
inverse regression based on this \code{mnlm} fit.  If \code{cat} is
specified, plots the objective near solution for that multinomial
category coefficient. }

\value{ A fabulous plot.}
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\note{\code{coefplot} is an internal function for when \code{cat} is specified.}
\seealso{
 mnlm, congress109, we8there
}
