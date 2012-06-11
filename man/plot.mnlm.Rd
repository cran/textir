\name{plot.mnlm}
\alias{plot.mnlm}
\title{
Multinomial logistic regression Plots
}
\description{ Plot function for \code{mnlm} objects, the output of multinomial logistic regression. }
\usage{
\method{plot}{mnlm}(x, type=c("response","reduction","roc"), covar=NULL, v=NULL, xlab=NULL, ylab=NULL, col=NULL, ...)
}
\arguments{
  \item{x}{An output object from the \code{mnlm} function.}
  \item{type}{ Under \code{"response"}, plot the fitted count expectations
  	       against observed non-zero counts.
	       Under \code{"reduction"},  plot the sufficient reduction scores 
	       \code{freq(counts)\%*\%loadings} from
	       inverse regression based on this \code{mnlm} fit.  
	       Under \code{"roc"},  plot the receiver operating characteristic 
	       for classification based on the fitted model
	       [Note: the roc plot only applies for data with \code{max(counts)==1} and can be slow if \code{ncol(counts)} is very large].  }
\item{covar}{For \code{type="reduction"}. The covariate direction to plot.  Defaults to 1.}
\item{v}{For \code{type="reduction"}. Optional argument for the fitted reduction to be plotted against (if, e.g., you 
wish to plot against unnormalized response).}
\item{xlab}{ The x-axis label; will be automatically set if NULL. }
\item{ylab}{ The y-axis label; will be automatically set if NULL. For binary data, this becomes the legend title. }
\item{col}{ The color(s). Usage changes depending on plot type (for roc, it must be a \code{ncol(predict(x,x$covars))} length vector). }
 \item{...}{Additional \code{plot} arguments}
}

\value{ A fabulous plot.}
\references{
   
Taddy (2012), \emph{Multinomial Inverse Regression for Text Analysis}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
 mnlm, congress109, we8there
}
