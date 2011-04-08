\name{plot.mnlm}
\alias{plot.mnlm}
\alias{binaryplot}
\title{
Multinomial logistic regression Plots
}
\description{ Plot function for \code{mnlm} objects, the output of multinomial logistic regression. }
\usage{
\method{plot}{mnlm}(x, type=c("response","reduction"), covar=NULL, v=NULL, xlab=NULL, ylab=NULL, col=NULL, ...)
}
\arguments{
  \item{x}{An output object from the \code{mnlm} function.}
  \item{type}{ Under \code{"response"}, plot the fitted count expectations
  	       against observed non-zero counts.
	       Under \code{"reduction"}, we plot the fitted reduction \eqn{F\phi}.  }
\item{covar}{For \code{type="reduction"}. The covariate direction to plot.  Defaults to 1.}
\item{v}{For \code{type="reduction"}. Optional argument for the fitted reduction to be plotted against (if, e.g., you 
wish to plot against unnormalized response).}
\item{xlab}{ The x-axis label; will be automatically set if NULL. }
\item{ylab}{ The y-axis label; will be automatically set if NULL. For binary data, this becomes the legend title. }
\item{col}{ The color(s). Usage changes depending on plot type. }
 \item{...}{Additional \code{plot} arguments}
}

\details{ Option \code{type="reduction"} provides a plot of fitted and observed counts. 
	  Option \code{type="reduction"}, 
	  plots the fitted reductions \code{freq(counts)\%*\%loadings} from
	  inverse regression based on this \code{mnlm} fit.  }

\value{ A fabulous plot.}
\note{ \code{binaryplot} is an internal function for binary choice data plots. }
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
 mnlm, congress109, we8there
}
