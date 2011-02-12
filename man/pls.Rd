\name{pls}
\alias{pls}
\title{ Partial Least Squares }
\description{ A simple partial least squares procedure. }
\usage{ pls(F, y, K=1, scale=TRUE, verb=TRUE) }
\arguments{
  \item{F}{ The covariate matrix, in either \code{simple_triplet_matrix} or \code{matrix} format. }
   \item{y}{ The response vector.  }
   \item{K}{ The number of desired PLS directions. }
   \item{scale}{ An indicator for whether to standardize \code{F}; usually a good idea. 
   		 If \code{scale=TRUE}, F will be scaled to have variance-one columns.   }
   \item{verb}{ Whether or not to print a small progress script. }
}
\value{ A pls object list with the following entries
	\item{y}{The response vector. }
	\item{F}{The covariate matrix. If \code{scale=TRUE}, scaled to have variance-one columns. }
	\item{z}{The pls directions \code{F\%*\%phi}. }
	\item{v}{Response factors. }
	\item{yhat}{\code{K} columns of fitted values for each number of directions. }
	\item{fwdmod}{ The \code{lm} object from forward regression \code{lm(as.numeric(y)~z)}. }
	\item{scale}{ Indicator for whether input \code{F} was scaled. }
}
\details{ Fits the Partial Least Squares algorithm described in Taddy (2011; Section 3.1). 
	In particular, 	we obtain loadings \code{phi[,k]} as the correlation between
	\code{F} and factors \code{v[,k]}, where \code{v[,1]} is initialized      
	at \code{normalize(as.numeric(y))} and subsequent factors are orthogonal 
	to the k'th pls direction, \code{z[,k]=F\%*\%phi[,k]}. }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}

Wold, H. (1975), \emph{Soft modeling by latent variables: The nonlinear iterative partial least squares approach}. 
In Perspectives in Probability and Statistics, Papers in Honour of M.S. Bartlett.

}
\author{ Matt Taddy \email{taddy@chicagobooth.edu} }
\seealso{
  plot.pls, normalize, freq, corr, we8there, congress109
}
\examples{
data(congress109)
fit <- pls(freq(congress109Counts), congress109Ideology$repshare, K=4)
plot(fit, pch=21, bg=c(4,3,2)[congress109Ideology$party])

data(we8there)
fit <- pls(freq(we8thereCounts),as.factor(we8thereRatings$Overall))
plot(fit, col=c(2,2,2,3,3))
}