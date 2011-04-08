\name{pls}
\alias{pls}
\title{ Partial Least Squares }
\description{ A simple partial least squares procedure. }
\usage{ pls(X, y, K=1, scale=TRUE, verb=TRUE) }
\arguments{
  \item{X}{ The covariate matrix, in either \code{simple_triplet_matrix} or \code{matrix} format. }
   \item{y}{ The response vector.  }
   \item{K}{ The number of desired PLS directions. }
   \item{scale}{ An indicator for whether to standardize \code{X}; usually a good idea. 
   		 If \code{scale=TRUE}, X will be scaled to have variance-one columns.   }
   \item{verb}{ Whether or not to print a small progress script. }
}
\value{ A pls object list with the following entries
	\item{y}{The response vector. }
	\item{X}{The covariate matrix. If \code{scale=TRUE}, scaled to have variance-one columns. }
	\item{directions}{The pls directions \code{X\%*\%loadings}. }
	\item{factors}{Response factors. }
	\item{phi}{The pls loadings. }
	\item{fitted}{\code{K} columns of fitted \code{y} values for each number of directions. }
	\item{fwdmod}{ The \code{lm} object from forward regression \code{lm(as.numeric(y)~directions)}. }
	\item{scale}{ If \code{scale=TRUE} on input, the standard deviations used to scale \code{X}.  }
}
\details{ Fits the Partial Least Squares algorithm described in Taddy (2011; Section 3.1). 
	In particular, 	we obtain loadings \code{loadings[,k]} as the correlation between
	\code{X} and factors \code{factors[,k]}, where \code{factors[,1]} is initialized      
	at \code{normalize(as.numeric(y))} and subsequent factors are orthogonal 
	to the k'th pls direction, \code{directions[,k]=X\%*\%loadings[,k]}. }
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
summary( fit <- pls(freq(congress109Counts), congress109Ideology$repshare, K=3) )
plot(fit, pch=21, bg=c(4,3,2)[congress109Ideology$party])
predict(fit, newdata=freq(congress109Counts[c(68,388),]))

data(we8there)
summary( fit <- pls(freq(we8thereCounts),as.factor(we8thereRatings$Overall)) )
plot(fit, col=c(2,2,2,3,3))
}