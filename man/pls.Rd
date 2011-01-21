\name{pls}
\alias{pls}
\title{ Partial Least Squares }
\description{ A simple partial least squares procedure. }
\usage{ pls(F, y, K=1) }
\arguments{
  \item{F}{ The covariate matrix.  It is usually best that this is first standardized 
   with the \code{normalize} function to have mean-zero and variance-one columns}
   \item{y}{ The response vector.  }
   \item{K}{ The number of desired PLS directions. }
}
\value{ A pls object list with the following entries
	\item{y}{ The response vector. }
	\item{F}{The covariate matrix.}
	\item{z}{The pls directions \code{F\%*\%phi}.}
	\item{v}{Response factors.}
	\item{yhat}{\code{K} columns of fitted values for each number of directions.}
	\item{fwdmod}{ The \code{lm} object from forward regression \code{lm(y~z)}. }
}
\details{ Fits the Partial Least Squares algorithm described in Taddy (2011; Section 3.1). 
	In particular, 	we obtain loadings \code{phi[,k]} as the correlation between
	\code{F} and factors \code{v[,k]}, where \code{v[,1]} is initialized      
	at \code{y} and subsequent factors are orthogonal 
	to the k'th pls direction, \code{z[,k]=F\%*\%phi[,k]}. }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}

Wold, H. (1975), \emph{Soft modeling by latent variables: The nonlinear iterative partial least squares approach}. 
In Perspectives in Probability and Statistics, Papers in Honour of M.S. Bartlett.

}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  plot.pls, normalize, freq
}
\examples{
data(congress109)
F <- normalize(freq(congress109Counts))
y <- normalize(congress109Ideology$repshare)
fit <- pls(F, y, K=4)
plot(fit, pch=21, bg=c(4,3,2)[congress109Ideology$party])
 }