\name{topics}
\alias{topics}
\title{
Estimation for Topic Models
}
\description{ MAP estimation of Topic models  }
\usage{
topics(counts, K, alpha=NULL, initheta=NULL, tol=0.1, 
bf=FALSE, kill=2, ord=TRUE, verb=1, ...)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} phrases/categories 
    for \code{nrow(counts)} documents/observations. Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.}
  \item{K}{The number of latent topics.  If \code{length(K)>1}, 
  \code{topics} will find the Bayes factor (vs a null single topic model) for each element and return parameter
  estimates for the highest probability K.}
  \item{alpha}{ Optional Dirichlet prior concentration parameter for topic-phrase probabilities.  Defaults to \code{1/(K*ncol(counts))}.}
  \item{initheta}{ Optional start-location for \eqn{[\theta_1 ... \theta_K]}, the topic-phrase probabilities.  
  Dimensions must accord with the smallest element of \code{K}.  
  If \code{initheta=NULL}, the initial estimates are built by incrementally adding topics onto a simple two-topic base. }
  \item{tol}{Convergence tolerance.  The optimization terminates when log posterior increase over a full parameter-set update
  is less than \code{tol}.}
  \item{bf}{An indicator for whether or not to calculate the Bayes factor for univariate \code{K}.  
  If \code{length(K)>1}, this is ignored and Bayes factors are always calculated. }
  \item{kill}{For choosing from multiple \code{K} numbers of topics (evaluated in increasing order),
   the search will stop after \code{kill} consecutive drops in the corresponding Bayes factor.  Specify \code{kill=0}
   if you want Bayes factors for all elements of \code{K}. }
 \item{ord}{If \code{TRUE}, the returned topics (columns of \code{theta}) will be ordered by decreasing usage (i.e., 
 by decreasing \code{colSums(omega)}).}
  \item{verb}{A switch for controlling printed output.  \code{verb > 0} will print something, with the level of detail
  increasing with \code{verb}.}
  \item{...}{Additional arguments to the undocumented internal \code{tpx*} functions. }
 }
\details{A latent topic model represents each i'th document's term-count vector \eqn{X_i} 
(with \eqn{\sum_{j} x_{ij} = m_i} total phrase count)
as having been drawn from a mixture of \code{K} multinomials, each parameterized by topic-phrase
probabilities \eqn{\theta_i}, such that 
\deqn{X_i \sim MN(m_i, \omega_1 \theta_1 + ... + \omega_K\theta_K).}
We assign a K-dimensional Dirichlet(1/K) prior to each document's topic weights 
\eqn{[\omega_{i1}...\omega_{iK}]}, and the prior on each \eqn{\theta_k} is Dirichlet with concentration \eqn{\alpha}.
The \code{topics} function uses quasi-newton accelerated EM, augmented with sequential quadratic programming 
for conditional \eqn{\Omega | \Theta} updates, to obtain MAP estimates for the topic model parameters. 
We also provide Bayes factor estimation, from marginal likelihood
calculations based on a Laplace approximation around the converged MAP parameter estimates.  If  input \code{length(K)>1}, these
Bayes factors are used for model selection. Full details are in Taddy (2011). }
\note{ Estimates
for each \eqn{\omega_i} are actually MAP in the (K-1)-dimensional logit transformed natural exponential family parameters space, while 
\eqn{\Theta} is MAP in the original probability space.  }
\value{
 An \code{topics} object list with entries
\item{K}{The number of latent topics estimated.  If input \code{length(K)>1}, 
on output this is a single value corresponding to the model with the highest Bayes factor. }
\item{theta}{The \code{ncol{counts}} by \code{K} matrix of estimated topic-phrase probabilities.}
\item{omega}{The \code{nrow{counts}} by \code{K} matrix of estimated document-topic weights. }
\item{BF}{The log Bayes factor for each number of topics in the input \code{K}, against a null single topic model.}
\item{residuals}{Standardized residuals (approximately \eqn{\sim N(0,1)}) \emph{for only nonzero count entries}. } 
\item{X}{The input count matrix, in \code{simple_triplet_matrix} format.}
}
\references{ Taddy (2011), \emph{Estimation of Topic Models}.}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  plot.topics, summary.topics, predict.topics, wsjIBM, congress109, we8there
}
\examples{
## see wsjibm, congress109, and we8there for data examples

## Simulation Parameters
K <- 10
n <- 150
p <- 200
omega <- t(rdir(n, rep(1/K,K)))
theta <- rdir(K, rep(1/p,p))

## Simulated counts
Q <- omega\%*\%t(theta)
counts <- matrix(ncol=p, nrow=n)
totals <- rpois(n, 250)
for(i in 1:n){ counts[i,] <- rmultinom(1, size=totals[i], prob=Q[i,]) }

## Bayes Factor model selection (should choose K or nearby)
simselect <- topics(counts, K=5:15, tol=.01, verb=1) 
print(simselect$K)

## MAP fit for given K=10
simfit <- topics(counts, K=K, verb=2)

## Adjust for label switching and plot the fit (color by topic)
toplab <- rep(0,K)
for(k in 1:K){ toplab[k] <- which.min(colSums(abs(simfit$omega-omega[,k]))) }
par(mfrow=c(1,2))
tpxcols <- matrix(rainbow(K), ncol=nrow(theta), byrow=TRUE)
plot(theta,simfit$theta[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
plot(omega,simfit$omega[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
title("True vs Fitted Values (color by topic)", outer=TRUE, line=-2)

## The S3 method plot functions
par(mfrow=c(1,2))
plot(simfit, lgd.K=2)
plot(simfit, type="resid")
}
