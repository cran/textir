\name{topics}
\alias{topics}
\title{
Estimation for Topic Models
}
\description{ MAP estimation of Topic models  }
\usage{
topics(counts, K, shape=NULL, initopics=5, tol=0.01, admix=TRUE, grp=NULL,
bf=FALSE, kill=2, ord=TRUE, verb=1, ...)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} phrases/categories 
    for \code{nrow(counts)} documents/observations. Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.}
  \item{K}{The number of latent topics.  If \code{length(K)>1}, 
  \code{topics} will find the Bayes factor (vs a null single topic model) for each element and return parameter
  estimates for the highest probability K.}
  \item{shape}{ Optional argument to specify the Dirichlet prior concentration parameter as \code{1+shape} for topic-phrase probabilities.  Defaults to \code{1/(K*ncol(counts))}.  For fixed single \code{K}, this can also be a \code{ncol(counts)} 
  		by \code{K} matrix of unique shapes for each topic element.}
  \item{initopics}{ Either an integer, or an optional start-location for \eqn{[\theta_1 ... \theta_K]}, the topic-phrase probabilities.  
  		    In the latter case, dimensions must accord with the smallest element of \code{K}.  
  		    If an integer, the initial estimates are built in \code{initopics} steps by incrementally adding topics. }
  \item{tol}{ Convergence tolerance.  The optimization terminates when absolute change over a full parameter-set update
  is less than \code{tol}.}
  \item{admix}{ An indicator for whether to fit the standard mixed-membership `admixture' topic model (default), 
                 or a simpler finite mixture of multinomials. }
  \item{grp}{ If \code{!admix}, an optional vector of group labels for each document, corresponding to blocked mixture proportions.  Interpreted as a \code{factor}. If \code{NULL}, every document is assumed to have the same prior mixture component probabilities. }  
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
for each \eqn{\omega_i} are actually functions of the MAP (K-1)-dimensional 
logit transformed natural exponential family parameter, while 
\eqn{\Theta} is MAP in the original probability space.  }
\value{
 An \code{topics} object list with entries
\item{K}{The number of latent topics estimated.  If input \code{length(K)>1}, 
on output this is a single value corresponding to the model with the highest Bayes factor. }
\item{theta}{The \code{ncol{counts}} by \code{K} matrix of estimated topic-phrase probabilities.}
\item{omega}{The \code{nrow{counts}} by \code{K} matrix of estimated document-topic weights. }
\item{BF}{The log Bayes factor for each number of topics in the input \code{K}, against a null single topic model.}
\item{dispersion}{The mean of squared positive residuals, an estimate of the multinomial dispersion parameter, 
for each number of topics in the input \code{K}.}
\item{residuals}{Standardized residuals (approximately \eqn{\sim N(0,1)}) \emph{for only nonzero count entries}. } 
\item{X}{The input count matrix, in \code{simple_triplet_matrix} format.}
\item{admix}{The type of model.  Unchanged from input.}
\item{grp}{Optional grouping under \code{!admix}.  Unchanged from input.}
}
\references{ Taddy (2011), \emph{Estimation of Topic Models}.}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  plot.topics, summary.topics, predict.topics, wsjibm, congress109, we8there
}
\examples{
## see wsjibm, congress109, and we8there for data examples

## Simulation Parameters
K <- 10
n <- 150
p <- 250
omega <- t(rdir(n, rep(1/K,K)))
theta <- rdir(K, rep(1/p,p))

## Simulated counts
Q <- omega\%*\%t(theta)
counts <- matrix(ncol=p, nrow=n)
totals <- rpois(n, 200)
for(i in 1:n){ counts[i,] <- rmultinom(1, size=totals[i], prob=Q[i,]) }

## Bayes Factor model selection (should choose K or nearby)
summary(simselect <- topics(counts, K=K+c(-5:5)), nwrd=0)

## MAP fit for given K
simfit <- topics(counts,  K=K)

## Adjust for label switching and plot the fit (color by topic)
toplab <- rep(0,K)
for(k in 1:K){ toplab[k] <- which.min(colSums(abs(simfit$omega-omega[,k]))) }
par(mfrow=c(1,2))
tpxcols <- matrix(rainbow(K), ncol=ncol(theta), byrow=TRUE)
plot(theta,simfit$theta[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
plot(omega,simfit$omega[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
title("True vs Fitted Values (color by topic)", outer=TRUE, line=-2)

## The S3 method plot functions
par(mfrow=c(1,2))
plot(simfit, lgd.K=2)
plot(simfit, type="resid")

## Same thing, but now with a simple mixture model 
## (ie, not the mixed membership admixture model)

## mixture simulation
n <-5000
theta <- rdir(K, rep(1/p,p))
omega <- (1:K)/sum(1:K)
grp <- apply(rmultinom(n, 1, omega),2,which.max)
counts <- matrix(ncol=p, nrow=n)
totals <- rpois(n,200)
for(i in 1:n){ counts[i,] <- rmultinom(1, size=totals[i], prob=theta[,grp[i]]) }

## fit the mixture model
summary( mxt <- topics(counts, K=10, admix=FALSE), nwrd=0)

## Adjust for label switching and plot the fit (color by topic)
par(mfrow=c(1,1))
toplab <- rep(0,K)
for(k in 1:K){ toplab[k] <- which.min(colSums(abs(mxt$theta-theta[,k]))) }
tpxcols <- matrix(rainbow(K), ncol=ncol(theta), byrow=TRUE)
plot(theta, mxt$theta[,toplab], ylab="fitted values", pch=21, 
     bg=tpxcols, main="True vs Fitted Values (color by topic)")



}
