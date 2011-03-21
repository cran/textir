\name{mnlm}
\alias{mnlm}
\alias{mncheck}
\title{
Estimation for high-dimensional Multinomial Logistic Regression
}
\description{
MAP estimation of Multinomial logistic regression models.  
  }
\usage{
mnlm(counts, covars, normalize=FALSE, lambda=NULL, start=NULL,
                  tol=0.1, tmax=1000, delta=1, dmin=0, bins=0, verb=TRUE)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} categories 
    for \code{nrow(counts)} individuals/observations.  This can be a standard \code{matrix}, 
    but for most text-analysis applications should be a \code{simple_triplet_matrix}, as defined in the \code{slam} package.
}
  \item{covars}{
    A matrix of \code{ncol(covars)} covariate values for each of the \code{nrow(counts)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.   
}
 \item{normalize}{
    Whether or not to normalize the covariate matrix to have mean zero and variance one.  
}
\item{lambda}{
  Either a single fixed value, or a vector of length 2 giving the gamma hyperprior prior shape and rate parameters 
  (e.g., \code{c(s=2, r=2)} ).  Here, lambda (>0) is a joint 
  scale parameter for the Laplace prior on each non-intercept 
  regression coefficient, parametrized
   \deqn{p(b) =  (\lambda/2)exp[-\lambda|b|].}  
  Refer to Taddy (2010) for details.
  Under the default \code{lambda=NULL}, 
  we use \code{shape=rate=p/4}.
}
\item{start}{
    An initial guess for the full \code{ncol(counts)} by \code{ncol(covars)+1}
    matrix of regression coefficients.   
    Under the default \code{start=NULL},
    the intercept is a logit transform of mean phrase frequencies and 
    coefficients are the correlation between \code{covars} 
    and  \code{freq(counts)}.
}
  \item{tol}{
  Optimization convergence tolerance for the improvement on the 
  un-normalized negative log posterior over a single full parameter sweep.
}
  \item{tmax}{
    The maximum number of optimization iterations. 
}
\item{delta}{
   An initial step size for the least upper bound approximation to parameter information;
   implies a starting trust region of \code{2*delta}.  
}
\item{dmin}{ Minimum trust region delta.}
\item{bins}{ For faster inference on large data sets (or just to collapse observations across levels for factor covariates), 
	     you can specify the number of \code{bins} for step-function 
	     approximations to the columns of \code{covars}.  Counts are then collapsed 
	     across levels of the interaction between columns of the resulting 
	     discrete covariate matrix, typically resulting in a smaller number of observations for estimation.}
\item{verb}{
   Level of print-statement output.  \code{TRUE} prints some initial info 
   and updates every iteration.
}
}
\details{
 Finds the posterior mode for multinomial logistic regression parameters using cyclic coordinate descent. 
This is designed to be useful for inverse regression analysis of sentiment in 
text, where the multinomial response is quite large, but should be generally useful 
for any large-scale multinomial logistic regression. We allow for joint estimation of
regression coefficients  and a Laplace regularization penalty.  Regression coefficients are
identified by augmenting each response vector with a null count of 0.01
and assuming zero coefficients for this category.
Full details are available in Taddy (2011).   }
\value{
 An \code{mnlm} object list with entries
\item{intercept}{The intercept estimates for each phrase (\eqn{\alpha}).}
\item{loadings}{The intercept estimates for each phrase (\eqn{\phi}).}
\item{counts}{\code{simple_triplet_matrix} form of the \code{counts} input matrix}
\item{X}{The design matrix used for analysis; includes an added null column and may have merged observations from \code{counts}.}
\item{covars}{The input covariates, possibly normalized.}
\item{V}{The covariate matrix used for analysis; possibly normalized or binned, and including the intercept}
\item{covarMean}{If \code{normalize=TRUE}, the original covariate means.  Otherwise empty.}
\item{covarSD}{If \code{normalize=TRUE}, the original covariate standard deviations.  
Otherwise empty.}
\item{maplam}{An indicator for whether the regularization penalty was estimated.}
\item{lampar}{Parameters (init, shape, rate) for the regularization penalty.}
\item{lambda}{The path of \code{lambda} estimates.}
\item{delta}{The trust region deltas upon convergence.}
\item{L}{The unnormalized negative log posterior at each iteration.}
\item{niter}{The number of iterations}
\item{tol,tmax}{Convergence parameters, unchanged from input.}
\item{start}{The initial coefficient estimates.}
}
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}

}
\note{ \code{mncheck} is an undocumented utility for discretizing \code{covars} as
specified by the \code{bins} argument}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  congress109, we8there, plot.mnlm, pls
}
\examples{

## See congress109 and we8there for real data examples

## Binomial simulation; re-run to see sampling variability
n <- 20
size <- 10
v <- rnorm(n)
p <- (1+exp(-(1 + v*2)))^{-1} 
y <- rbinom(n, size=size, prob=p)
counts <- cbind(size-y, y)

## fit the logistic model
fit <- mnlm(counts, v)

## extract fitted probabilities
eta <- fit$intercept + fit$loadings\%*\%t(v)
q0 <- 1/(1+colSums(exp(eta))) # null category
phat <- t(exp(eta))*( q0/(1-q0) )
plot(p, phat[,2], pch=21, bg=rainbow(n), 
	xlab="true", ylab="fitted", main="binomial probability")
 }
