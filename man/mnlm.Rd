\name{mnlm}
\alias{mnlm}
\title{
Estimation for high-dimensional Multinomial Logistic Regression
}
\description{
MAP estimation of Multinomial logistic regression models.  
  }
\usage{
mnlm(counts, covars, normalize=FALSE, lambda=NULL, start=NULL,
                  tol=0.1, tmax=1000, delta=1, dmin=0, verb=TRUE)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} categories 
    for \code{nrow(counts)} individuals/observations.
}
  \item{covars}{
    A matrix of \code{ncol(covars)} covariate values for each of the \code{nrow(counts)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.   
}
 \item{normalize}{
    Whether or not to normalize the covariate matrix to have mean zero and variance one.  
}
\item{lambda}{
  Either a single fixed value, or a vector of length 3 giving initial lambda followed 
  by gamma prior shape and rate parameters 
  (e.g., \code{c(init=sqrt(2), shape=2, rate=2)} ).  Here, lambda (>0) is a joint 
  scale parameter for the Laplace prior on each non-intercept 
  regression coefficient, parametrized
   \deqn{p(b) =  (\lambda/2)exp[-\lambda|b|].}  
  Refer to Taddy (2010) for details.
  Under the default \code{lambda=NULL}, 
  we initialize at \code{lambda=sqrt(2)} with \code{shape=rate=(ncol(counts)-1)/4}.
}
\item{start}{
    An initial guess for the full \code{ncol(counts)} by \code{ncol(covars)+1}
    matrix of regression coefficients.  The first row (usually set to zero) corresponds to
    our "null" category and is not updated during the model fit.  
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
regression coefficients  and a Laplace regularization penalty.
Full details are available in Taddy (2011). }
\value{
 An \code{mnlm} object list with entries
\item{coef}{Regression coefficient estimates, including the intercept column}
\item{L}{The unnormalized negative log posterior at each iteration.}
\item{niter}{The number of iterations}
\item{X}{\code{simple_triplet_matrix} form of the \code{counts} input matrix}
\item{V}{Covariate matrix (including the intercept)}
\item{covarMean}{If \code{normalize=TRUE}, the original covariate means.  Otherwise empty.}
\item{covarSD}{If \code{normalize=TRUE}, the original covariate standard deviations.  
Otherwise empty.}
\item{maplam}{An indicator for whether the regularization penalty was estimated.}
\item{lampar}{Parameters (init, shape, rate) for the regularization penalty.}
\item{lambda}{The path of \code{lambda} estimates.}
\item{delta}{The trust region deltas upon convergence.}
\item{tol,tmax}{Convergence parameters, unchanged from input.}
\item{start}{The initial coefficient estimates.}
}
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}

}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  congress109, we8there, plot.mnlm, pls
}
\examples{

## See congress109 and we8there for data examples

n <- 1000
size <- 10

v <- rnorm(n)
covars <- cbind(v)

p1 <- (1+exp(-(1 + v*2)))^{-1} 
p2 <- (1+exp(-(1 - v*2)))^{-1}
y1 <- rbinom(n, size=size, prob=p1)
y2 <- rbinom(n, size=size, prob=p2)

# binomial
counts <- cbind(size-y1, y1)
print(mnlm(counts, covars)$coef)

# 3 category multinomial
counts <- cbind(2*size-y1-y2,y1,y2)
print(mnlm(counts, covars)$coef)

 }
