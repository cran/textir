\name{mnlm}
\alias{mnlm}
\title{
Estimation for high-dimensional Multinomial Logistic Regression
}
\description{
MAP estimation of multinomial logistic regression models.  
  }
\usage{
mnlm(counts, covars, normalize=TRUE, penalty=c(shape=1,rate=1/2), 
                     start=NULL, tol=1e-2, bins=0, verb=FALSE, ...)
1}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} or \code{nlevel(counts)} categories 
    for \code{nrow(counts)} observations.  This can be a \code{matrix}, a vector of response factors, 
   or a \code{simple_triplet_matrix} (as defined in the \code{slam} package).  
   Refer to the details for a model identification note.
}
  \item{covars}{
    A \code{matrix} or \code{simple_triplet_matrix} of \code{ncol(covars)} covariate values for each of the \code{nrow(counts)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.   
}
 \item{normalize}{
    Whether or not to normalize the covariates. Default is \code{TRUE}.  If \code{covars} is a \code{matrix}, 
    this will scale the inputs to have mean zero and standard deviation of one.   If \code{covars} is a \code{simple_triplet_matrix}, 
    we assume that you want to stay in sparse format; hence the inputs are scaled to have sd = 1 but left unshifted.
}
\item{penalty}{
  This input argument is a vector of length 2 containing \eqn{[s, r]} 
  -- shape "\eqn{s}" and rate "\eqn{r}" -- 
  parameters for the Gamma prior on L1 (lasso) penalty \eqn{\lambda},
  such that \eqn{E\lambda = s/r}.   Refer to the details section for additional 
  information on this gamma-lasso specification. 
  The default is appropriate for normalized covariates. 

  Additionally, you can specify a normal (ridge) prior with variance
  1/\code{rate} by setting the shape to zero (i.e. with
  \code{penalty=c(0,rate)}), set \code{penalty} to a single fixed
  value of \eqn{\lambda>0}, or fix coefficients at \code{start} by
  giving a penalty of -1.

  Finally, \code{penalty} can also be defined as a list with  elements containing unique specification
  for each column of the design matrix (including the intercept).   
 }
\item{start}{
    An optional initial guess for the full  \code{ncol(covars)+1} by \code{ncol(counts)} 
    matrix of regression coefficients (including the intercept).   
    Under the default \code{start=NULL},
    the intercept is a logit transform of mean phrase frequencies and 
    coefficients are the correlation between \code{covars} 
    and  \code{freq(counts)}.
}
\item{tol}{
  Optimization convergence tolerance for the improvement on the 
  un-normalized negative log posterior over a single full parameter sweep.
}
\item{bins}{ For faster inference on large data sets 
	     (or just to collapse observations across levels for factor covariates), 
	     you can specify the number of \code{bins} for step-function 
	     approximations to the columns of \code{covars}.  Counts are then collapsed 
	     across levels of the interaction between columns of the resulting 
	     discrete covariate matrix, typically resulting in a smaller number of observations for estimation.}
\item{verb}{
   Control for print-statement output.  \code{TRUE} prints some initial info 
   and updates every iteration.
}
\item{...}{
   Additional undocumented arguments to internal functions.   }
}
\details{
 Finds the posterior mode for multinomial logistic regression
 parameters using cyclic coordinate descent.  This is designed to be
 useful for inverse regression analysis of sentiment in text, where
 the multinomial response is quite large, but should be 
 useful for any large-scale logistic regression.

 For binomial response, the first category is assumed null.  For
 multinomial response, the model is identified by placing a Normal(0,1)
 prior on the intercepts (this can be changed via the list
 specification for \code{penalty}).

 Coefficient penalization is based upon the precision parameters
 \eqn{\lambda} of independent Laplace priors on each non-intercept
 regression coefficient.  Here, the Laplace density is \eqn{p(z) =
 (\lambda/2)exp[-\lambda|z|]}, with variance \eqn{2/\lambda}. Via the
 \code{penalty} argument, this precision is either fixed, which
 corresponds to the L1 penalty \eqn{\lambda|z|}, or it is assigned a
 \eqn{Gamma(s, r)} prior and estimated jointly with the coefficient,
 which corresponds to the `gamma-lasso' non-convex penalty \eqn{s*log[1 + |z|/r]}.

 In the case of gamma-lasso estimation, prior variance
 \eqn{s/r^2 = E\lambda/r} controls the degree of penalty curvature.
 In the case that the variance is large relative to the amount of
 information in the likelihood, the posterior can become multimodal.
 Since this leads to unstable optimization and less meaningful MAP
 estimates, \code{mnlm} will warn and automatically double \eqn{r} and
 \eqn{s} until obtaining a concave posterior.  If the resulting prior
 precision is higher than you would like, it may be worth the
 computational effort to integrate over penalty uncertainty in mean,
 rather than MAP, estimation; the \code{reglogit} package is available
 for such inference in binomial regression settings.

 Additional details are available in Taddy (2012).   }
\value{
 An \code{mnlm} object list with entries
\item{intercept}{The intercept estimates for each phrase (\eqn{\alpha}).}
\item{loadings}{A \code{simple_triplet_matrix} of estimates for coefficients (\eqn{\Phi}) on the scale fitted (possibly normalized) covariates.}
\item{counts}{\code{simple_triplet_matrix} form of the \code{counts} input matrix}
\item{X}{If \code{bins>0}, the binned counts matrix used for analysis.}
\item{covars}{The input covariates, possibly normalized.}
\item{V}{If \code{bins>0}, the binned (and possibly normalized) covariate \code{simple_triplet_matrix} used for analysis.}
\item{penalty}{The penalty specification upon convergence.}
\item{normalized}{The input \code{normalize} indicator.}
\item{binned}{An indicator for whether the observations was binned.}
\item{covarMean}{If \code{normalize=TRUE}, the amount covariates were shifted (original means for matrix \code{covars}, 0 for sparse stm \code{covars}).  Otherwise empty.}
\item{covarSD}{If \code{normalize=TRUE}, the original covariate standard deviations.  Otherwise empty.}
\item{prior}{The penalty prior (gamma hyperparameters, or fixed laplace scale, or normal precision).}
\item{fitted}{Fitted count expectations.  With binomial response, this is a vector of fitted probabilities. 
	For multinomial response, it is a simple triplet matrix if of fitted probabilities 
	ONLY for non-zero count observations (and with empty entries for zero count observations).}
}
\references{
Taddy (2012), \emph{Multinomial Inverse Regression for Text Analysis}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  congress109, we8there, plot.mnlm, summary.mnlm, predict.mnlm
}
\examples{

### See congress109 and we8there for more real data examples

### Bernoulli simulation; re-run to see sampling variability ###
n <- 100
v <- rnorm(n)
p <- (1+exp(-(v*2)))^{-1} 
y <- rbinom(n, size=1, prob=p)

## fit the logistic model
summary( fit <- mnlm(y, v, verb=TRUE) )
par(mfrow=c(1,2))
plot(fit)

## use predict to see fitted probabilities (could also just use fit$fitted)
phat <-  predict(fit, newdata=matrix(v,ncol=1))
plot(p, phat, pch=21, bg=c(2,4)[y+1], xlab="true probability", ylab="fitted probability")

### Ripley's Cushing Data ###

## see help(Cushings) for data
library(MASS)
data(Cushings)
train <- Cushings[Cushings$Type != "u",]
newdata <- as.matrix(Cushings[Cushings$Type == "u", 1:2])

## fit, summarize, predict, and plot
fit <- mnlm(counts=factor(train$Type), covars=train[,1:2])
summary(fit)
round(coef(fit),2)
predict(fit, newdata)
par(mfrow=c(1,1))
plot(fit)

 }
