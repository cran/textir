\name{mnlm}
\alias{mnlm}
\title{
Estimation for high-dimensional Multinomial Logistic Regression
}
\description{
MAP estimation of multinomial logistic regression models.  
  }
\usage{
mnlm(counts, covars, normalize=FALSE, penalty=c(shape=1,rate=1), 
                     start=NULL, tol=0.01, bins=0, verb=FALSE, ...)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} categories 
    for \code{nrow(counts)} individuals/observations.  This can be a \code{matrix}, or a vector of response factors, 
    but for most text-analysis applications should be a \code{simple_triplet_matrix}, as defined in the \code{slam} package.  
    See the details for a note on null-category parameterization.
}
  \item{covars}{
    A matrix of \code{ncol(covars)} covariate values for each of the \code{nrow(counts)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.   
}
 \item{normalize}{
    Whether or not to normalize the covariates to have mean zero and standard deviation of 1/2.  
}
\item{penalty}{
  This input argument is either a single fixed value of the L1 (lasso) penalty \eqn{\lambda>0}, 
  or a vector of length 2 containing \eqn{[s, r]} -- shape '\eqn{s}' and rate '\eqn{r}' parameters for the Gamma prior on \eqn{\lambda}, 
  such that \eqn{E\lambda =s/r}.  Refer to the details section for additional information.
  \code{penalty} 
  can also be a list containing penalty specification
  for each column of the 
  design matrix (i.e. including the intercept).  Using this list specification, 
  you can also fix a column at the \code{start} values by giving a penalty of -1 for that coefficient.
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
\item{bins}{ For faster inference on large data sets (or just to collapse observations across levels for factor covariates), 
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
 Finds the posterior mode for multinomial logistic regression parameters using cyclic coordinate descent. 
 This is designed to be useful for inverse regression analysis of sentiment in 
 text, where the multinomial response is quite large.  It should be generally useful 
 for any large-scale multinomial logistic regression, but is optimized for the large-response setting 
 (e.g., \code{counts} are treated as sparse while \code{covars} are dense). 

 The model is identified by fixing coefficients at zero for a
 specified null category.  With binomial response, the first category
 is assumed null.  For multinomial response dimension greater than
 two, each response vector is augmented with a null category count of
 zero, such that the linear model equations can then be interpreted as
 log odds of each response category against a very rare null event
 with covariate independent probability. This specification is
 designed to work well for high-dimension response (e.g. text), our
 motivating application, but should work in a variety of
 settings. Fitted probabilities and those obtained using
 \code{predict.mnlm} are corrected to condition on response not coming
 from this null category.

 Coefficient penalization is based upon the precision parameters
 \eqn{\lambda} of independent Laplace priors on each non-intercept
 regression coefficient.  Here, the Laplace density is \eqn{p(z) = (\lambda/2)exp[-\lambda|z|]},
 with variance \eqn{2/\lambda}. Via the \code{penalty} argument, this precision is either fixed, 
 which corresponds to the L1 penalty \eqn{\lambda|z|}, or it is assigned a \eqn{Gamma(s, r)} 
 prior and estimated jointly with the coefficient, which corresponds to the non-convex penalty \eqn{s*log[1 + |z|/r]}.  

 In the case of joint penalty-coefficient estimation, prior precision \eqn{r^2/s = r/E\lambda} controls the degree of penalty curvature.  
 In the case that the precision is small relative to the amount of information in the likelihood, the posterior can become multimodal.  
 Since this leads to unstable optimization and less meaningful MAP estimates, 
 \code{mnlm} will warn and automatically double \eqn{r} and \eqn{s}  until obtaining a concave posterior. 
 If the resulting prior precision is higher than you would like, 
 it may be worth the computational effort to integrate over penalty uncertainty in mean, 
 rather than MAP, estimation; the \code{reglogit} package is available for such inference in binomial regression settings.

 Additional details are available in Taddy (2011).   }
\value{
 An \code{mnlm} object list with entries
\item{intercept}{The intercept estimates for each phrase (\eqn{\alpha}).}
\item{loadings}{The intercept estimates for each phrase (\eqn{\Phi}).}
\item{counts}{\code{simple_triplet_matrix} form of the \code{counts} input matrix}
\item{X}{If \code{bins>0}, the binned counts matrix used for analysis.}
\item{covars}{The input covariates, possibly normalized.}
\item{V}{If \code{bins>0}, the binned (and possibly normalized) covariate matrix used for analysis.}
\item{penalty}{The penalty specification upon convergence.}
\item{normalized}{An indicator for whether the covariates were normalized to have mean 0 and standard deviation 1/2.}
\item{binned}{An indicator for whether the observations was binned.}
\item{covarMean}{If \code{normalize=TRUE}, the original covariate means.  Otherwise empty.}
\item{covarSD}{If \code{normalize=TRUE}, the original covariate standard deviations.  Otherwise empty.}
\item{prior}{The penalty prior (gamma distribution parameters, or fixed values).}
\item{lambda}{Posterior MAPs (or just the fixed input values) for each coefficient's penalty.}
\item{gradient}{Loading parameter gradients at convergence (set to 0 for zero coefficients that lie at the point of a local minimum).}
\item{fitted}{Fitted count expectations.  With binary response, this is a vector of fitted probabilities. 
	For binomial or multinomial response, it is a simple triplet matrix
	 with empty entries for zero count observations.}
}
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
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
summary( fit <- mnlm(y, v) )
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
fit <- mnlm(counts=factor(train$Type), covars=train[,1:2], normalize=TRUE)
summary(fit)
predict(fit, newdata)
par(mfrow=c(1,1))
plot(fit)

 }
