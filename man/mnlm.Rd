\name{mnlm}
\alias{mnlm}
\title{
Estimation for high-dimensional Multinomial Logistic Regression
}
\description{
MAP estimation of multinomial logistic regression models.  
  }
\usage{
mnlm(counts, covars, normalize=FALSE, penalty=c(shape=1,rate=0.2), start=NULL, 
                  tol=0.1, tmax=1000, delta=1, dmin=0, bins=0, verb=FALSE)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} categories 
    for \code{nrow(counts)} individuals/observations.  This can be a \code{matrix}, or a vector of response factors, 
    but for most text-analysis applications should be a \code{simple_triplet_matrix}, as defined in the \code{slam} package.
}
  \item{covars}{
    A matrix of \code{ncol(covars)} covariate values for each of the \code{nrow(counts)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.   
}
 \item{normalize}{
    Whether or not to normalize the covariate matrix to have mean zero and variance one.  
}
\item{penalty}{
  Either a single fixed value, or a vector of length 2 giving the gamma hyperprior
 shape and rate parameters (for expected penalty of \code{shape/rate}).
  Here, the penalty (\eqn{\lambda>0}) is a 
  scale parameter for the Laplace prior on each non-intercept 
  regression coefficient, parametrized
   \deqn{p(b) =  (\lambda/2)exp[-\lambda|b|].}  
  Refer to Taddy (2010) for details.  Alternatively, \code{penalty} 
  can be a list with with individual penalties
  (single value or length-2 prior vector) for each column of the 
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
   Control for print-statement output.  \code{TRUE} prints some initial info 
   and updates every iteration.
}
}
\details{
 Finds the posterior mode for multinomial logistic regression parameters using cyclic coordinate descent. 
This is designed to be useful for inverse regression analysis of sentiment in 
text, where the multinomial response is quite large, but should be generally useful 
for any large-scale multinomial logistic regression. We allow for joint estimation of
regression coefficients  and Laplace regularization penalties.  
For response dimension greater than two, the regression coefficients are
identified by augmenting each response vector with a null count of 1/1000 of each observation's total
and assuming zero coefficients for this category.  With binomial response, the first category is assumed null.
Full details are available in Taddy (2011).   }
\value{
 An \code{mnlm} object list with entries
\item{intercept}{The intercept estimates for each phrase (\eqn{\alpha}).}
\item{loadings}{The intercept estimates for each phrase (\eqn{\phi}).}
\item{counts}{\code{simple_triplet_matrix} form of the \code{counts} input matrix}
\item{X}{If \code{bins>0}, the binned counts matrix used for analysis.}
\item{covars}{The input covariates, possibly normalized.}
\item{V}{If \code{bins>0}, the binned (and possibly normalized) covariate matrix used for analysis.}
\item{normalized}{An indicator for whether the covariates were normalized.}
\item{binned}{An indicator for whether the observations was binned.}
\item{covarMean}{If \code{normalize=TRUE}, the original covariate means.  Otherwise empty.}
\item{covarSD}{If \code{normalize=TRUE}, the original covariate standard deviations.  
Otherwise empty.}
\item{L}{The unnormalized negative log posterior at each iteration.}
\item{residuals}{Standardized Pearson residuals, for only nonzero count entries.  
			  In simple triplet matrix format, with empty entries for zero count observations.}
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
fit <- mnlm(counts=factor(train$Type), covars=as.matrix(train[,1:2]))
summary(fit)
predict(fit, newdata)
par(mfrow=c(1,1))
plot(fit)

 }
