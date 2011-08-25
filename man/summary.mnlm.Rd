\name{summary.mnlm}
\alias{summary.mnlm}
\alias{print.mnlm}
\title{
mnlm summary
}
\description{ Summary function for Multinomial Logistic Regression }
\usage{
\method{summary}{mnlm}( object, y=NULL, ... )
}
\arguments{
  \item{object}{An output object from the \code{mnlm} function.}
  \item{y}{A possible response (sentiment) variable of interest in an inverse regression setting.}
  \item{...}{Additional unused arguments. }
}
\details{ A short summary function for mnlm objects.
	 }
\value{ A printout describing the regression coefficients (dimension and sparsity) along with 
	some within-sample correlations or error rates, depending on response and covariate formats.
	If \code{y} is specified and code{ncol{object$counts}>2},
	 we also print the goodness-of-fit R2 for least-squares linear regression of
	\code{y} onto the sufficient reduction \code{freq(X)\%*\%loadings}. 
}
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  mnlm, congress109, we8there
}
