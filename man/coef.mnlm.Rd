\name{coef.mnlm}
\alias{coef.mnlm}
\alias{coefficients.mnlm}
\title{
mnlm coefficients
}
\description{ Coefficients for Multinomial Regression }
\usage{
\method{coef}{mnlm}( object, origscale=TRUE, ... )
}
\arguments{
  \item{object}{An output object from the \code{pls} function.}
  \item{origscale}{Whether to ouput coefficients on the original covariate scale (i.e. before possible normalization). 
  Default is \code{TRUE}, and \code{origscale=FALSE} just outputs raw loadings for the fitted model}
  \item{...}{Additional unused arguments. }
}
\value{ A \code{ncol(object$covars)+1} by \code{ncol(object$counts)} (or by 1 for binary response) matrix of coefficients, including the intercept.  }
\references{
Taddy (2012), \emph{Multinomial Inverse Regression for Text Analysis}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  mnlm
}
