\name{predict.pls}
\alias{predict.pls}
\title{
pls predict
}
\description{ Predict function for Partial Least Squares }
\usage{
\method{predict}{pls}( object, newdata, type="response", ... )
}
\arguments{
  \item{object}{An output object from the \code{pls} function.}
  \item{newdata}{    An \code{nrow(object$loadings)}-column matrix of 
  		     multinomial phrase/category counts for new documents/observations. 
		     Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.   }
  \item{type}{  If "response", predictions scaled to the original response.  If "reduction", fitted partial least squares directions. }
  \item{...}{Additional unused arguments. }
}
\details{ This function returns the pls projection \eqn{X\Phi} for new covariates, or \eqn{\alpha + \beta*X\phi} if \code{type="response"} with regression coefficients taken from \code{object$fwdmod}.
	 }
\value{ Output is either a vector of predicted resonse or an \code{nrow(newcounts)} by \code{ncol(object$loadings)} matrix of pls directions for each new observation. }
\references{
Taddy (2012), \emph{Multinomial Inverse Regression for Text Analysis}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  pls, congress109
}
