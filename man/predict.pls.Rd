\name{predict.pls}
\alias{predict.pls}
\title{
pls predict
}
\description{ Predict function for Partial Least Squares }
\usage{
\method{predict}{pls}( object, newdata, response=TRUE, ... )
}
\arguments{
  \item{object}{An output object from the \code{pls} function.}
  \item{newdata}{    An \code{nrow(object$loadings)}-column matrix of 
  		     multinomial phrase/category counts for new documents/observations. 
		     Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.   }
  \item{response}{  An indicator for whether to scale predictions to the original response. }
  \item{...}{Additional unused arguments. }
}
\details{ This function returns the pls projection \eqn{X\Phi} for new covariates, or \eqn{\alpha + \beta*X\phi} if \code{response=TRUE} with regression coefficients taken from \code{object$fwdmod}.
	 }
\value{ Output is an \code{nrow(newcounts)} by \code{ncol(object$loadings)} matrix of pls directions for each new observation. }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  pls, congress109
}
