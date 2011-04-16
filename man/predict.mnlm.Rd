\name{predict.mnlm}
\alias{predict.mnlm}
\title{
mnlm predict
}
\description{ Predict function for Multinomial Logistic Regression }
\usage{
\method{predict}{mnlm}( object, newdata, type=c("response","reduction"), ... )
}
\arguments{
  \item{object}{An output object from the \code{mnlm} function.}
  \item{type}{ Under "reduction", provide the fitted reduction \eqn{F\phi}.  Under "response", provide the fitted multinomial probabilities.}
  \item{newdata}{   Under \code{"response"}, an \code{ncol(object$loadings)}-column matrix of new covariates. 
		     Under \code{"reduction"}, an \code{nrow(object$loadings)}-column matrix of 
  		     multinomial phrase/category counts for new documents/observations. 
		     Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.   }
  \item{...}{Additional unused arguments. }
}
\details{ This function returns either the inverse regression sufficient reduction \eqn{F\Phi} for new documents, or
	  fitted multinomial probabilities given new covariate vectors.
	 }
\value{ Under \code{type="response"}, output is an  \code{nrow(newcounts)} by \code{nrow(object$loadings)} matrix of 
	predicted probabilities for each response category. 
	Under \code{type="reduction"}, output is an \code{nrow(newcounts)} by \code{ncol(object$loadings)} matrix of document 
	scores in each factor (\code{object$covars}) direction. }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  mnlm, congress109
}
