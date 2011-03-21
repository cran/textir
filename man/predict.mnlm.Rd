\name{predict.mnlm}
\alias{predict.mnlm}
\title{
topic predict
}
\description{ Predict function for Logistic Inverse Regression }
\usage{
\method{predict}{mnlm}( object, newcounts, ... )
}
\arguments{
  \item{object}{An output object from the \code{mnlm} function.}
  \item{newcounts}{  An \code{nrow(object$loadings)}-column matrix of multinomial phrase/category counts for new documents/observations. 
Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}. }
  \item{...}{Additional arguments to the undocumented internal \code{tpx*} functions. }
}
\details{ This function returns the sufficient reduction \eqn{F\Phi} for new documents, based on loadings from a fitted mnlm object.}
\value{ The output is an \code{nrow(newcounts)} by \code{ncol(object$loadings)} matrix of document scores in each factor (\code{object$covars}) direction. }
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
