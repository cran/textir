\name{predict.topics}
\alias{predict.topics}
\title{
topic predict
}
\description{ Predict function for Topic Models }
\usage{
\method{predict}{topics}( object, newcounts, ... )
}
\arguments{
  \item{object}{An output object from the \code{topics} function.}
  \item{newcounts}{  An \code{nrow(object$theta)}-column matrix of multinomial phrase/category counts for new documents/observations. 
Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}. }
  \item{...}{Additional arguments to the undocumented internal \code{tpx*} functions. }
}
\details{ This function uses sequential quadratic programming to fit topic weights \eqn{\Omega} for new documents.   
Estimates for each new \eqn{\omega_i} are, conditional on \code{object$theta}, 
MAP in the (K-1)-dimensional logit transformed parameter space.}
\value{ The output is an \code{nrow(newcounts)} by \code{object$K} matrix of document topic weights. }
\references{
Taddy (2011), \emph{Estimation of Topic Models}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  topics, plot.topics, summary.topics, we8there, congress109, wsjibm
}
