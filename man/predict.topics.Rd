\name{predict.topics}
\alias{predict.topics}
\title{
topic predict
}
\description{ Predict function for Topic Models }
\usage{
\method{predict}{topics}( object, newcounts, grp=NULL, ... )
}
\arguments{
  \item{object}{An output object from the \code{topics} function.}
  \item{newcounts}{  An \code{nrow(object$theta)}-column matrix of multinomial phrase/category counts for new documents/observations. 
Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}. }
  \item{grp}{ If \code{!object$admix}, this is an optional group membership vector. }
  \item{...}{Additional arguments to the undocumented internal \code{tpx*} functions. }
}
\details{ Under the default mixed-membership topic model, this function uses sequential quadratic programming to fit topic weights \eqn{\Omega} for new documents.   
Estimates for each new \eqn{\omega_i} are, conditional on \code{object$theta}, 
MAP in the (K-1)-dimensional logit transformed parameter space. If \code{!object$admix}, this instead returns the highest posterior probability topic for each document. }
\value{ The output is an \code{nrow(newcounts)} by \code{object$K} matrix of document topic weights, or a membership vector if  \code{!object$admix}. }
\references{
Taddy (2011), \emph{Estimation of Topic Models}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  topics, plot.topics, summary.topics, we8there, congress109, wsjibm
}
