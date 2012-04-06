\name{normalize}
\alias{normalize}
\title{
Normalize
}
\description{
Normalize matrix columns. }
\usage{
normalize(x, m=NULL, s=NULL, undo=FALSE)
}
\arguments{
  \item{x}{
    A matrix.
}
\item{m}{Optional column shifts.}
\item{s}{Optional column scalings.}
\item{undo}{If \code{undo=TRUE} this will undo a previous normalization.  Otherwise, just normalize.}
}
\value{
 Under default, a matrix with mean-zero and variance-one columns.  
 If shift and scale are specified, a matrix with columns shifted by \code{-m}  and divided by \code{s}.
 If \code{undo=TRUE}, and shift and scale are specified, an un-normalized matrix with column means \code{m} and standard deviations \code{s}.
 In the special case where \code{x} is a \code{simple_triplet_matrix} and \code{m=0}, columns are scaled by \code{s} but left unshifted and the function returns a \code{simple_triplet_matrix}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
normalize( matrix(1:9, ncol=3) )
 }
\seealso{ freq, corr, sdev, pls  }
