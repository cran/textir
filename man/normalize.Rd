\name{normalize}
\alias{normalize}
\title{
Normalize
}
\description{
Normalize matrix columns. }
\usage{
normalize(x, m=NULL, s=NULL)
}
\arguments{
  \item{x}{
    A matrix.
}
\item{m}{Optional column shifts.}
\item{s}{Optional column scalings.}
}
\value{
 Under default, a matrix with mean-zero and variance-one columns.  
 If shift and scale are specified, a matrix with columns shifted by \code{-m}  and divided by \code{s}. 
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
normalize( matrix(1:9, ncol=3) )
 }
\seealso{ freq, pls }
