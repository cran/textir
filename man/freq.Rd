\name{freq}
\alias{freq}
\title{
Frequency Matrix Conversion
}
\description{
Convert a count matrix to the corresponding frequency matrix.}
\usage{
freq(x, byrow=TRUE)
}
\arguments{
  \item{x}{A matrix or \code{simple_triplet_matrix} with count entries. }
  \item{byrow}{An indicator for whether you have observation-rows and category-columns, or vice versa. }
}
\value{
 A matrix with row (\code{byrow=TRUE}) or column (\code{byrow=FALSE}) sums of one.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
freq( t(rmultinom(10, 20, c(1/2,1/4,1/8,1/8))) )
 }
\seealso{ corr, pls }
