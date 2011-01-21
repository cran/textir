\name{freq}
\alias{freq}
\title{
Frequency Matrix Conversion
}
\description{
Convert a count matrix with observation-rows and category-columns to the corresponding frequency matrix.}
\usage{
freq(x)
}
\arguments{
  \item{x}{
    A matrix with count entries.  For use with the \code{tm} package, \code{freq} will first transpose a \code{TermDocumentMatrix}.
}
}
\value{
 A matrix with row sums of one.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
freq( t(rmultinom(10, 20, c(1/2,1/4,1/8,1/8))) )
 }
\seealso{ normalize }
