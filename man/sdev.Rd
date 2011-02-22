\name{sdev}
\alias{sdev}
\title{Sparse Matrix Standard Deviation}
\description{ Standard deviation for columns of a \code{simple_triplet_matrix}.}
\usage{ sdev(x) }
\arguments{
  \item{x}{A \code{simple_triplet_matrix} (or a \code{matrix}, in which case the 
  function returns \code{apply(x,2,sd)}). }
}
\value{
 An \code{ncol(x)}-length vector containing standard deviations of the columns of x.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
data(congress109)
sdev(congress109Counts)[1:20]
 }
\seealso{ sd, freq, corr, congress109 }
