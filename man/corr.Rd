\name{corr}
\alias{corr}
\title{Sparse Matrix Correlation}
\description{
Correlation calculation for a \code{simple_triplet_matrix} and a \code{matrix}.}
\usage{
corr(x, y)
}
\arguments{
  \item{x}{A \code{simple_triplet_matrix} (or a \code{matrix}, in which case the 
  function returns \code{cor(x,y)}). }
  \item{y}{A \code{matrix} with \code{nrow(y)=nrow(x)}. }
}
\value{
 An \code{ncol(x)} by \code{ncol(y)} matrix containing correlation between x and y.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
data(congress109)
r <- corr(congress109Counts, congress109Ideology$repshare)
## 20 terms for Democrats
sort(r[,1])[1:20]
## 20 terms for Republicans
sort(r[,1], decreasing=TRUE)[1:20]

 }
\seealso{ cor, sdev, freq, congress109 }
