\name{tfidf}
\alias{tfidf}
\title{
Term Frequency * Inverse Document Frequency
}
\description{
Convert a count matrix to the corresponding tfidf matrix.}
\usage{
tfidf(x, freq=FALSE)
}
\arguments{
  \item{x}{A matrix or \code{simple_triplet_matrix}. }
  \item{freq}{An indicator for whether \code{x} is already a frequency matrix. }
}
\value{
 A matrix with entries \eqn{f_{ij} \log(n/d_j)}, where \eqn{f_{ij}} is term-j frequency in document-i, 
and \eqn{d_j} is the number of documents containing term-j.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
## 20 important terms
data(congress109)
sort(sdev(tfidf(congress109Counts)), decreasing=TRUE)[1:20]
 }
\seealso{ freq }
