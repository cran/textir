\name{plot.pls}
\alias{plot.pls}
\title{
pls plot
}
\description{ Plot function for Partial Least Squares }
\usage{
\method{plot}{pls}(x, K=NULL, ...)
}
\arguments{
  \item{x}{An output object from the \code{pls} function.}
\item{K}{The number of pls directions to be used. Can be a vector.
 If \code{K}, plot fitted values for \code{1:fit$K} directions.}
 \item{...}{Additional \code{plot} arguments}
}

\details{
Plots response versus fitted values for least-squares fit onto the K pls directions. }

\value{ A fabulous plot. }
\references{
   
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}

}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  pls, we8there
}
