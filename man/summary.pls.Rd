\name{summary.pls}
\alias{summary.pls}
\alias{print.pls}
\title{
pls summary
}
\description{ Summary function for Partial Least Squares }
\usage{
\method{summary}{pls}( object, ... )
}
\arguments{
  \item{object}{An output object from the \code{pls} function.}
  \item{...}{Additional unused arguments. }
}
\details{ A short summary function for pls objects.
	 }
\value{ A printout of the number of pls directions and the input dimension, followed by a summary of the 
	corresponding forward regression \code{lm(as.numeric(y)~directions)}.  }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  pls
}
