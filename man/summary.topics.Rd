\name{summary.topics}
\alias{summary.topics}
\title{
topic summary
}
\description{ Summary function for Topic Models }
\usage{
\method{summary}{topics}( object, nwrd=5, tpk=NULL, verb=TRUE, ... )
}
\arguments{
  \item{object}{ An output object from the \code{topics} function.}
  \item{nwrd}{ The number of phrases to output for each topic.}
  \item{tpk}{ Optional list of topics to summarize.  Defaults to \code{1:x$K}.  }
  \item{verb}{ Whether or not to print the summary.  }
  \item{...}{ Unused arguments from other functions, for S3 compatibility.}
}
\details{ This summary orders phrases for each topic according to the lift
\eqn{\theta_{kj}/q_{j}}, where \eqn{q_j} is the null-model probability estimate
\eqn{\sum_i x_{ij} / \sum_{i} m_i }.  This ordering of term relevance can be used to 
identify representative vocabulary for each topic.}
\value{ The function prints usage percentages (i.e. \code{colMeans(omega)}) 
and possible model selection results, along with the 
top \code{nwrd} phrases by term-lift for each topic in \code{tpk}.  The 
matrix of top \code{nwrd} phrases and their associated lift is returned invisibly. }
\references{
Taddy (2011), \emph{Estimation of Topic Models}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  topics, plot.topics, we8there, congress109, wsjibm
}
