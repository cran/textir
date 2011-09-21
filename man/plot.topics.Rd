\name{plot.topics}
\alias{plot.topics}
\title{
topic plots
}
\description{ Plot function for Topic Models }
\usage{
\method{plot}{topics}(x, type=c("weight","resid"), group=NULL, labels=NULL, 
col=NULL, xlab=NULL, ylab=NULL, main=NULL, tpk=NULL, lgd.K=NULL, 
cex.lgdc = 1, cex.lgdt = 1, cex.rmar= 1, ... )
}
\arguments{
  \item{x}{An output object from the \code{topics} function.}
\item{type}{If "weight", the default, provide an image plot of document-topic weights. 
	    If "resid", just show a simple histogram of standardized residuals for the positive count entries. }
\item{group}{ Optional \code{logical} vector containing membership in some group for each document;
this will be used to color the topic-weight shadings.  
See the \code{textir} dataset examples, which color by good reviews for the we8there data or by republicans in congress109.}
\item{labels}{ Optional length-two \code{character} vector of labels for the membership specified in \code{groups}.
\code{labels[1]} corresponds to \code{group=FALSE} and \code{labels[2]} to \code{group=TRUE}. }
\item{col}{ If \code{type="weight"}, a number from \code{1:4} specifying the shade color (grey, followed by red, green blue).  
If \code{group} is specified, \code{col[1]} corresponds to \code{group=FALSE}, and \code{col[2]} to \code{group=TRUE}.
If \code{type="resid"}, this is just standard R coloring for the histogram bars.} 
\item{xlab}{ Optional x-axis label.}
\item{ylab}{ Optional y-axis label.}
\item{main}{ Optional title.}
\item{tpk}{ Optional list of topics to plot.  Defaults to \code{1:x$K}.  }
\item{lgd.K}{ Optional number of topic-increments (along the X-axis) 
outside of the plot region at which the legend is centered. }
\item{cex.lgdc}{ Magnification factor for legend color-boxes. }
\item{cex.lgdt}{ Magnification factor for legend text. }
\item{cex.rmar}{ Magnification factor for the right plot margin. }
\item{...}{ Additional arguments to the \code{image} function. }
}

\value{ A fabulous plot. }
\references{
Taddy (2011), \emph{Estimation and Selection for Topic Models}.
\url{http://arxiv.org/abs/1109.4518}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  topics, summary.topics, we8there, congress109, wsjibm
}
