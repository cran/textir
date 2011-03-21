\name{wsjibm}
\alias{wsjibm}
\alias{wsjibmReturns}
\alias{wsjibmCounts}

\title{ WSJ Stories on IBM }

\description{ Word counts for Wall Street Journal story abstracts with IBM in the title,
 along with the concurrent returns on IBM stock. }

\details{  Headlines and one-sentence abstracts for Wall Street Journal (WSJ) stories with IBM in the headline, dating from August 1988 to August 2010, were retrieved from the ProQuest database. Each article is accompanied by two-day return and return-over-market for shares in IBM listed on the New York Stock Exchange, calculated from the opening of the previous day to market close on the day of publication.  Full details are available in Taddy (2011). }

\value{
   \item{wsjibmCounts}{ A \code{simple_triplet_matrix} of counts indexed by article-rows and word-columns.}
   \item{wsjibmReturns}{ A \code{matrix} containing the corresponding publication \code{DATE} along with IBM's two-day holding returns (\code{RET}) and return over the S&P500 (\code{ROM}). }
}

\references{ 
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{ 
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\seealso{ topics, plot.topics }

\examples{
data(wsjibm)
## fit a small topic model
newstpx <- topics(wsjibmCounts, K=5, tol=1) 
plot(newstpx, group=wsjibmReturns$ROM>=0, labels=c("loss","gain"), col=c(2,1))
summary(newstpx, nwrd=10)
}