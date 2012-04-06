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
## fit a simple topic model
summary( newstpx <- topics(wsjibmCounts, K=10, tol=100), nwrd=10 )

\dontrun{
## fit topics over years, using prior shape to allow them to change in time
year <- factor(1900 + as.POSIXlt(wsjibmReturns$DATE)$year)
Y <- nlevels(year)
annualtopics <- vector(length=Y, mode="list")
topwords <- c()
shape=NULL
for(i in 1:Y){ 
      annualtopics[[i]] <- topics(wsjibmCounts[year==levels(year)[i],], K=5, shape=shape, ord=FALSE)
      topwords <- cbind(topwords, as.character(summary(annualtopics[[i]], verb=FALSE)$phrase))
      delta <- 10000 # weight of the previous year in number of words observed per topic		
      shape <- annualtopics[[i]]$theta*delta }
## top 5 words by topic in past 4 years
dimnames(topwords) <- list(topic=rep(1:5,each=5), year=levels(year))
print(topwords[,Y - 3:0]) }	  	         
}