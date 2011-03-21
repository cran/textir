\name{congress109}
\alias{congress109}
\alias{congress109Ideology}
\alias{congress109Counts}

\title{ Ideology in Political Speeches }

\description{ Phrase counts and ideology scores by speaker for members of the 109th US congress. }

\details{ This data originally appears in Gentzkow and Shapiro
(GS; 2010) and considers text of the 2005 Congressional Record,
containing all speeches in that year for members of the United States
House and Senate. In particular, GS record the number times each of
529 legislators used terms in a list of 1000 phrases (i.e., each
document is a year of transcripts for a single speaker). Associated
sentiments are repshare -- the two-party vote-share from each
speaker's constituency (congressional district for representatives;
state for senators) obtained by George W. Bush in the 2004
presidential election -- and the speaker's first and second
common-score values (from http://voteview.com). Full
parsing and sentiment details are in Taddy (2011; Section 2.1). }

\value{
   \item{congress109Counts}{ A \code{simple_triplet_matrix} of phrase counts indexed by speaker-rows and phrase-columns.}
   \item{congress109Ideology}{ A \code{matrix} containing the associated \code{repshare} and common scores \code{[cs1,cs2]}, as well as speaker
   characteristics: \code{party} (`R'epublican, `D'emocrat, or `I'ndependent), \code{state}, and \code{chamber} (`H'ouse or `S'enate). }
}

\references{ 
Gentzkow, M. and J. Shapiro (2010), \emph{What drives media slant? Evidence from U.S. daily newspapers}. Econometrica 78, 35-7.  The full dataset is at \url{http://dx.doi.org/10.3886/ICPSR26242}.

Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{ 
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\seealso{ mnlm, pls, we8there, plot.mnlm, normalize }

\examples{
data(congress109)

## Inverse Regression Sentiment Modeling 
fit <- mnlm(congress109Counts, congress109Ideology[,6:7], normalize=TRUE, tol=1)
par(mfrow=c(1,2))
plot(fit, v=congress109Ideology$repshare, xlab="Republican Vote-Share",
	  covar=1, pch=21, bg=c(4,3,2)[congress109Ideology$party], main="1st common score")
plot(fit, v=congress109Ideology$repshare, xlab="Republican Vote-Share", 
	  covar=2, pch=21, bg=c(4,3,2)[congress109Ideology$party], main="2nd common score")
## example usage of the predict method
predict(fit, newcounts=congress109Counts[c(68,388),])

## A small topic model 
par(mfrow=c(1,1))
tpx <- topics(congress109Counts, K=15)
plot(tpx, group=congress109Ideology$party=="R", col=c(4,2), labels=c("Dem","GOP"))
summary(tpx)

}