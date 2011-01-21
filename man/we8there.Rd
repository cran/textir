\name{we8there}
\alias{we8there}
\alias{we8thereRatings}
\alias{we8thereCounts}

\title{ On-Line Restaurant Reviews }

\description{
Counts for 2978 bigrams in 6147 restaurant reviews from the site www.we8there.com. }

\details{ The short user-submitted reviews are accompanied by a five-star rating on four specific aspects of restaurant quality - food, service, value, and atmosphere - as well as the overall experience.  The reviews originally appear in Maua and Cozman (2009), and the parsing details behind these specific counts are in Taddy (2011). }

\value{
   \item{we8thereCounts}{ A \code{simple_triplet_matrix} of phrase counts indexed by review-rows and bigram-columns.}
   \item{we8thereRatings}{ A \code{matrix} containing the associated review ratings. }
}

\references{ 
Maua, D.D. and Cozman, F.G. (2009), \emph{Representing and classifying user reviews}. In ENIA '09: VIII Enconro Nacional de Inteligencia Artificial, Brazil.

Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}

\author{ 
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\seealso{ pls, mnlm, congress109 }

\examples{
data(we8there)
y <- normalize(we8thereRatings$Overall)
F <- normalize(freq(we8thereCounts))
fit <- pls(F, y)
plot(fit, y=as.factor(we8thereRatings$Overall), col=c(2,2,2,3,3))

## Uncomment the following to get multinomial inverse regression.  
## The MAP estimation can take up to 5 minutes, depending on your machine.
# fit <- mnlm(we8thereCounts, y)
# par(mfrow=c(1,2))
# plot(fit, y= as.factor(we8thereRatings$Overall), col=c(2,2,2,3,3))
# plot(fit, cat=2842, delta=0.05)


}