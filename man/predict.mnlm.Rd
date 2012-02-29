\name{predict.mnlm}
\alias{predict.mnlm}
\title{
mnlm predict
}
\description{ Predict function for Multinomial Logistic Regression }
\usage{
\method{predict}{mnlm}( object, newdata, type=c("response","reduction"), ... )
}
\arguments{
  \item{object}{An output object from the \code{mnlm} function.}
  \item{type}{ Under "reduction", provide the fitted reduction \eqn{F\phi}.  Under "response", provide the fitted multinomial probabilities.}
  \item{newdata}{   Under \code{"response"}, an \code{ncol(object$loadings)}-column matrix of new covariates. 
		     Under \code{"reduction"}, an \code{nrow(object$loadings)}-column matrix of 
  		     multinomial phrase/category counts for new documents/observations. 
		     Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.   }
  \item{...}{Additional unused arguments. }
}
\details{ Under \code{'response'}, this returns fitted multinomial probabilities given new covariate vectors.  Under \code{'reduction'}, we provide the sufficient reduction \eqn{F\Phi} for new documents, with \eqn{F} a document-term frequency matrix (i.e., the counts divided by document totals).
	 }
\value{ Under \code{type="response"}, output is an  \code{nrow(newcounts)} by \code{nrow(object$loadings)} matrix of 
	predicted probabilities for each response category. 
	Under \code{type="reduction"}, output is an \code{nrow(newcounts)} by \code{ncol(object$loadings)} matrix of document 
	scores in each factor (\code{object$covars}) direction. }
\references{
Taddy (2011), \emph{Inverse Regression for Analysis of Sentiment in Text}.
\url{http://arxiv.org/abs/1012.2098}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  mnlm, congress109
}

\examples{
## fit a congress 109 mnlm model using a random 300 members
data(congress109)
train <- sample(1:529, 300)
counts <- congress109Counts[,col_sums(congress109Counts[train,])>0]
fitRep <- mnlm(counts[train,], congress109Ideology$repshare[train], normalize=TRUE, bins=10)

## extract the reduced dimension text score
Z <- predict(fitRep, newdata=counts[train,], type="reduction")

## use this to build a forward regression model
fwdRep <- lm(repshare ~ Z, data=data.frame(repshare=congress109Ideology$repshare[train], Z=Z[,1]) )

## predict scores for the out-of-sample members
Znew <- predict(fitRep, newdata=counts[-train,], type="reduction")
predicted <- predict(fwdRep, newdata=data.frame(Z=Znew[,1]))
plot(congress109Ideology$repshare[-train], predicted, 
     pch=21, bg=c(4,3,2)[congress109Ideology$party[-train]], xlab="repshare")
abline(a=0,b=1, col=8) 
}