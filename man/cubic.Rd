\name{cubic}
\alias{cubic}
\title{
Cubic function solver
}
\description{
Finds analytical roots to the function y = x^3 + ax^2 + bx + c }
\usage{
cubic(a, b, c, quiet=FALSE, plot=FALSE)
}
\arguments{
  \item{a,b,c}{
    Cubic function coefficients (MONIC FORM).
}
  \item{quiet}{
  If false, the solution is printed to screen.   
}
 \item{plot}{
  If true, the function and real root(s) are plotted.
}
}
\details{
	Finds roots to the cubic function \eqn{y = x^3 + ax^2 + bx + c}.}
\value{
 A list with entries for the coefficients, roots, and solution characterization.  In particular, 
\item{type}{The solution characterization is either "one real", "three real", or "one real, two complex".}
\item{coef}{The input coefficients \code{c(a,b,c)}.}
\item{roots}{A vector for the cubic roots.  The first root \code{roots[1]} is always real.  
If there are complex roots, \code{roots[2]} is the real part and \code{roots[3]} 
is the imaginary part (i.e., complex roots are \code{roots[2] +- roots[3]*i}).  
Otherwise, \code{roots[2:3]} are (possibly identical) real roots.}
}
\references{
   Gerolamo Cardano, Ars Magna (1545)
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}


\seealso{
  'polyroot' for numerical solutions.
}
\examples{
cubic(0,-15,-4, plot=TRUE)
 }
