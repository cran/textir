\name{polynomial roots}
\alias{cubic}
\alias{quadratic}
\title{
Cubic and quadratic function solvers
}
\description{
Find analytical roots to qubic and quadratic polynomials. }
\usage{
quadratic(b, c, quiet=FALSE, plot=FALSE)
cubic(a, b, c, quiet=FALSE, plot=FALSE)
}
\arguments{
  \item{a,b,c}{
    Polynomial function coefficients (MONIC FORM).
}
  \item{quiet}{
  If false, the solution is printed to screen.   
}
 \item{plot}{
  If true, the function and real root(s) are plotted.
}
}
\details{
	Finds roots to the cubic function \eqn{y = x^3 + ax^2 + bx + c} or quadratic function \eqn{y = x^2 + bx + c}.}
\value{
 A list with entries for the coefficients, roots, and solution characterization.  In particular, 
\item{type}{The solution characterization: number of complex and real roots.}
\item{coef}{The input coefficients.}
\item{roots}{A vector of the equation roots.    

For the quadratic equation, if there are complex roots, \code{roots[1]} is the real part and \code{roots[2]} 
is the imaginary part (i.e., complex roots are \code{roots[1] +- roots[2]*i}).  
Otherwise, \code{roots} are (possibly identical) real roots.

For the cubic equation, the first root \code{roots[1]} is always real.  
If there are complex roots, \code{roots[2]} is the real part and \code{roots[3]} 
is the imaginary part (i.e., complex roots are \code{roots[2] +- roots[3]*i}).  
Otherwise, \code{roots[2:3]} are (possibly identical) real roots.}
}
\references{
   Abramovitz and Stegun, Handbook of Mathematical Functions, 1972.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}


\seealso{
  'polyroot' for numerical solutions.
}
\examples{
quadratic(1,-2, plot=TRUE)
cubic(0,-15,-4, plot=TRUE)
 }
