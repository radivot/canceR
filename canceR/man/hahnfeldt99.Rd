\name{hahnfeldt99}
\alias{hahnfeldt99}
\title{Antiangiogenesis model of Hahnfeldt et al. 1999}
\usage{
  hahnfeldt99(Time, State, Pars)
}
\arguments{
  \item{Time}{The parameters of this model do not depend on
  time: this argument exists here as a dummy only because
  deSolve expects it to exist.}

  \item{State}{Vector of current states. The elements of
  this vector must be named because within this function's
  definition is \code{with(as.list(State,Pars), {XXX} )}
  and code in XXX presupposes that it can call named
  elements as variables with those names.}

  \item{Pars}{Vector of parameter values with \emph{named}
  (see above) elements.}
}
\value{
  A list of length 2 (as expected by deSolve) where the
  first list element is the vector of derivatives, i.e. the
  right hand side of the ODEs and the second element of the
  list is a vector of auxiliary variables that one may wish
  to track over time.
}
\description{
  This function returns the right hand side of an ordinary
  differential equation model of antiangiogenesis and tumor
  growth that was published by Hahnfeldt et al in Cancer
  Research in 1999.  The intended use of this function is
  as an argument to \code{ode()} of the \code{deSolve}
  package.
}
\details{
  The model captures tumor carrying capacity dynamics.
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(canceR)
out1=lsoda(y=c(V=180,K= 625,g=0),times=seq(0,20,1),hahnfeldt99,
	parms=c(lambda1=.192,lambda2=0,b=5.85,d=0.00873,e=0,clr=0),rtol=1e-4, atol= rep(1e-4,3))
plot(out1)
}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  Philip Hahnfeldt, Dipak Panigrahy, Judah Folkman, and
  Lynn Hlatky, Tumor Development under Angiogenic
  Signaling: A Dynamical Theory of Tumor Growth, Treatment
  Response, and Postvascular Dormancy, \emph{Cancer
  Research}, \bold{59}, 4770-4775 (1999)
}
\seealso{
  \code{\link{canceR-package}}
}
\keyword{IO}

