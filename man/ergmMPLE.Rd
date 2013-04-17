%  File man/ergmMPLE.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergmMPLE}
\title{ERGM Predictors and response for logistic regression calculation of MPLE}
\usage{
ergmMPLE (formula, fitmodel=FALSE, control=control.ergm(),
          verbose=FALSE, \dots) 
}
\alias{ergmMPLE}
\description{
  Return the predictor matrix, response vector, and vector of weights
  that can be used to calculate the MPLE for an ERGM.
}
\arguments{
  \item{formula}{An ergm formula.  See \code{\link{ergm}}.}
  \item{fitmodel}{Logical:  Should the MPLE be found?  
      If \code{fitmodel==TRUE}, then \code{ergmMPLE} simply calls the 
      \code{\link{ergm}} function with the \code{estimate="MPLE"} option set.
      If \code{fitmodel==FALSE} (the default), then only the response, predictor,
      and weights are returned; thus, the MPLE may be found by hand or the vector
      of change statistics may be used in some other way.}
  \item{control}{A list of control parameters for tuning the fitting of an 
      ERGM.  Most of these parameters are irrelevant in this context.
      See \code{\link{control.ergm}} for details about all of 
      the control parameters.}
  \item{verbose}{Logical; if \code{TRUE}, the program will print out some
      additional information.}
  \item{\dots}{Additional arguments, to be passed to lower-level functions.}
  }
\details{
  The MPLE for an ERGM is calculated by first finding the matrix of change 
  statistics.  Each row of this matrix is associated with a particular pair 
  (ordered or unordered, depending on whether the network is directed
  or undirected) of nodes, and the row equals the change in the
  vector of network statistics (as defined in \code{formula}) when that
  pair is toggled from a 0 (no edge) to a 1 (edge), holding all the rest
  of the network fixed.  The MPLE results if we perform a logistic 
  regression in which the predictor matrix is the matrix of change statistics and
  the response vector is the observed network (i.e., each entry is either 0 or
  1, depending on whether the corresponding edge exists or not).
  
  To save space, the algorithm will automatically search for any duplicated 
  rows in the predictor matrix (and corresponding response values).  Thus,
  the \code{ergmMPLE} function will return three items:  A response
  vector, a predictor matrix, and a vector of weights, which are really counts
  that tell how many times each corresponding response, predictor pair is
  repeated.

  Using \code{fitmodel=FALSE},
  note that the result of the fit may be obtained from the 
  \code{\link{glm}} function, as shown in the examples below.  
}
\value{
  If \code{fitmodel==FALSE}, a list with three elements, called "response",
  "predictor", and "weights".  If \code{fitmodel==TRUE}, an object
  of class \code{ergm} that gives the fitted pseudolikelihood model.
}
\seealso{
  \code{\link{ergm}}, \code{\link{glm}}
}
\examples{
data(faux.mesa.high)
formula <- faux.mesa.high ~ edges + nodematch("Sex") + nodefactor("Grade")
mplesetup <- ergmMPLE(formula)

# Obtain MPLE coefficients "by hand":
glm(mplesetup$response ~ . - 1, data = data.frame(mplesetup$predictor), 
    weights = mplesetup$weights, family="binomial")$coefficients

# Check that the coefficients agree with the output of the ergm function:
ergmMPLE(formula, fitmodel=TRUE)$coef
}
\keyword{regression}
\keyword{models}
