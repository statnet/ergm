%  File man/ergmMPLE.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergmMPLE}
\title{ERGM Predictors and response for logistic regression calculation of MPLE}
\usage{
ergmMPLE(formula, fitmodel=FALSE, output=c("matrix","array", "fit"),
         as.initialfit = TRUE, control=control.ergm(),
         verbose=FALSE, \dots) 
}
\alias{ergmMPLE}
\description{
  Return the predictor matrix, response vector, and vector of weights
  that can be used to calculate the MPLE for an ERGM.
}
\arguments{
  \item{formula}{An ERGM formula. See \code{\link{ergm}}.}
  \item{fitmodel}{Deprecated. Use \code{output="fit"} instead.}
  \item{output}{Character, partially matched. See Value. }
  \item{as.initialfit}{Logical. Specifies whether terms are initialized
    with argument \code{initialfit==TRUE} (the default). Generally, if
    \code{TRUE}, all curved ERGM terms will be treated as having their
    curved parameters fixed. See Example. }
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

  Using \code{output="matrix"}, note that the result of the fit may be
  obtained from the \code{\link{glm}} function, as shown in the examples
  below.
  
  When \code{output="array"}, the \code{MPLE.max.dyad.types} control
  parameter must be greater than \code{network.dyadcount(.)} of the
  response network, or not all elements of the array that ought to be
  filled in will be.  
}
\value{
  
  If \code{output=="matrix"} (the default), then only the response,
  predictor, and weights are returned; thus, the MPLE may be found by
  hand or the vector of change statistics may be used in some other
  way. To save space, the algorithm will automatically search for any
  duplicated rows in the predictor matrix (and corresponding response
  values). \code{ergmMPLE} function will return a list with three
  elements, \code{response}, \code{predictor}, and \code{weights},
  respectively the response vector, the predictor matrix, and a vector
  of weights, which are really counts that tell how many times each
  corresponding response, predictor pair is repeated.

  If \code{output=="array"}, a list with similarly named three elements
  is returned, but \code{response} is formatted into a sociomatrix;
  \code{predictor} is a 3-dimensional array of with cell
  \code{predictor[t,h,k]} containing the change score of term
  \code{k} for dyad (\code{t},\code{h}); and \code{weights} is also
  formatted into a sociomatrix, with an element being 1 if it is to be
  added into the pseudolikelihood and 0 if it is not.

  In particular, for a unipartite network, cells corresponding to
  self-loops, i.e., \code{predictor[i,i,k]} will be \code{NA} and
  \code{weights[i,i]} will be 0; and for a unipartite undirected
  network, lower triangle of each \code{predictor[,,k]} matrix will be
  set to \code{NA}, with the lower triangle of \code{weights} being set
  to 0.

  If \code{output=="fit"}, then \code{ergmMPLE} simply calls the
  \code{\link{ergm}} function with the \code{estimate="MPLE"} option
  set, returning an object of class \code{ergm} that gives the fitted
  pseudolikelihood model.
  
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
ergmMPLE(formula, output="fit")$coef

# We can also format the predictor matrix into an array:
mplearray <- ergmMPLE(formula, output="array")

# The resulting matrices are big, so only print the first 5 actors:
mplearray$response[1:5,1:5]
mplearray$predictor[1:5,1:5,]
mplearray$weights[1:5,1:5]

formula2 <- faux.mesa.high ~ gwesp(0.5,fix=FALSE)

# The term is treated as fixed: only the gwesp term is returned:
colnames(ergmMPLE(formula2, as.initialfit=TRUE)$predictor)

# The term is treated as curved: individual esp# terms are returned:
colnames(ergmMPLE(formula2, as.initialfit=FALSE)$predictor)
}
\keyword{regression}
\keyword{models}
