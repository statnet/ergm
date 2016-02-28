%  File man/coef.length.model.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{coef.length.model}
\alias{coef.length.model}
\alias{coef.sublength.model}
\title{Extract Number of parameters in ergm Model}
\usage{

\method{coef.length}{model}(object, ...)

\method{coef.sublength}{model}(object, ...)
}
\arguments{
  \item{object}{an ergm model object} 
  \item{\dots}{other arguments.}
}
\value{
  \code{coef.sublength.model} returns a vector containing the number of
  model parameters corresponding to each model term.  
  \code{coef.length.model} returns the sum of the values returned by
  \code{coef.sublength.model}.
}
\description{
  \code{coef.sublength} and \code{coef.length} are methods that extract
the numbers of parameters for ergm model objects.
}
\keyword{models}
