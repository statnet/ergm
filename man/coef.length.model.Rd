%  File man/coef.length.model.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{coef.length.model}
\alias{coef.length.model}
\alias{coef.sublength.model}
\alias{eta.length.model}
\alias{eta.sublength.model}
\title{Extract Number of parameters in a \code{ergm} Model}
\usage{

\method{coef.length}{model}(object, offset=NA, ...)

\method{coef.sublength}{model}(object, offset=NA, ...)

\method{eta.length}{model}(object, offset=NA, ...)

\method{eta.sublength}{model}(object, offset=NA, ...)

}
\arguments{
  \item{object}{an ergm model object}
  \item{offset}{If \code{NA} (the default), all model terms are
    counted; if \code{TRUE}, only offset terms are counted; and if
    \code{FALSE}, offset terms are skipped.}
  \item{\dots}{other arguments.}
}
\value{
  \code{coef.sublength.model} returns a vector containing the number of
  model parameters corresponding to each model term.  
  \code{coef.length.model} returns the total number of parameters.

  \code{eta.sublength.model} returns a vector containing the number of
  canonical parameters (also the number of sufficient statistics) corresponding to each model term.  
  \code{eta.length.model} returns the total number of canonical parameters.

}
\description{
  \code{coef.sublength}, \code{coef.length}, \code{eta.sublength}, and \code{eta.length} are methods that extract
the numbers of parameters for ergm model objects.
}
\keyword{models}
