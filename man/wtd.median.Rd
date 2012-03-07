%  File ergm/man/wtd.median.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{wtd.median}
\alias{wtd.median}
\title{Weighted Median}
\description{
  Compute weighted median.
  }
\usage{
wtd.median (x, na.rm = FALSE, weight=FALSE)
} 
\arguments{
  \item{x}{ Vector of data, same length as \code{weight} }
  \item{na.rm}{Logical:  Should NAs be stripped before computation proceeds?}
  \item{weight}{ Vector of weights }
}
\details{
 Uses a simple algorithm based on sorting.
 }
\value{
  Returns an empirical .5 quantile from a weighted sample.
}
\keyword{robust}
