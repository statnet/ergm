%  File man/ergm.init.methods.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.init.methods}
\alias{ergm.init.methods}
\title{Set up the initial fitting methods for reference measure and
  query available methods for that reference measure}
\description{
  This is a low-level function not intended to be called directly by end users. This function sets up the available initial fitting methods 
    for each reference measure and queries them.
}
\usage{
 ergm.init.methods(reference, new.methods)
	}
\arguments{
  \item{reference}{The reference measure used in the model.}
  \item{new.methods}{If passed, prepends the new initial fitting methods
    to the list for that reference measure.}
}
\value{
  A character vector listing initial methods for the reference measure
  specified. (If \code{new.methods} is passed, does so invisibly.)
}
