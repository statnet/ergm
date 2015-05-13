%  File man/ergm.allstats.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.init.methods}
\alias{ergm.init.methods}
\title{Set up the initial fitting methods for reference measure}
\description{
  This is a low-level function not intended to be called directly by end users. This function sets up the available initial fitting methods 
    for each reference measure. (It is defined in the scope of local environment)  
}
\usage{
 ergm.init.methods(reference, new.methods)
	}
\arguments{
  \item{reference}{The reference measure used in the model.}
  \item{new.methods}{Add new initial fitting methods to the reference measure.}
}
