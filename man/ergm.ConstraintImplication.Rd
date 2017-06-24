%  File man/ergm.ConstraintImplication.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.ConstraintImplications}
\alias{ergm.ConstraintImplications}
\title{Set up the implied constraints from the current constraint}
\description{
  This is a low-level function not intended to be called directly by end users. For information on constraints, see the \code{\link{ergm-constraints}} page.  This function set up the implied constraints from the current
  constraint. (It is defined in the scope of local environment)  
}
\usage{
 ergm.ConstraintImplications(implier, implies)
	}
\arguments{
  \item{implier}{The current constraint specified in the model. 
  For the list of constraints, see \code{\link{ergm-constraints}}}
  \item{implies}{Implied constraints from the current constraint 
  (based on the user's knowledge).}
}
