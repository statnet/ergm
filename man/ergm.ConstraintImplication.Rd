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
