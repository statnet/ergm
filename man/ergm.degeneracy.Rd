\name{ergm.degeneracy}
\alias{ergm.degeneracy}
\title{
Checks an ergm Object for Degeneracy
}
\description{
 The \code{ergm.degeneracy} function checks a given ergm object for degeneracy by
 computing and returning the instability value of the model and the value
 of the log-likelihood function at the maximized theta values
}
\usage{
ergm.degeneracy(object, control = object$control, 
                fast = TRUE, test.only = FALSE, 
                verbose = FALSE)
}
\arguments{
  \item{object}{
an \code{\link{ergm}} object
}
  \item{control}{
the list of control parameters as returned by \code{control.ergm};
               default=control.ergm()
}
  \item{fast}{
whether the degeneracy check should be "fast", i.e to sample
               changeobs(?) when there are > 100, rather than use all changeobs;
               default=TRUE
}
  \item{test.only}{
whether to silence printing of the model instability calculation
               (T or F); this parameter is ignored if the instability > 1;
               default=FALSE
}
  \item{verbose}{
whether to print a notification when 'object' is deemed degenerate
               (T or F); default=FALSE
}
}
\value{
   returns the original ergm object with 2 additional components:
   \itemize{
     \item degeneracy.value:  the instability of the model
     \item degeneracy.type :  a 2-element vector containing  
     \itemize{
     \item{loglikelihood: the value of the log-likelihood function corresponding to 'theta';
                    if degenerate, this is a vector of Inf}
     \item{theta        : the vector of theta values found through maximixing the log-
                   likelihood; if degenerate, this is 'guess'}
}
     }
}

