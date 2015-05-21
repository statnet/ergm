\name{ergm.MCMLE}
\alias{ergm.MCMLE}
\title{
Estimation Techniques for Finding the MLE
}
\description{
 The \code{ergm.MCMLE} function provides one of the styles of maximum
 likelihood estimation that can be used. This one is the default and uses
 optimization of an MCMC estimate of the log-likelihood.  Most of the parameters are documented on the \link{control.ergm} page. 
}
\usage{
ergm.MCMLE(init, nw, model, initialfit, control, MHproposal, 
    MHproposal.obs, verbose = FALSE, sequential = control$MCMLE.sequential, 
    estimate = TRUE, response = NULL, ...)
}

\arguments{
  \item{init}{
the initial theta values
}
  \item{nw}{
the network 
}
  \item{model}{
the model, as returned by \code{\link{ergm.getmodel}}
}
  \item{initialfit}{
an ergm object, as the initial fit, possibly returned
                   by \code{\link{ergm.initialfit}}
}
  \item{control}{
a list of \code{\link{control.ergm}} parameters for controlling the MCMC sampling; most parameters have prefix "MCMLE"
}
  \item{MHproposal}{
an MHproposal object for 'nw', as returned by  \code{\link{MHproposal}}
}
  \item{MHproposal.obs}{
an MHproposal object for the observed network of'nw', as returned by \code{\link{MHproposal}}
}
  \item{verbose}{
logical, whether the MCMC sampling should be verbose. default=FALSE
}
  \item{sequential}{
whether to update the network returned in
                    'v$newnetwork'; if the network has missing edges,
                   this is ignored; default=control$MCMLE.sequential
}
  \item{estimate}{
whether to optimize the init coefficients (via ergm.estimate); default=TRUE
}
  \item{response}{
name of an edge attribute
}
  \item{\dots}{
additional parameters that may be passed from within; all are ignored
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 returns an ergm object as a list containing several items; for details see
      the return list for \code{\link{ergm}}
      note that if the model is degenerate, only 'coef' and 'sample' are
      returned; if 'estimate'=FALSE, the MCMC and se variables will be
      NA or NULL
}


