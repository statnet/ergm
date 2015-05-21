\name{ergm.initialfit}
\alias{ergm.initialfit}
\title{
Internal function to fit an initial ergm object using either ML or MPL
estimation
}
\description{
The \code{ergm.initialfit} function fits an initial ergm object using either ML or MPL
estimation.  If initial parameters are provided in 'init' and 'MLestimate' is 
TRUE, the number of parameters in 'init' is checked for correctness.
}


\usage{
ergm.initialfit(init,
                initial.is.final, 
                formula, 
                nw, 
                m, 
                response = NULL, 
                reference = ~Bernoulli, 
                method = NULL, 
                MPLEtype = "glm", 
                conddeg = NULL, 
                control = NULL, 
                MHproposal = NULL, 
                MHproposal.obs = NULL, 
                verbose = FALSE, 
                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
either a vector whose first entry is "MPLE" or a vector
                    of initial coefficients
}
  \item{initial.is.final}{
%%     ~~Describe \code{initial.is.final} here~~
}
  \item{formula}{
a formula of the form (nw ~ term(s)) 
}
  \item{nw}{
a \code{\link{network}} object, presumably that of 'formula'
}
  \item{m}{
the model as returned by \code{\link{ergm.getmodel}}
}
  \item{response}{
optional character string giving the name of an edge attribute (only used with CD)
}
  \item{reference}{
not used?
}
  \item{method}{
one of "MPLE", "zeros", "CD"
}
  \item{MPLEtype}{
the method for MPL estimation as either "glm", "penalized",
                    or "logitreg"; this is ignored if ML estimation is used;
                    default="glm". Passed to \code{\link{ergm.mple}}
}
  \item{conddeg}{
a formula for the conditional degree terms
}
  \item{control}{
a list of parameters for tuning the MCMC sampling; (see \code{\link{control.ergm}})
                   the only recognized component is 'samplesize'
}
  \item{MHproposal}{
an MHproposal object as returned by \code{\link{MHproposal}}
}
  \item{MHproposal.obs}{
%%     ~~Describe \code{MHproposal.obs} here~~
}
  \item{verbose}{
whether the MPL estimation should be verbose (T or F); 
                   default=FALSE
}
  \item{\dots}{
addtional parameters that are used with MPL estimation;
                    the only recognized parameeter is 'compressflag' which
                   compresses the design matrix used by \code{\link{ergm.mple}} 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns  an \code{\link{ergm}} object as one of the following lists
\itemize{
  \item{ if MLE: a list with 2 components
   \itemize{
                \item coef   : 'init'
                \item mle.lik:  the MLE likelihood}}
  \item { if MPLE -- the list returned by \code{\link{ergm.mple}}
}
}



\keyword{internal}
