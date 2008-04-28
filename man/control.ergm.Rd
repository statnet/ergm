\name{control.ergm}
\alias{control.ergm}
\title{ Auxiliary for Controlling ERGM Fitting }
\description{
  Auxiliary function as user interface for fine-tuning 'ergm' fitting.
}
\usage{
control.ergm(prop.weights = "default", prop.args = NULL,
             nr.maxit = 100, nr.reltol = sqrt(.Machine$double.eps), 
             calc.mcmc.se = TRUE, hessian = FALSE, compress = TRUE,
             maxNumDyadTypes = 1e+6, maxedges = 20000, maxchanges = 1e+06,
             maxMPLEsamplesize = 100000, MPLEtype=c("glm", "penalized"), trace = 0,
             steplength = 0.5, drop = TRUE, force.mcmc = FALSE, check.degeneracy=TRUE, mcmc.precision =
             0.05, metric = c("Likelihood", "raw"), method = c("BFGS", "Nelder-Mead"),
             trustregion = 20, initial.loglik = NULL, initial.network = NULL,
             style = c("Newton-Raphson", "Robbins-Monro",
             "Stochastic-Approximation"), phase1_n = NULL, initial_gain = NULL,
             nsubphases = "maxit", niterations = NULL, phase3_n = NULL,
             RobMon.phase1n_base = 7, RobMon.phase2n_base = 7, RobMon.phase2sub
             = 4, RobMon.init_gain = 0.4, RobMon.phase3n = 500, dyninterval =
             1000, packagenames="ergm", parallel = 0, returnMCMCstats = TRUE)
             }
\details{
    This function is only used within a call to the \code{\link{ergm}} function.
    See the \code{usage} section in \code{\link{ergm}} for details.
}
\arguments{
  \item{prop.weights}{ Specifies the method to allocate probabilities of
    being proposed to dyads. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints (in the \code{\link{ergm}} function).
    }
  \item{prop.args}{ An alternative, direct way of specifying additional arguments to proposal. }
  \item{nr.maxit}{ count; The maximum number of iterations in the
    Newton-Raphson optimization.  Defaults to \code{100}.
    \code{maxit} gives the total number of likelihood
    function evaluations. }
  \item{nr.reltol}{Relative convergence tolerance.  The Newton-Raphson optimization
    stops if it is unable to increase the log-likelihood by a factor of
    \code{reltol * (abs(log-likelihood) + reltol)} at a step.  Defaults to
    \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.}
  \item{calc.mcmc.se}{logical; should the contribution to the 
    standard errors of the estimator incurred by the MCMC sampling
    be computed. Default is \code{TRUE}.}
  \item{hessian}{logical; Should the Hessian matrix
    of the likelihood function be computed. 
    Default is \code{TRUE}.}
  \item{compress}{logical; Should the matrix of sample statistics
    returned be compressed to the set of unique statistics with a 
    column of frequencies post-pended.  This also uses a compression
    algorithm in the computation of the maximum psuedo-likelihood
    estimate that will dramatically speed it for large networks.
    Default is \code{FALSE}.}
  \item{maxNumDyadTypes}{count; The maximum number of unique
    pseudolikelihood change statistics to be allowed if \code{compress=TRUE}.
    It is only relevant in that case.
    Default is \code{1e+6}.}
  \item{maxedges}{ Maximum number of edges for which to allocate space. }
  \item{maxchanges}{ Maximum number of changes in dynamic network
    simulation for which to allocate space. }
  \item{maxMPLEsamplesize}{count; the determines the
    sample size to use for endogenous
    sampling in the pseudolikelihood computation. The sample size is the number
    of dyads in a network with \code{maxMPLEsamplesize} nodes.
    Default is \code{100000}.}
  \item{MPLEtype}{one of "glm" or "penalized"; method to use
    for psuedolikelihood. "glm" is the usual
    formal logistic regression. "penalized" uses the bias-reduced method
    of Firth (1993) as originally implemented by 
    Meinhard Ploner, Daniela Dunkler, Harry
    Southworth, and Georg Heinze in the "logistf" package. 
    Default is "glm".}
  \item{trace}{non-negative integer; If positive,
    tracing information on the
    progress of the optimization is produced. Higher values may
    produce more tracing information: for method \code{"L-BFGS-B"}
    there are six levels of tracing.  (To understand exactly what
    these do see the source code for \code{\link[stats]{optim}}: higher levels 
    give more detail.)}
  \item{steplength}{ Multiplier for step length, to make fitting more
    stable at the cost of efficiency. }
  \item{drop}{logical; Should the degenerate terms in the model be
    dropped from the fit?
    If statistics occur on the extreme of their range they
    correspond to infinite parameter estimates.
    Default is \code{TRUE}.}
  \item{force.mcmc}{logical; should MCMC maximum likelihood be used?  Only 
    relevant for dyadic independent networks, in which the MLE could be found
    using MPLE instead.}
  \item{check.degeneracy}{Logical:  Should the diagnostics include a 
  check for model degeneracy?}
  \item{mcmc.precision}{vector; upper bounds on the precision of the 
    standard errors induced by the MCMC algorithm.
    Defaults to \code{0.05}.}
  \item{metric}{character; The name of the optimization metric
    to use. Defaults to \code{"Likelihood"}.}
  \item{method}{character; The name of the optimization method
    to use. See \code{\link[stats]{optim}} for the options. The default method
    \code{"BFGS"} is a quasi-Newton method (also known as a variable
    metric algorithm). It is attributed to
    Broyden, Fletcher, Goldfarb and Shanno. This uses function values
    and gradients to build up a picture of the surface to be optimized.}
  \item{trustregion}{numeric; The maximum amount the algorithm will
    allow the approximated likelihood to be increased at a given iteration.
    Defaults to 20.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{initial.loglik}{Initial value of loglikelihood, if known.}
  \item{initial.network}{Initial network for MCMC, if different from observed
  network.}
  \item{style}{character; The style of maximum 
    likelihood estimation to use. The default is optimization of an
    MCMC estimate of the log-likelihood. An alternative is to use 
    a form of stochastic approximation (\code{"Robbins-Monro"}).
    The direct use of the likelihood function has many theoretical
    advantages over stochastic approximation, but the choice will
    depend on the model and data being fit. See
    Hunter and Handcock (2006) for details.}
  \item{phase1_n}{count; The number of MCMC samples to draw
    in Phase 1 of the stochastic approximation algorithm.
    Defaults to 7 plus 3 times the number of terms in the model.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{initial_gain}{numeric; The initial gain to
    Phase 2 of the stochastic approximation algorithm.
    Defaults to 0.1.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{nsubphases}{count; The number of sub-phases 
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to \code{maxit}.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{niterations}{count; The number of MCMC samples to draw
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to 7 plus the number of terms in the model.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{phase3_n}{count; The sample size for the MCMC sample
    in Phase 3 of the stochastic approximation algorithm.
    Defaults to 1000.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{RobMon.phase1n_base}{Robbins-Monro control parameter}
  \item{RobMon.phase2n_base}{Robbins-Monro control parameter}
  \item{RobMon.phase2sub}{Robbins-Monro control parameter}
  \item{RobMon.init_gain}{Robbins-Monro control parameter}
  \item{RobMon.phase3n}{Robbins-Monro control parameter}
  \item{returnMCMCstats}{logical; If this is \code{TRUE} (the
    default) the matrix of change 
    statistics from the MCMC run is returned as component \code{sample}.
    This matrix is actually an object of class \code{mcmc} and can be 
    used directly in the \code{CODA} package to assess MCMC
    convergence.}
  \item{dyninterval}{ Number of Metropolis-Hastings proposal for each
    phase in the dynamic network simulation. }
  \item{packagenames}{Names of packages in which changestatistics are found.}
  \item{parallel}{ Number of threads in which to run the sampling. }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{ergm}}. The \code{\link{control.simulate}} 
function performs a 
similar function for
\code{\link{simulate.ergm}}; 
\code{\link{control.gof}} performs a
similar function for \code{\link{gof}}.   }
\references{
\item
Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003),
StOCNET User\'s Manual. Version 1.4.

\item Firth (1993), Bias Reduction in Maximum Likelihood Estimates.  Biometrika, 80: 27-38.

\item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
exponential family models for networks. Journal of Computational
and Graphical Statistics, 15: 565-583.
}
\keyword{ models }

