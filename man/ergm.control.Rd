\name{ergm.control}
\alias{ergm.control}
\title{ Auxiliary for Controlling ERGM Fitting }
\description{
  Auxiliary function as user interface for fine-tuning 'ergm' fitting.
}
\usage{
ergm.control(prop.weights = "default", prop.args = NULL, prop.weights.diss = "default", prop.args.diss = NULL, nr.maxit = 100, calc.mcmc.se = TRUE, hessian = FALSE, compress = FALSE, maxNumDyadTypes = 10000, maxedges = 20000, maxchanges = 1e+06, MPLEsamplesize = 50000, trace = 0, boundDeg = NULL, steplength = 0.5, drop = TRUE, force.mcmc = FALSE, mcmc.precision = 0.05, metric = "Likelihood", method = "BFGS", trustregion = 20, style = "Newton-Raphson", phase1_n = NULL, initial_gain = NULL, nsubphases = "maxit", niterations = NULL, phase3_n = NULL, dyninterval = 1000, parallel = 0, returnMCMCstats = TRUE)
}
\arguments{
  \item{prop.weights}{ Specifies the method to allocate probabilities of
    being proposed to dyads. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. }
  \item{prop.args}{ An alternative, direct way of specifying additional arguments to proposal. }
  \item{prop.weights.diss}{ As \code{prop.weights}, for dissolution model. }
  \item{prop.args.diss}{ As \code{prop.args}, for dissolution model. }
  \item{nr.maxit}{ count; The maximum number of iterations in the
    Newton-Raphson optimization.  Defaults to \code{100}.
    \code{maxit} gives the total number of likelihood
    function evaluations. }
  \item{\code{calc.mcmc.se}}{logical; should the contribution to the 
    standard errors of the estimator incurred by the MCMC sampling
    be computed. Default is \code{TRUE}.}
  \item{\code{hessian}}{logical; Should the Hessian matrix
    of the likelihood function be computed. 
    Default is \code{TRUE}.}
  \item{\code{compress}}{logical; Should the matrix of sample statistics
    returned be compressed to the set of unique statistics with a 
    column of frequencies post-pended.  This also uses a compression
    algorithm in the computation of the maximum psuedo-likelihood
    estimate that will dramatically speed it for large networks.
    Default is \code{FALSE}.}
  \item{\code{maxNumDyadTypes}}{count; The maximum number of unique
    pseudolikelihood change statistics to be allowed if \code{compress=TRUE}.
    It is only relevant in that case.
    Default is \code{10000}.}
  \item{maxedges}{ Maximum number of edges for which to allocate space. }
  \item{maxchanges}{ Maximum number of changes in dynamic network
    simulation for which to allocate space. }
  \item{\code{MPLEsamplesize}}{count; the sample size to use for endogenous
    sampling in the pseudolikelihood computation.
    Default is \code{50000}.}
  \item{\code{trace}}{non-negative integer; If positive,
    tracing information on the
    progress of the optimization is produced. Higher values may
    produce more tracing information: for method \code{"L-BFGS-B"}
    there are six levels of tracing.  (To understand exactly what
    these do see the source code for \code{\link[stats]{optim}}: higher levels 
    give more detail.)}
  \item{\code{boundDeg}}{list; comprising up to five matrices and 
    one or more arguments, which are used to constrain the number of in 
    or out degrees that a node can have. By default no constraining 
    is applied.  This is explained more fully under the 
    \code{Conditioning within models} section of documentation for \code{\link{ergm}}. Defaults to no conditioning.}
  \item{steplength}{ Multiplier for step length, to make fitting more
    stable at the cost of efficiency. }
  \item{\code{drop}}{logical; Should the degenerate terms in the model be
    dropped from the fit?
    If statistics occur on the extreme of their range they
    correspond to infinite parameter estimates.
    Default is \code{TRUE}.}
  \item{force.mcmc}{ ~~Describe \code{force.mcmc} here~~ }
  \item{\code{mcmc.precision}}{vector; upper bounds on the precision of the 
    standard errors induced by the MCMC algorithm.
    Defaults to \code{0.05}.}
  \item{\code{metric}}{character; The name of the optimization metric
    to use. Defaults to \code{"Likelihood"}.}
  \item{\code{method}}{character; The name of the optimization method
    to use. See \code{\link[stats]{optim}} for the options. The default method
    \code{"BFGS"} is a quasi-Newton method (also known as a variable
    metric algorithm). It is attributed to
    Broyden, Fletcher, Goldfarb and Shanno. This uses function values
    and gradients to build up a picture of the surface to be optimized.}
  \item{\code{style}}{character; The style of maximum 
    likelihood estimation to use. The default is optimization of an
    MCMC estimate of the log-likelihood. An alternative is to use 
    a form of stochastic approximation (\code{"Robbins-Monro"}).
    The direct use of the likelihood function has many theoretical
    advantages over stochastic approximation, but the choice will
    depend on the model and data being fit. See Handcock (2000) and 
    Hunter and Handcock (2006) for details.}
  \item{\code{trustregion}}{numeric; The maximum amount the algorithm will
    allow the approximated likelihood to be increased at a given iteration.
    Defaults to 20.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{\code{phase1_n}}{count; The number of MCMC samples to draw
    in Phase 1 of the stochastic approximation algorithm.
    Defaults to 7 plus 3 times the number of terms in the model.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{\code{nsubphases}}{count; The number of sub-phases 
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to \code{maxit}.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{\code{initial_gain}}{numeric; The initial gain to
    Phase 2 of the stochastic approximation algorithm.
    Defaults to 0.1.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  
  
  \item{nsubphases}{ ~~Describe \code{nsubphases} here~~ }
  
  \item{\code{niterations}}{count; The number of MCMC samples to draw
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to 7 plus the number of terms in the model.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{\code{phase3_n}}{count; The sample size for the MCMC sample
    in Phase 3 of the stochastic approximation algorithm.
    Defaults to 1000.
    See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
  \item{\code{returnMCMCstats}}{logical; If this is \code{TRUE} (the
    default) the matrix of change 
    statistics from the MCMC run is returned as component \code{sample}.
    This matrix is actually an object of class \code{mcmc} and can be 
    used directly in the \code{CODA} package to assess MCMC
    convergence.}
  \item{dyninterval}{ Number of Metropolis-Hastings proposal for each
    phase in the dynamic network simulation. }
  \item{parallel}{ Number of threads in which to run the sampling. }
  \item{returnMCMCstats}{ Whether network statistics from the sampled
    networks should be returned. }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{ergm}}, \code{\link{glm.control}} performs a
  similar function for \code{\link{glm}} }

\keyword{ models }

