\name{control.stergm}
\alias{control.stergm}
\title{Auxiliary for Controlling Separable Temporal ERGM Fitting
}
\description{  Auxiliary function as user interface for fine-tuning 'stergm' fitting.

}
\usage{
control.stergm(prop.weights.form = "default", prop.args.form = NULL,
prop.weights.diss = "default", prop.args.diss = NULL, compress = FALSE,
SAN.burnin = 10000, SAN.interval = 1000, maxNumDyadTypes = 1e+06,
maxedges = 20000, maxchanges = 1e+06, initialfit = "MPLE", maxMPLEsamplesize = 1e+05, MPLEtype = c("glm", "penalized"), trace = 0, sequential = TRUE, style = c("Robbins-Monro", "SPSA", "SPSA2", "Nelder-Mead"), RM.phase1n_base = 7, RM.phase2n_base = 100, RM.phase2sub = 4, RM.init_gain = 0.5, RM.interval = 100, RM.burnin = 1000, SPSA.a = 1, SPSA.alpha = 0.602, SPSA.A = 100, SPSA.c = 1, SPSA.gamma = 0.101, SPSA.iterations = 1000, SPSA.interval = 1000, SPSA.burnin = 1000, NM.abstol = 0, NM.reltol = sqrt(.Machine$double.eps), NM.alpha = 1, NM.beta = 0.5, NM.gamma = 2, NM.maxit = 500, NM.interval = 1000, NM.burnin = 1000, packagenames = "ergm", parallel = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prop.weights.form}{Specifies the method to allocate probabilities of
    being proposed to dyads in the formation phase. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints (in the \code{\link{ergm}} function).
}
\item{prop.args.form}{
  An alternative, direct way of specifying additional arguments to
  the proposal in the formation phase.
}
  \item{prop.weights.diss}{As \code{prop.weights.form}, for the dissolution phase.
}
  \item{prop.args.diss}{As \code{prop.args.form}, for the dissolution phase.
}
  \item{compress}{logical; Should the matrix of sample statistics
    returned be compressed to the set of unique statistics with a 
    column of frequencies post-pended.  This also uses a compression
    algorithm in the computation of the maximum psuedo-likelihood
    estimate that will dramatically speed it for large networks.
    Default is \code{FALSE}.
}
  \item{SAN.burnin}{Burnin used for calling SAN routine.  If NULL,
  \code{burnin} is used.
}
  \item{SAN.interval}{
???
}
  \item{maxNumDyadTypes}{count; The maximum number of unique
    pseudolikelihood change statistics to be allowed if \code{compress=TRUE}.
    It is only relevant in that case.
    Default is \code{10000}.
}
  \item{maxedges}{Maximum number of edges for which to allocate space.
}
  \item{maxchanges}{Maximum number of changes in dynamic network
    simulation for which to allocate space.
  }

  \item{initialfit}{Method used to find the initial formation
  parameter. Defaults to \code{"MPLE"}.
  }
  \item{maxMPLEsamplesize}{count; the sample size to use for endogenous
    sampling in the pseudolikelihood computation.
    Default is \code{10^11}.
}
  \item{MPLEtype}{one of "glm" or "penalized"; method to use
    for psuedolikelihood. "glm" is the usual
    formal logistic regression. "penalized" uses the bias-reduced method
    of Firth (1993) as originally implemented by 
    Meinhard Ploner, Daniela Dunkler, Harry
    Southworth, and Georg Heinze in the "logistf" package. 
    Default is "glm".
}
  \item{trace}{non-negative integer; If positive,
    tracing information on the
    progress of the optimization is produced. Higher values may
    produce more tracing information: for method \code{"L-BFGS-B"}
    there are six levels of tracing.  (To understand exactly what
    these do see the source code for \code{\link[stats]{optim}}: higher levels 
    give more detail.)
}
  \item{sequential}{Should the next iteration of the fit use the last network 
    sampled as the starting network or always use the initially passed network?
    The results should be similar (stochastically), but the
    \code{sequential=TRUE} option if \code{meanstats} is far from the
    passed network.
}
  \item{style}{character; The style of method of moments estimation
    to use. The default is a form of stochastic approximation
    (\code{"Robbins-Monro"}), but it should only be used if it is known
    a priori that the derivative of each element of the equilibrium
    expected values of statistics of interest with respect to the
    corresponding formation phase parameter is positive. The other
    option, (\code{"SPSA"}) is less precise but does not make this
    assumption. \code{"SPSA2"} is multithreaded, but has more overhead,
    and may give wrong results, since \R's routines are not threadsafe.
}
  \item{RM.phase1n_base}{this helps define the 'phase1n' param, which in turn
                     multiplies 'RM.interval' to control the number of
                     phase1 iterations; this is the base portion of 'phase1n',
                     which is added to 3*(the number of formation coefficients)
}
  \item{RM.phase2n_base}{phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
}
\item{RM.phase2sub}{
  phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
}
\item{RM.init_gain}{
  this is only used to adjust 'aDdiaginv'in phase1,
                      in particular:
                             aDdiaginv = gain/sqrt(aDdiaginv)
                      default=0.5
}

  \item{RM.interval}{like the SPSA.interval, this seems a little more like
                      a sample size, than an interval, it helps control the
                      number of MCMCsteps used in phase1 and phase2; in
                      phase2, this limits the innermost loop counter; default=100
}
  \item{RM.burnin}{the number of MCMC steps to disregard for the burn-in
                      period; default=1000
}
  \item{SPSA.a,SPSA.alpha,SPSA.A,SPSA.c,SPSA.gamma}{ Controls the rate
                      of decay of SPSA gain as
                           SPSA.c/(i+1)^(SPSA.gamma)
                      where i is indexed from 0 to SPSA.iterations
}
  \item{SPSA.iterations}{the number of iterations to use in the SPSA optimization;
                      default=1000
}
  \item{SPSA.interval}{this is eventually received as 'S' and looks like a
                      a sample size, rather than an interval, since 'S' controls
                      the number of MCMC steps that contribute to the stats vector
}
  \item{SPSA.burnin}{the number of MCMC steps to disregard for the burnin
                      period; default=1000
}
  \item{NM.abstol}{
??
}
  \item{NM.reltol}{
??
}
  \item{NM.alpha}{
    ??
  }
  \item{NM.beta}{
??
}
  \item{NM.gamma}{
??
}
  \item{NM.maxit}{
??
}
  \item{NM.interval}{
??
}
  \item{NM.burnin}{
??
}
  \item{packagenames}{the packages in which change statistics are found; default="ergm"
}
  \item{parallel}{the number of threads in which to run sampling; default=0
}
}
\details{
      This function is only used within a call to the \code{\link{stergm}} function.
    See the \code{usage} section in \code{\link{stergm}} for details.
}
\value{
  A list with arguments as components.
}
\references{
    \itemize{
    \item
    Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003),
    StOCNET User\'s Manual. Version 1.4.
    
    \item Firth (1993),
    Bias Reduction in Maximum Likelihood Estimates.
    Biometrika,
    80: 27-38.

    \item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
    exponential family models for networks. Journal of Computational
    and Graphical Statistics, 15: 565-583.
    
    \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2010),
    A Steplength Algorithm for Fitting ERGMs, Penn State Department
    of Statistics Technical Report.
  }
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{stergm}}. The \code{\link{control.simulate.stergm}} 
function performs a 
similar function for
\code{\link{simulate.stergm}}.  
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (prop.weights.form = "default", prop.args.form = NULL, 
    prop.weights.diss = "default", prop.args.diss = NULL, compress = FALSE, 
    SAN.burnin = 10000, SAN.interval = 1000, maxNumDyadTypes = 1e+06, 
    maxedges = 20000, maxchanges = 1e+06, maxMPLEsamplesize = 1e+05, 
    MPLEtype = c("glm", "penalized"), trace = 0, sequential = TRUE, 
    style = c("Robbins-Monro", "SPSA", "SPSA2", 
        "Nelder-Mead"), RM.phase1n_base = 7, RM.phase2n_base = 100, 
    RM.phase2sub = 4, RM.init_gain = 0.5, RM.interval = 100, 
    RM.burnin = 1000, 
    SPSA.a = 1, SPSA.alpha = 0.602, SPSA.A = 100, SPSA.c = 1, 
    SPSA.gamma = 0.101, SPSA.iterations = 1000, SPSA.interval = 1000, 
    SPSA.burnin = 1000, NM.abstol = 0, NM.reltol = sqrt(.Machine$double.eps), 
    NM.alpha = 1, NM.beta = 0.5, NM.gamma = 2, NM.maxit = 500, 
    NM.interval = 1000, NM.burnin = 1000, packagenames = "ergm", 
    parallel = 0) 
{
    control <- list()
    for (arg in names(formals(sys.function()))) control[[arg]] <- get(arg)
    control$MPLEtype <- match.arg(MPLEtype)
    control$style <- match.arg(style)
    control
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
