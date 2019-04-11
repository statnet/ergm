#  File R/ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
###############################################################################
# The <ergm> function fits ergms from a specified formula returning either
# MPLEs or approximate MLE's based on MCMC estimation.
#
#
# --PARAMETERS--
#   formula       :  a formula of the form 'nw ~ model term(s)'
#   init        :  a vector of starting values for estimation or offset values, or optionally
#                    if these are to be estimated, NULL (the default);
#   constraints   :  a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     idegreedist
#                         observed  odegreedist
#                    default="~ ."
#   target.stats     :  a vector of the mean value parameters;
#                    default=the observed statistic from the 'nw' in formula
#   control       :  a list of control parameters returned from <control.ergm>;
#                    default=control.ergm()
#   verbose       :  whether ergm should be verbose (T or F); default=FALSE
#
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:                               
#   because an ergm object is the return type of several functions, and
#   because this is a rather lengthy list, this list represents the return
#   type for all funtions returning an ergm. The symbol preceding each
#   component indicates which function returns it:
#       <ergm>             = $
#       <ergm.mainfitloop> = *
#       <ergm.mple>        = !
#       <ergm.stepping>    = @
#       <ergm.stocapprox>  = %
#       <ergm.estimate>    = ^
#       <ergm.robmon>      = &
#       <ergm.mapl>        = #
#       <ergm.maple>       = ~
#
#   the components include:
#
#    $*!@%^&#~+  coef            :  the vector of estimated model coefficients
#    $* @%^&  +  sample          :  the row-binded matrix of network statistics from
#                                   each sample; 'sample' will also have 2 attributes:
#                   mcpar        : the following vector taken from control:
#                                        c(burnin+1, endrun, interval)
#                                  where 'endrun' is defined as
#                                     burnin+interval*(samplesize-1)
#                   class        : "mcmc"
#    $*!@%^&#~+  iterations      :  the number of Newton-Raphson iterations required
#                                   before convergence
#    $*!@%^&#~+  MCMCtheta       :  the vector of natural parameters used to produce
#                                   the MCMC sample
#    $* @%^&  +  loglikelihood   :  the estimated change in log-likelihood in the last
#                                   iteration
#    $*!@%^&#~+  gradient        :  the value of the gradient of the approximated log-
#                                   likelihood function at the maximizing value
#    $*!@%^&#~+  covar           :  the approximated covariance matrix for the MLE
#    $*!@%^&#~+  samplesize      :  the size of the MCMC sample
#    $*!@%^&#~+  failure         :  whether estimation failed (T or F)
#    $*!@%^&#~+  mc.se           :  the standard error estimates
#    $* @% &# +  newnetwork      :  the final network sampled; in the ergm returned from
#                                   <ergm.robmom>, this='network'
#    $* @% &# +  network         :  the 'nw' inputted to <ergm> via the 'formula'
#    $* @% &  +  theta.original  :  the theta values at the start of the MCMC sampling
#    $* @        mplefit         :  the MPLE fit as a glm object, and returned by
#                                   <ergm.mple>
#    $*!@% &#~+  mle.lik         :  the approximate log-likelihood for the MLE, if computed
#    $* @        etamap          :  the set of function mapping theta -> eta;
#                                   see <etamap>? for the components of this list
#    $           degeneracy.value:  the degeneracy value assigned by <ergm.degeneracy>
#    $           degeneracy.type :  a vector of length 2, as returned by
#                                   <ergm.compute.degeneracy> (found in the
#                                   <ergm.degeracy> file)
#    $      #    formula         :  the 'formula' value inputted to <ergm>
#    $      #    constraints     :  the 'constraints' value inputted to <ergm>
#           #    prop.args       :  the list of arguments that were passed onto the
#                                   <InitErgmProposal> routines
#    $      #    prop.weights    :  the MCMC proposal weights inputted to <ergm> via
#                                  'control'
#    $      #    offset          :  a vector of whether each model parameter was set at
#                                  a fixed value (not estimated)
#    $      #    drop            :  list of dropped terms
#     * @%^&  +  sample.obs      :  the matrix of sample network statistics for observed
#                                   data
#     * @        parallel        :  the number of additional threads used when sampling
#      !    #~   glm             :  the fit established by MPL estimation and returned
#                                   by <ergm.logitreg>, <ergm.pen.glm> or <glm>
#                                   depending on the 'MPLEtype';
#      !    #~   glm.null        :  the null fit established by MPL estimation and
#                                   returned by <ergm.logitreg>, <ergm.pen.glm> or <glm>
#                                   depending on the 'MPLEtype';
#      !   #~    theta1          :  the vector of ??
#         &      rm.coef         :  the robmon coefficients used as 'init' in the final
#                                   estimation
#      !   #~   loglikelihoodratio: the log-likelihood corresponding to
#                                   'coef'
#
#####################################################################################    

#' Exponential-Family Random Graph Models
#'
#' \code{\link{ergm}} is used to fit exponential-family random graph
#' models (ERGMs), in which
#' the probability of a given network, \eqn{y}, on a set of nodes is 
#' \eqn{h(y) \exp\{\eta(\theta) \cdot
#' g(y)\}/c(\theta)}, where
#' \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}),
#' \eqn{g(y)} is a vector of network statistics for \eqn{y},
#' \eqn{\eta(\theta)} is a natural parameter vector of the same 
#' length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
#' normalizing constant for the distribution.
#' \code{\link{ergm}} can return a maximum pseudo-likelihood
#' estimate, an approximate maximum likelihood estimate based on a Monte
#' Carlo scheme, or an approximate contrastive divergence estimate based
#' on a similar scheme.
#' (For an overview of the package, see \code{\link{ergm-package}}.)
#' 
#' @param formula {An \R \code{\link{formula}} object, of the form
#' \code{y ~ <model terms>},
#' where \code{y} is a \code{\link[network]{network}} object or a matrix that can be
#' coerced to a \code{\link[network]{network}}  object.  For the details on the possible
#' \code{<model terms>}, see \code{\link{ergm-terms}} and Morris, Handcock and
#' Hunter (2008) for binary ERGM terms and
#' Krivitsky (2012) for valued ERGM
#' terms (terms for weighted edges).  To create a
#' \code{\link[network]{network}} object in \R, use the \code{network()} function,
#' then add nodal attributes to it using the \code{\%v\%}
#' operator if necessary. Enclosing a model term in \code{offset()}
#' fixes its value to one specified in \code{offset.coef}.
#' }
#' @template response
#' @param reference {A one-sided formula specifying
#' the reference measure (\eqn{h(y)}) to be used. (Defaults to \code{~Bernoulli}.)
#' See help for [ERGM reference measures][ergm-references] implemented in the
#' **[ergm][ergm-package]** package.}
#' 
#' @param constraints {A formula specifying one or more constraints
#' on the support of the distribution of the networks being modeled,
#' using syntax similar to the \code{formula} argument, on the
#' right-hand side. Multiple constraints
#' may be given, separated by \dQuote{+} and \dQuote{-} operators. (See
#' [ERGM constraints][ergm-constraints] for the explanation of
#' their semantics.)
#' Together with the model terms in the formula and the reference measure, the constraints
#' define the distribution of networks being modeled.
#' 
#' It is also possible to specify a proposal function directly
#' either by passing a string with the function's name (in which case,
#' arguments to the proposal should be specified through the
#' \code{prop.args} argument to \code{\link{control.ergm}}) or by
#' giving it on the LHS of the constraints formula, in which case it
#' will override the one chosen automatically.
#' 
#' The default is \code{~.}, for an unconstrained model.
#' 
#' See the [ERGM constraints][ergm-constraints] documentation for
#' the constraints implemented in the **[ergm][ergm-package]**
#' package. Other packages may add their own constraints.
#' 
#' Note that not all possible combinations of constraints and reference
#' measures are supported. However, for relatively simple constraints
#' (i.e., those that simply permit or forbid specific dyads or sets of
#' dyads from changing), arbitrary combinations should be possible.
#' }
#' @param offset.coef {A vector of coefficients for the offset terms.}
#' @param target.stats {vector of "observed network statistics,"
#' if these statistics are for some reason different than the 
#' actual statistics of the network on the left-hand side of
#' \code{formula}.
#' Equivalently, this vector is the mean-value parameter values for the
#' model.  If this is given, the algorithm finds the natural
#' parameter values corresponding to these mean-value parameters.
#' If \code{NULL}, the mean-value parameters used are the observed
#' statistics of the network in the formula.
#' }
#' @param eval.loglik {Logical:  For dyad-dependent models, if TRUE, use bridge
#' sampling to evaluate the log-likelihoood associated with the
#' fit. Has no effect for dyad-independent models.
#' Since bridge sampling takes additional time, setting to FALSE may
#' speed performance if likelihood values (and likelihood-based
#' values like AIC and BIC) are not needed. Can be set globally via `option(ergm.eval.loglik=...)`, which is set to `TRUE` when the package is loaded.
#' }
#' @param estimate {If "MPLE," then the maximum pseudolikelihood estimator
#' is returned.  If "MLE" (the default), then an approximate maximum likelihood
#' estimator is returned.  For certain models, the MPLE and MLE are equivalent,
#' in which case this argument is ignored.  (To force MCMC-based approximate
#' likelihood calculation even when the MLE and MPLE are the same, see the
#' \code{force.main} argument of \code{\link{control.ergm}}. If "CD" (\emph{EXPERIMENTAL}),
#' the Monte-Carlo contrastive divergence estimate is returned. )
#' }
#' @param control {A list of control parameters for algorithm
#' tuning. Constructed using \code{\link{control.ergm}}. 
#' }
#' @param verbose {logical; if this is
#' \code{TRUE}, the program will print out additional
#' information, including goodness of fit statistics.
#' }
#' @param \dots {Additional
#' arguments, to be passed to lower-level functions.
#' }
#' 
#' @return
#' \code{\link{ergm}} returns an object of class \code{\link{ergm}} that is a list
#' consisting of the following elements:
#' \item{coef}{The Monte Carlo maximum likelihood estimate
#' of \eqn{\theta}, the vector of coefficients for the model
#' parameters.}
#' \item{sample}{The \eqn{n\times p} matrix of network statistics, 
#' where \eqn{n} is the                               
#' sample size and \eqn{p} is the number of network statistics specified in the
#' model, generated by the last iteration of the MCMC-based likelihood maximization routine. These statistics are centered with respect to the observed statistics or `target.stats`, unless missing data MLE is used.}
#' \item{sample.obs}{As \code{sample}, but for the constrained sample.}
#' \item{iterations}{The number of Newton-Raphson iterations required
#' before convergence.}
#' \item{MCMCtheta}{The value of \eqn{\theta} used to produce the Markov chain
#' Monte Carlo sample.  As long as the Markov chain mixes sufficiently
#' well, \code{sample} is roughly a random sample from the distribution
#' of network statistics specified by the model with the parameter equal
#' to \code{MCMCtheta}.  If \code{estimate="MPLE"} then 
#' \code{MCMCtheta} equals the MPLE.}
#' \item{loglikelihood}{The approximate change in log-likelihood 
#' in the last iteration.
#' The value is only approximate because it is estimated based 
#' on the MCMC random sample.}
#' \item{gradient}{The value of the gradient vector of the approximated
#' loglikelihood function, evaluated at the maximizer.  This vector
#' should be very close to zero.}
#' \item{covar}{Approximate covariance matrix for the MLE, based on the inverse
#' Hessian of the approximated loglikelihood evaluated at the maximizer.}
#' \item{failure}{Logical:  Did the MCMC estimation fail?}
#' \item{network}{Original network}
#' \item{newnetworks}{A list of the final networks at the end of the MCMC
#' simulation, one for each thread.}
#' \item{newnetwork}{The first (possibly only) element of \code{netwonetworks}.}
#' \item{coef.init}{The initial value of \eqn{\theta}.}
#' \item{est.cov}{The covariance matrix of the model statistics in the final MCMC sample.}
#' \item{coef.hist, steplen.hist, stats.hist, stats.obs.hist}{
#' For the MCMLE method, the history of coefficients, Hummel step lengths, and average model statistics for each iteration..
#' }
#' \item{control}{The control list passed to the call.}
#' \item{etamap}{The set of functions mapping the true parameter theta
#' to the canonical parameter eta (irrelevant except in a curved exponential
#' family model)}
#' \item{formula}{The original \code{\link{formula}} entered into the \code{\link{ergm}} function.}
#' \item{target.stats}{The target.stats used during estimation (passed through from the Arguments)}
#' \item{target.esteq}{Used for curved models to preserve the target mean values of the curved terms. It is identical to target.stats for non-curved models.}
#' \item{constrained}{The list of constraints implied by the constraints used by original \code{ergm} call}
#' \item{constraints}{Constraints used during estimation (passed through from the Arguments)}
#' \item{reference}{The reference measure used during estimation (passed through from the Arguments)}
#' \item{estimate}{The estimation method used (passed through from the Arguments).}
#' \item{offset}{vector of logical telling which model parameters are to be set
#' at a fixed value (i.e., not estimated).}
#' 
#' \item{drop}{If \code{\link[=control.ergm]{control$drop}=TRUE}, a numeric vector indicating which terms were dropped due to to extreme values of the
#' corresponding statistics on the observed network, and how:
#' \describe{
#' \item{\code{0}}{The term was not dropped.}
#' \item{\code{-1}}{The term was at its minimum and the coefficient was fixed at
#' \code{-Inf}.}
#' \item{\code{+1}}{The term was at its maximum and the coefficient was fixed at
#' \code{+Inf}.}
#' }
#' }
#' 
#' \item{estimable}{A logical vector indicating which terms could not be
#' estimated due to a \code{constraints} constraint fixing that term at a
#' constant value.
#' }
#' 
#' \item{null.lik}{Log-likelihood of the null model. Valid only for
#' unconstrained models.}
#' \item{mle.lik}{The approximate log-likelihood for the MLE.
#' The value is only approximate because it is estimated based 
#' on the MCMC random sample.}
#' 
#' \item{degeneracy.value}{Score calculated to assess the degree of 
#' degeneracy in the model. Only shows when MCMLE.check.degeneracy is TRUE in \code{control.ergm}. }
#' \item{degeneracy.type}{Supporting output for \code{degeneracy.value}. Only shows when MCMLE.check.degeneracy is TRUE in \code{control.ergm}. Mainly for internal use.}
#' 
#' See the method \code{\link{print.ergm}} for details on how
#' an \code{\link{ergm}} object is printed.  Note that the
#' method \code{\link{summary.ergm}} returns a summary of the
#' relevant parts of the \code{\link{ergm}} object in concise summary
#' format.
#' 
#' @section Notes on model specification:
#' Although each of the statistics in a given model is a summary
#' statistic for the entire network, it is rarely necessary to
#' calculate statistics for an entire network
#' in a proposed Metropolis-Hastings step.
#' Thus, for example, if the triangle term is included in the model,
#' a census of all triangles in the observed network is never
#' taken; instead, only the change in the number of triangles
#' is recorded for each edge toggle.
#' 
#' In the implementation of \code{\link{ergm}}, the model is initialized
#' in \R, then all the model information is passed to a C program
#' that generates the sample of network statistics using MCMC.
#' This sample is then returned to \R, which implements a
#' simple Newton-Raphson algorithm to approximate the MLE.
#' An alternative style of maximum likelihood estimation is to use a stochastic
#' approximation algorithm. This can be chosen with the 
#' \code{control.ergm(style="Robbins-Monro")} option.
#' 
#' The mechanism for proposing new networks for the MCMC sampling
#' scheme, which is a Metropolis-Hastings algorithm, depends on 
#' two things:  The \code{constraints}, which define the set of possible
#' networks that could be proposed in a particular Markov chain step,
#' and the weights placed on these possible steps by the 
#' proposal distribution.  The former may be controlled using the
#' \code{constraints} argument described above.  The latter may
#' be controlled using the \code{prop.weights} argument to the
#' \code{\link{control.ergm}} function.
#' 
#' The package is designed so that the user could conceivably add additional 
#' proposal types. 
#' 
#' @references
#' Admiraal R, Handcock MS (2007).
#' \pkg{networksis}: Simulate bipartite graphs with fixed
#' marginals through sequential importance sampling.
#' Statnet Project, Seattle, WA.
#' Version 1. \url{statnet.org}.
#' 
#' Bender-deMoll S, Morris M, Moody J (2008).
#' Prototype Packages for Managing and Animating Longitudinal
#' Network Data: \pkg{dynamicnetwork} and \pkg{rSoNIA}.
#' \emph{Journal of Statistical Software}, 24(7).
#' \url{https://www.jstatsoft.org/v24/i07/}.
#' 
#' 
#' Butts CT (2007).
#' \pkg{sna}: Tools for Social Network Analysis.
#' R package version 2.3-2. \url{https://cran.r-project.org/package=sna}.
#' 
#' Butts CT (2008).
#' \pkg{network}: A Package for Managing Relational Data in \R.
#' \emph{Journal of Statistical Software}, 24(2).
#' \url{https://www.jstatsoft.org/v24/i02/}.
#' 
#' Butts C (2015).
#' \pkg{network}: The Statnet Project (https://statnet.org). R package version 1.12.0, \url{https://cran.r-project.org/package=network}.
#' 
#' Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
#' A \pkg{statnet} Tutorial.
#' \emph{Journal of Statistical Software}, 24(8).
#' \url{https://www.jstatsoft.org/v24/i08/}.
#' 
#' Goodreau SM, Kitts J, Morris M (2008b).
#' Birds of a Feather, or Friend of a Friend? Using Exponential
#' Random Graph Models to Investigate Adolescent Social Networks.
#' \emph{Demography}, 45, in press.
#' 
#' Handcock, M. S. (2003)
#' \emph{Assessing Degeneracy in Statistical Models of Social Networks},
#' Working Paper \#39, 
#' Center for Statistics and the Social Sciences,
#' University of Washington.
#' \url{www.csss.washington.edu/Papers/wp39.pdf}
#' 
#' Handcock MS (2003b).
#' \pkg{degreenet}: Models for Skewed Count Distributions Relevant
#' to Networks.
#' Statnet Project, Seattle, WA.
#' Version 1.0, \url{statnet.org}.
#' 
#' Handcock MS and Gile KJ (2010). Modeling Social Networks from Sampled Data. \emph{Annals of Applied Statistics}, 4(1), 5-25. \doi{10.1214/08-AOAS221}
#' 
#' Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003a).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks.
#' Statnet Project, Seattle, WA.
#' Version 2, \url{statnet.org}.
#' 
#' Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
#' \pkg{statnet}: Software Tools for the Statistical Modeling of
#' Network Data.
#' Statnet Project, Seattle, WA.
#' Version 2, \url{statnet.org}.
#' 
#' Hunter, D. R. and Handcock, M. S. (2006)
#' \emph{Inference in curved exponential family models for networks},
#' Journal of Computational and Graphical Statistics.
#' 
#' Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks.
#' \emph{Journal of Statistical Software}, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#' 
#' Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#' 
#' Morris M, Handcock MS, Hunter DR (2008).
#' Specification of Exponential-Family Random Graph Models:
#' Terms and Computational Aspects.
#' \emph{Journal of Statistical Software}, 24(4).
#' \url{https://www.jstatsoft.org/v24/i04/}.
#' 
#' Snijders, T.A.B. (2002),
#' Markov Chain Monte Carlo Estimation of Exponential Random Graph Models.
#' Journal of Social Structure.
#' Available from 
#' \url{https://www.cmu.edu/joss/content/articles/volume3/Snijders.pdf}.
#' 
#' @seealso network, \%v\%, \%n\%, \code{\link{ergm-terms}}, \code{\link{ergmMPLE}},
#' \code{\link{summary.ergm}}, \code{\link{print.ergm}}
#' 
#' @examples
#' \donttest{
#' #
#' # load the Florentine marriage data matrix
#' #
#' data(flo)
#' #
#' # attach the sociomatrix for the Florentine marriage data
#' # This is not yet a network object.
#' #
#' flo
#' #
#' # Create a network object out of the adjacency matrix
#' #
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage
#' #
#' # print out the sociomatrix for the Florentine marriage data
#' #
#' flomarriage[,]
#' #
#' # create a vector indicating the wealth of each family (in thousands of lira) 
#' # and add it as a covariate to the network object
#' #
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#' flomarriage
#' #
#' # create a plot of the social network
#' #
#' plot(flomarriage)
#' #
#' # now make the vertex size proportional to their wealth
#' #
#' plot(flomarriage, vertex.cex=flomarriage %v% "wealth" / 20, main="Marriage Ties")
#' #
#' # Use 'data(package = "ergm")' to list the data sets in a
#' #
#' data(package="ergm")
#' #
#' # Load a network object of the Florentine data
#' #
#' data(florentine)
#' #
#' # Fit a model where the propensity to form ties between
#' # families depends on the absolute difference in wealth
#' #
#' gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
#' summary(gest)
#' #
#' # add terms for the propensity to form 2-stars and triangles
#' # of families 
#' #
#' gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#' summary(gest)
#' 
#' # import synthetic network that looks like a molecule
#' data(molecule)
#' # Add a attribute to it to mimic the atomic type
#' molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#' #
#' # create a plot of the social network
#' # colored by atomic type
#' #
#' plot(molecule, vertex.col="atomic type",vertex.cex=3)
#' 
#' # measure tendency to match within each atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"),
#' control=control.ergm(MCMC.samplesize=10000))
#' summary(gest)
#' 
#' # compare it to differential homophily by atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle
#'                         + nodematch("atomic type",diff=TRUE),
#' control=control.ergm(MCMC.samplesize=10000))
#' summary(gest)
#' }
#' @keywords models
#' @aliases is.ergm ergm.object
#' @export
ergm <- function(formula, response=NULL,
                 reference=~Bernoulli,
                 constraints=~.,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=getOption("ergm.eval.loglik"),
                 estimate=c("MLE", "MPLE", "CD"),
                 control=control.ergm(),
                 verbose=FALSE,...) {
  check.control.class("ergm", "ergm")
  control.toplevel(control,...)
  
  estimate <- match.arg(estimate)

  if(estimate=="CD"){
    control$init.method <- "CD"
    eval.loglik <- FALSE
  }

  if(estimate=="MPLE"){
    control$init.method <- "MPLE"
  }
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))
  if (verbose) message("Evaluating network in model.")
  
  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"
  
  
  # Missing data handling only needs to happen if the sufficient
  # statistics are not specified. If the sufficient statistics are
  # specified, the nw's dyad states are irrelevant.
  if(network.naedgecount(nw) && !is.null(target.stats)){
    warning("Target statistics specified in a network with missing dyads. Missingness will be overridden.")
    nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
  }
  
  proposal.obs <- if(network.naedgecount(nw)==0) NULL else append_rhs.formula(constraints, list(as.name("observed")), TRUE)

  if (verbose) message("Initializing Metropolis-Hastings proposal(s):",appendLF=FALSE) 
  
  ## FIXME: a more general framework is needed?
  if(!is.null(response) && reference==~Bernoulli){
    warn(paste0("The default Bernoulli reference distribution operates in the binary (",sQuote("response=NULL"),") mode only. Did you specify the ",sQuote("reference")," argument?"))
  }
  
  proposal <- ergm_proposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass,reference=reference,response=response)
  if (verbose) message(" ",proposal$pkgname,":MH_",proposal$name)
  
  
  if(!is.null(proposal.obs)){
    proposal.obs <- ergm_proposal(proposal.obs, weights=control$obs.MCMC.prop.weights, control$obs.MCMC.prop.args, nw, class=proposalclass, reference=reference, response=response)
    if (verbose) message(" ",proposal.obs$pkgname,":MH_",proposal.obs$name)
  }
  
  if (verbose) message("Initializing model.")
  
  # Construct the initial model.
  
  # The following kludge knocks out MPLE if the sample space
  # constraints are not dyad-independent. For example, ~observed
  # constraint is dyad-independent, while ~edges is not.
  #
  # TODO: Create a flexible and general framework to manage methods
  # for obtaining initial values.
  init.candidates <- proposal$reference$init_methods
  if("MPLE" %in% init.candidates && !is.dyad.independent(proposal$arguments$constraints,
                                                         proposal.obs$arguments$constraints)){
    init.candidates <- init.candidates[init.candidates!="MPLE"]
    if(verbose) message("MPLE cannot be used for this constraint structure.")
  }
  if("MPLE" %in% init.candidates && !is.null(target.stats) && is.curved(formula, response=response, term.options=control$term.options)){
    init.candidates <- init.candidates[init.candidates!="MPLE"]
    if(verbose) message("At this time, MPLE cannot be used for curved families when target.stats are passed.")
  }
  control$init.method <- ERRVL(try(match.arg(control$init.method, init.candidates), silent=TRUE), {
    message("Sepcified initial parameter method ", sQuote(control$init.method), " is not in the list of candidates. Use at your own risk.")
    control$init.method
  })
  if(verbose) message(paste0("Using initial method '",control$init.method,"'."))
  model.initial <- ergm_model(formula, nw, response=response, initialfit=control$init.method=="MPLE", term.options=control$term.options)
  
  ## Construct approximate response network if target.stats are given.
  
  if(!is.null(target.stats)){
    formula.no <- filter_rhs.formula(formula, function(x) (if(is.call(x)) x[[1]] else x)!="offset")
    nw.stats<-summary(formula.no,response=response, term.options=control$term.options)
    target.stats <- vector.namesmatch(target.stats, names(nw.stats))
    target.stats <- na.omit(target.stats)
    if(length(nw.stats)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    
    # no need to pass the offset term's init to SAN
    offset.terms <- model.initial$etamap$offsettheta
    san.control <- control$SAN.control
    
    if(verbose) message("Constructing an approximate response network.")
    ## If target.stats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    if(control$SAN.maxit > 0){
      TARGET_STATS<-san(formula.no, target.stats=target.stats,
                response=response,
                reference=reference,
                constraints=constraints,
                control=san.control,
                only.last=TRUE,
                output="pending_update_network",
                verbose=verbose)
      if(verbose) message("Finished SAN run.")
    }else{
      TARGET_STATS <- nw
    }
    nw <- TARGET_STATS <- as.network(TARGET_STATS)
    formula<-nonsimp_update.formula(formula,TARGET_STATS~., from.new="TARGET_STATS")
    offinfo <- offset.info.formula(formula,response=response,term.options=control$term.options)
    tmp <- rep(NA, length(offinfo$eta))
    tmp[!offinfo$eta] <- target.stats
    names(tmp)[!offinfo$eta] <- names(target.stats)
    nw.stats <- summary(formula,response=response, term.options=control$term.options)
    names(tmp)[offinfo$eta] <- names(nw.stats)[offinfo$eta]
    
    # From this point on, target.stats has NAs corresponding to the
    # offset terms.
    #
    # TODO: Only have target.stats contain non-offset terms'
    # statistics, and have the rest of the code handle it
    # intelligently.
    target.stats <- tmp
  } else {
    if (network.edgecount(nw) == 0) warning("Network is empty and no target stats are specified.")
  } 
  
  # If some control$init is specified...
  if(!is.null(control$init)){
    # Check length of control$init.
    if (length(control$init)!=length(model.initial$etamap$offsettheta)) {
      if(verbose){
        message("control$init =")
        message_print(control$init)
        message("number of statistics is ",length(model.initial$coef.names), "")
      }
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.",
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  }else control$init <- rep(NA, length(model.initial$etamap$offsettheta)) # Set the default value of control$init.
  
  if(!is.null(offset.coef)){
    # TODO: Names matching here?
    if(length(control$init[model.initial$etamap$offsettheta])!=length(offset.coef))
      stop("Invalid offset parameter vector offset.coef: ",
           "wrong number of parameters: expected ",
           length(control$init[model.initial$etamap$offsettheta]),
           " got ",length(offset.coef),".")
    control$init[model.initial$etamap$offsettheta]<-offset.coef
  }
  
  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model.initial$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(model.initial$coef.names[is.na(control$init)|model.initial$offsettheta]), ".", sep="")
  
  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model.initial, proposal, control$init)
  model.initial <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model.initial, nw=nw, init=control$init, response=response, target.stats=target.stats, drop=control$drop)
  model.initial <- extremecheck$model; control$init <- extremecheck$init
  
  
  # Construct the curved model, and check if it's different from the initial model. If so, we know that it's curved.
  model <- ergm_model(formula, nw, response=response, expanded=TRUE, silent=TRUE, term.options=control$term.options)
  # MPLE is not supported for curved ERGMs.
  if(estimate=="MPLE"){
    if(!is.null(response)) stop("Maximum Pseudo-Likelihood (MPLE) estimation for valued ERGMs is not implemented at this time. You may want to pass fixed=TRUE parameter in curved terms to specify the curved parameters as fixed.")
    if(length(model$etamap$offsetmap)!=length(model.initial$etamap$offsetmap)) stop("Maximum Pseudo-Likelihood (MPLE) estimation for curved ERGMs is not implemented at this time. You may want to pass fixed=TRUE parameter in curved terms to specify the curved parameters as fixed.")
    ## if(!is.dyad.independent(proposal$arguments$constraints,
    ##                         proposal.obs$arguments$constraints))
    ##   stop("Maximum Pseudo-Likelihood (MPLE) estimation for ERGMs with dyad-dependent constraints is only implemented for certain degree constraints at this time.")
  }
  
  if (verbose) { message("Fitting initial model.") }
  
  MPLE.is.MLE <- (proposal$reference$name=="Bernoulli"
                  && is.dyad.independent(model.initial)
                  && !is.curved(formula, response=response, term.options=control$term.options)
                  && !control$force.main
                  && is.dyad.independent(proposal$arguments$constraints,
                                         proposal.obs$arguments$constraints))
  
  # If all other criteria for MPLE=MLE are met, _and_ SAN network matches target.stats directly, we can get away with MPLE.
  if (!is.null(target.stats) && !isTRUE(all.equal(target.stats[!is.na(target.stats)],nw.stats[!is.na(target.stats)]))) message("Unable to match target stats. Using MCMLE estimation.")
  MCMCflag <- (estimate=="MLE" && (!MPLE.is.MLE
                                   || (!is.null(target.stats) && !isTRUE(all.equal(target.stats,nw.stats)))
  )
  || control$force.main)
  
  # Short-circuit the optimization if all terms are either offsets or dropped.
  if(all(model.initial$etamap$offsettheta)){
    # Note that this cannot be overridden with control$force.main.
    message("All terms are either offsets or extreme values. No optimization is performed.")
    return(structure(list(coef=control$init,
                          iterations=0,
                          loglikelihood=NA,
                          mle.lik=NULL,
                          gradient=rep(NA,length=length(control$init)),
                          failure=TRUE,
                          offset=model.initial$etamap$offsettheta,
                          drop=if(control$drop) extremecheck$extremeval.theta,
                          estimable=constrcheck$estimable,
                          network=nw,
                          reference=reference,
                          response=response,
                          newnetwork=nw,
                          formula=formula,
                          constrained=proposal$arguments$constraints,
                          constrained.obs=proposal.obs$arguments$constraints,
                          constraints=constraints,
                          target.stats=model.initial$target.stats,
                          target.esteq=if(!is.null(model.initial$target.stats)) ergm.estfun(rbind(model.initial$target.stats), initialfit$coef, model.initial),
                          estimate=estimate,
                          ergm_version=packageVersion("ergm"),
                          control=control
    ),
    class="ergm"))
    
  }
  
  model.initial$nw.stats <- summary(model.initial, nw=nw, response=response, initialfit=control$init.method=="MPLE", term.options=control$term.options)
  model.initial$target.stats <- NVL(target.stats, model.initial$nw.stats)
  
  if(control$init.method=="CD") if(is.null(names(control$init)))
    names(control$init) <- param_names(model.initial, FALSE)
  
  initialfit <- ergm.initialfit(init=control$init, initial.is.final=!MCMCflag,
                                formula=formula, nw=nw, reference=reference, 
                                m=model.initial, method=control$init.method,
                                MPLEtype=control$MPLE.type, 
                                control=control,
                                proposal=proposal,
                                proposal.obs=proposal.obs,
                                verbose=verbose, response=response,
                                maxNumDyadTypes=control$MPLE.max.dyad.types,
                                ...)
  
  if (!MCMCflag){ # Just return initial (non-MLE) fit and exit.
    message("Stopping at the initial estimate.")
    initialfit$MPLE_is_MLE <- MPLE.is.MLE
    initialfit$ergm_version <- packageVersion("ergm")
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- if(control$drop) extremecheck$extremeval.theta
    initialfit$estimable <- constrcheck$estimable
    initialfit$network <- nw
    initialfit$reference <- reference
    initialfit$response <- response
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constrained <- proposal$arguments$constraints
    initialfit$constrained.obs <- proposal.obs$arguments$constraints
    initialfit$constraints <- constraints
    initialfit$target.stats <- model.initial$target.stats
    initialfit$etamap <- model.initial$etamap
    initialfit$target.esteq <- if(!is.null(model.initial$target.stats)) ergm.estfun(rbind(model.initial$target.stats), initialfit$coef, model.initial)
    initialfit$estimate <- estimate
    
    initialfit$control<-control
    
    if(eval.loglik) initialfit$null.lik <- logLikNull.ergm(initialfit, verbose=verbose)
    if(any(!model.initial$etamap$offsettheta) && eval.loglik){
      message("Evaluating log-likelihood at the estimate. ",appendLF=FALSE)
      initialfit<-logLik(initialfit, add=TRUE, control=control$loglik.control, verbose=verbose)
      message("")
    }
    return(initialfit)
  }
  
  # Otherwise, set up the main phase of estimation:
  
  ergm.getCluster(control, max(verbose-1,0))
  
  # Revise the initial value, if necessary:
  init <- initialfit$coef
  init[is.na(init)] <- 0
  if(control$init.method=="MPLE"){ # Only MPLE requires these kludges.
    names(init) <- model.initial$coef.names
    # revise init to reflect additional parameters
    init <- ergm.reviseinit(model, init)
  }
  names(init) <- param_names(model, FALSE)
  
  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, proposal, init=init, silent=TRUE)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=init, response=response, target.stats=target.stats, drop=control$drop, silent=TRUE)
  model <- extremecheck$model; init <- extremecheck$init
  
  model$nw.stats <- summary(model, nw=nw, response=response, term.options=control$term.options)
  model$target.stats <- NVL(target.stats, model$nw.stats)
  
  mainfit <- switch(control$main.method,
                    "Robbins-Monro" = ergm.robmon(init, nw, model, 
                                                  proposal=proposal, verbose=verbose, control=control),
                    "Stochastic-Approximation" = ergm.stocapprox(init, nw, model, 
                                                                 control=control, proposal=proposal,
                                                                 verbose),
                    "Stepping" = ergm.stepping(init, nw, model, initialfit, constraints,
                                               #nstats=nstats, 
                                               #approx=lognormapprox, filename.prefix=NULL, 
                                               #control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
                                               #plots=FALSE,  # currently useless, but plots can be reimplemented
                                               control=control, 
                                               proposal=proposal, proposal.obs=proposal.obs, 
                                               verbose=verbose,...),
                    "MCMLE" = ergm.MCMLE(init, nw,
                                         model, 
                                         # no need to pass initialfit to MCMLE
                                         initialfit=(initialfit<-NULL),
                                         control=control, proposal=proposal,
                                         proposal.obs=proposal.obs,
                                         verbose=verbose,
                                         response=response,
                                         ...),
                    
                    
                    stop("Method ", control$main.method, " is not implemented.")
  )
  
  initialfit <- NULL
  
  # done with main fit
  
  if(!is.null(control$MCMLE.check.degeneracy) && control$MCMLE.check.degeneracy && (is.null(mainfit$theta1$independent) || !all(mainfit$theta1$independent))){
    if(verbose) {
      message("Checking for degeneracy.")
    }
    degeneracy <- ergm.degeneracy(mainfit, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  mainfit$ergm_version <- packageVersion("ergm")
  mainfit$MPLE_is_MLE <- MPLE.is.MLE
  mainfit$degeneracy.value <- degeneracy$degeneracy.value
  mainfit$degeneracy.type <- degeneracy$degeneracy.type
  
  mainfit$formula <- formula
  mainfit$target.stats <- model$target.stats
  mainfit$target.esteq <- if(!is.null(model$target.stats)) ergm.estfun(rbind(model$target.stats), mainfit$coef, model)
  
  mainfit$constrained <- proposal$arguments$constraints
  mainfit$constrained.obs <- proposal.obs$arguments$constraints
  mainfit$constraints <- constraints
  
  # unless the main fitting algorithm passes back a modified control
  if (is.null(mainfit$control)) mainfit$control<-control
  
  mainfit$response<-response
  mainfit$reference<-reference
  mainfit$estimate <- estimate
  
  mainfit$offset <- model$etamap$offsettheta
  mainfit$drop <- if(control$drop) extremecheck$extremeval.theta
  mainfit$estimable <- constrcheck$estimable
  mainfit$etamap <- model$etamap
  
  mainfit$null.lik<-logLikNull.ergm(mainfit, verbose=verbose)
  
  if (!control$MCMC.return.stats)
    mainfit$sample <- NULL
  
  if(eval.loglik){
    message("Evaluating log-likelihood at the estimate. ", appendLF=FALSE)
    mainfit<-logLik(mainfit, add=TRUE, control=control$loglik.control, verbose=verbose)
  }
    
  if (MCMCflag) {
    message(paste(strwrap("This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function."),collapse="\n"))
  }
  
  mainfit
}
