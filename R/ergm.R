#  File R/ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
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
#' @param formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a matrix that can be
#'   coerced to a \code{\link[network]{network}} object.  For the
#'   details on the possible \code{<model terms>}, see
#'   \code{\link{ergmTerm}} and Morris, Handcock and Hunter (2008)
#'   for binary ERGM terms and Krivitsky (2012) for valued ERGM terms
#'   (terms for weighted edges).  To create a
#'   \code{\link[network]{network}} object in \R, use the
#'   \code{network()} function, then add nodal attributes to it using
#'   the \code{\%v\%} operator if necessary. Enclosing a model term in
#'   \code{offset()} fixes its value to one specified in
#'   \code{offset.coef}.  (A second argument---a logical or numeric
#'   index vector---can be used to select *which* of the parameters
#'   within the term are offsets.)
#'
#' @template response
#' @template reference
#' @param constraints {A formula specifying one or more constraints
#' on the support of the distribution of the networks being modeled,
#' using syntax similar to the \code{formula} argument, on the
#' right-hand side. Multiple constraints
#' may be given, separated by \dQuote{+} and \dQuote{-} operators. (See
#' [ERGM constraints][ergmConstraint] for the explanation of
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
#' See the [ERGM constraints][ergmConstraint] documentation for
#' the constraints implemented in the **[ergm][ergm-package]**
#' package. Other packages may add their own constraints.
#' 
#' Note that not all possible combinations of constraints and reference
#' measures are supported. However, for relatively simple constraints
#' (i.e., those that simply permit or forbid specific dyads or sets of
#' dyads from changing), arbitrary combinations should be possible.
#' }
#'
#' @param obs.constraints { A one-sided formula specifying one or more
#'   constraints or other modification \emph{in addition} to those
#'   specified by \code{constraints} that had affected the observation
#'   process for the network, using syntax similar to the
#'   \code{formula} argument. Multiple constraints may be given,
#'   separated by \dQuote{+} operators.
#' 
#'   This allows the domain of the integral in the numerator of the
#'   partially obseved network face-value likelihoods of Handcock and
#'   Gile (2010) and Karwa et al. (2017) to be specified explicitly.
#' 
#'   The default is \code{~observed}, to constrain the integral to
#'   only integrate over the missing dyads. (It is dropped
#'   automatically if the network is completely observed.)
#'     
#'   It is also possible to specify a proposal function directly by
#'   passing a string with the function's name. In that case,
#'   arguments to the proposal should be specified through the
#'   \code{obs.prop.args} argument to \code{\link{control.ergm}}.
#' 
#'   See the [ERGM constraints][ergmConstraint] documentation for
#'   the constraints implemented in the **[ergm][ergm-package]**
#'   package. Other packages may add their own constraints.
#'     
#'   Note that not all possible combinations of constraints and
#'   reference measures are supported.
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
#' values like AIC and BIC) are not needed. Can be set globally via `option(ergm.eval.loglik=...)`, which is set to `TRUE` when the package is loaded. (See [`options?ergm`][ergm-options].)
#' }
#' @param estimate {If "MPLE," then the maximum pseudolikelihood estimator
#' is returned.  If "MLE" (the default), then an approximate maximum likelihood
#' estimator is returned.  For certain models, the MPLE and MLE are equivalent,
#' in which case this argument is ignored.  (To force MCMC-based approximate
#' likelihood calculation even when the MLE and MPLE are the same, see the
#' \code{force.main} argument of \code{\link{control.ergm}}. If "CD" (\emph{EXPERIMENTAL}),
#' the Monte-Carlo contrastive divergence estimate is returned. )
#' }
#'
#' @templateVar mycontrol control.ergm
#' @template control
#' @template verbose
#'
#' @param \dots Additional
#' arguments, to be passed to lower-level functions.
#'
#' @template basis
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
#' \item{network}{Network passed on the left-hand side of `formula`. If `target.stats` are passed, it is replaced by the network returned by [san()].}
#' \item{newnetworks}{A list of the final networks at the end of the MCMC
#' simulation, one for each thread.}
#' \item{newnetwork}{The first (possibly only) element of \code{newnetworks}.}
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
#' Version 1. \url{https://statnet.org}.
#' 
#' Bender-deMoll S, Morris M, Moody J (2008).
#' Prototype Packages for Managing and Animating Longitudinal
#' Network Data: \pkg{dynamicnetwork} and \pkg{rSoNIA}.
#' \emph{Journal of Statistical Software}, 24(7).
#' \doi{10.18637/jss.v024.i07}
#' 
#' 
#' Butts CT (2007).
#' \pkg{sna}: Tools for Social Network Analysis.
#' R package version 2.3-2. \url{https://cran.r-project.org/package=sna}.
#' 
#' Butts CT (2008).
#' \pkg{network}: A Package for Managing Relational Data in \R.
#' \emph{Journal of Statistical Software}, 24(2).
#' \doi{10.18637/jss.v024.i02}
#' 
#' Butts C (2015).
#' \pkg{network}: The Statnet Project (https://statnet.org). R package version 1.12.0, \url{https://cran.r-project.org/package=network}.
#' 
#' Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
#' A \pkg{statnet} Tutorial.
#' \emph{Journal of Statistical Software}, 24(8).
#' \doi{10.18637/jss.v024.i08}
#' 
#' Goodreau SM, Kitts J, Morris M (2008b).
#' Birds of a Feather, or Friend of a Friend? Using Exponential
#' Random Graph Models to Investigate Adolescent Social Networks.
#' \emph{Demography}, 45, in press.
#' 
#' Handcock, M. S. (2003)
#' \emph{Assessing Degeneracy in Statistical Models of Social Networks},
#' Working Paper #39,
#' Center for Statistics and the Social Sciences,
#' University of Washington.
#' \url{https://csss.uw.edu/research/working-papers/assessing-degeneracy-statistical-models-social-networks}
#' 
#' Handcock MS (2003b).
#' \pkg{degreenet}: Models for Skewed Count Distributions Relevant
#' to Networks.
#' Statnet Project, Seattle, WA.
#' Version 1.0, \url{https://statnet.org}.
#' 
#' Handcock MS and Gile KJ (2010). Modeling Social Networks from Sampled Data. \emph{Annals of Applied Statistics}, 4(1), 5-25. \doi{10.1214/08-AOAS221}
#' 
#' Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003a).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks.
#' Statnet Project, Seattle, WA.
#' Version 2, \url{https://statnet.org}.
#' 
#' Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
#' \pkg{statnet}: Software Tools for the Statistical Modeling of
#' Network Data.
#' Statnet Project, Seattle, WA.
#' Version 2, \url{https://statnet.org}.
#' 
#' Hunter, D. R. and Handcock, M. S. (2006)
#' \emph{Inference in curved exponential family models for networks},
#' Journal of Computational and Graphical Statistics.
#' 
#' Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks.
#' \emph{Journal of Statistical Software}, 24(3).
#' \doi{10.18637/jss.v024.i03}
#'
#' Karwa V, Krivitsky PN, and Slavkovi\'{c} AB (2017). Sharing Social Network
#' Data: Differentially Private Estimation of Exponential-Family Random
#' Graph Models. \emph{Journal of the Royal Statistical Society, Series
#' C}, 66(3):481--500. \doi{10.1111/rssc.12185}
#' 
#' Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#' 
#' Morris M, Handcock MS, Hunter DR (2008).
#' Specification of Exponential-Family Random Graph Models:
#' Terms and Computational Aspects.
#' \emph{Journal of Statistical Software}, 24(4).
#' \doi{10.18637/jss.v024.i04}
#' 
#' Snijders, T.A.B. (2002),
#' Markov Chain Monte Carlo Estimation of Exponential Random Graph Models.
#' Journal of Social Structure.
#' Available from 
#' \url{https://www.cmu.edu/joss/content/articles/volume3/Snijders.pdf}.
#' 
#' @seealso [`network`], [`%v%`], [`%n%`], [`ergmTerm`], [`ergmMPLE`],
#' [summary.ergm()]
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
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"))
#' summary(gest)
#' 
#' # compare it to differential homophily by atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle
#'                         + nodematch("atomic type",diff=TRUE))
#' summary(gest)
#' }
#' @keywords models
#' @aliases is.ergm ergm.object
#' @export
ergm <- function(formula, response=NULL,
                 reference=~Bernoulli,
                 constraints=~.,
                 obs.constraints=~.-observed,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=getOption("ergm.eval.loglik"),
                 estimate=c("MLE", "MPLE", "CD"),
                 control=control.ergm(),
                 verbose=FALSE,..., basis=ergm.getnetwork(formula)) {
  check.control.class("ergm", "ergm")
  handle.control.toplevel("ergm", ...)

  ergm_call <- match.call(ergm)

  reference <- trim_env_const_formula(reference)
  constraints <- trim_env_const_formula(constraints)
  obs.constraints <- trim_env_const_formula(obs.constraints)
  
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
  
  nw <- basis
  ergm_preprocess_response(nw,response)

  proposalclass <- "c"

  if(!is(constraints, "ergm_proposal")){
    # Handle the observation process and other "automatic" constraints.
    tmp <- .handle.auto.constraints(nw, constraints, obs.constraints, target.stats)
    nw <- tmp$nw
    constraints.obs <- tmp$constraints.obs
    constraints <- tmp$constraints
  }else if(!is(obs.constraints, "ergm_proposal")){
  # Handle the observation process and other "automatic" constraints.
    tmp <- .handle.auto.constraints(nw, trim_env(~.), obs.constraints, target.stats)
    nw <- tmp$nw
    constraints.obs <- tmp$constraints.obs
    constraints <- tmp$constraints
  }
  
  if(!is(constraints, "ergm_proposal")){
    if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
    
  ## FIXME: a more general framework is needed?
  if(is.valued(nw) && reference==trim_env(~Bernoulli)){
    warn(paste0("The default Bernoulli reference distribution operates in the binary (",sQuote("response=NULL"),") mode only. Did you specify the ",sQuote("reference")," argument?"))
  }
    
    proposal <- ergm_proposal(constraints, hints=control$MCMC.prop, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass,reference=reference, term.options=control$term.options)
  }else proposal <- constraints
  
  if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")
  
  if (verbose) message("Initializing model...")
  model <- ergm_model(formula, nw, extra.aux=NVL3(proposal$auxiliaries,list(proposal=.)), term.options=control$term.options)
  proposal$aux.slots <- model$slots.extra.aux$proposal
  if (verbose) message("Model initialized.")
  
  if(!is(obs.constraints, "ergm_proposal")){
    if(!is.null(constraints.obs)){
      if (verbose) message("Initializing constrained Metropolis-Hastings proposal: ", appendLF=FALSE)
      proposal.obs <- ergm_proposal(constraints.obs, hints=control$obs.MCMC.prop, weights=control$obs.MCMC.prop.weights, control$obs.MCMC.prop.args, nw, class=proposalclass, reference=reference, term.options=control$term.options)
      if (verbose) message(sQuote(paste0(proposal.obs$pkgname,":MH_",proposal.obs$name)), appendLF=FALSE)
      
      if(!is.null(proposal.obs$auxiliaries)){
        if(verbose) message(" (requests auxiliaries: updating model).")
        model$obs.model <- c(model, ergm_model(trim_env(~.), nw, extra.aux=list(proposal=proposal.obs$auxiliaries), term.options=control$term.options))
        proposal.obs$slots.extra.aux <- model$model.obs$slots.extra.aux$proposal
        if(verbose) message("Model reinitialized.")
      }else if(verbose) message(".")
    }else proposal.obs <- NULL
  }else proposal.obs <- obs.constraints
  
  ## Construct approximate response network if target.stats are given.
  if(!is.null(target.stats)){
    nw.stats <- summary(model, nw)[!model$etamap$offsetmap]
    target.stats <- vector.namesmatch(target.stats, names(nw.stats))
    target.stats <- na.omit(target.stats)
    if(length(nw.stats)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    
    # no need to pass the offset term's init to SAN
    san.control <- control$SAN
    
    if(verbose) message("Constructing an approximate response network.")
    ## If target.stats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    if(control$SAN.maxit > 0){
      TARGET_STATS<-san(formula, target.stats=target.stats,
                reference=reference,
                constraints=constraints,
                control=san.control,
                only.last=TRUE,
                output="ergm_state",
                verbose=verbose,
                basis=nw,
                offset.coef=NVL(offset.coef,control$init[model$etamap$offsettheta]))
      if(verbose) message("Finished SAN run.")
    }else{
      TARGET_STATS <- nw
    }
    
    # From this point on, target.stats has NAs corresponding to the
    # offset terms.
    #
    # TODO: Only have target.stats contain non-offset terms'
    # statistics, and have the rest of the code handle it
    # intelligently.
    target.stats <- .align.target.stats.offset(model, target.stats)   

    nw <- TARGET_STATS <- as.network(TARGET_STATS)
    formula<-nonsimp_update.formula(formula,TARGET_STATS~., from.new="TARGET_STATS")
  } else {
    if (network.edgecount(nw) == 0) warning("Network is empty and no target stats are specified.")
  }

  # TODO: SAN has this information, so maybe we should grab it from there if SAN does get run.
  nw.stats <- summary(model, nw, term.options=control$term.options)

  # conddeg MPLE has been superceded, but let the user know:
  if(!is.directed(nw) && ("degrees" %in% names(proposal$arguments$constraints) ||
                                           all(c("b1degrees","b2degrees") %in% names(proposal$arguments$constraints)))) message("Note that degree-conditional MPLE has been removed in version 4.0, having been superceded by Contrastive Divergence.")  
  
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

  control$init.method <- ERRVL(try(match.arg(control$init.method, init.candidates), silent=TRUE), {
    message("Specified initial parameter method ", sQuote(control$init.method), " is not in the list of candidates. Use at your own risk.")
    control$init.method
  })
  if(verbose) message(paste0("Using initial method '",control$init.method,"'."))

  if(is.curved(model) ||
     !all(model$etamap$mintheta==-Inf) ||
     !all(model$etamap$maxtheta==+Inf)){ # Curved or constrained model: use ergm.logitreg() rather than glm().
      control$MPLE.type <- "logitreg"
  }
  
  # If some control$init is specified...
  if(!is.null(control$init)){
    # Check length of control$init.
    if (length(control$init)!=length(model$etamap$offsettheta)) {
      if(verbose){
        message("control$init =")
        message_print(control$init)
        message("number of statistics is ", nparam(model, canonical=TRUE), "")
      }
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters."))
    }
  }else control$init <- rep(NA, length(model$etamap$offsettheta)) # Set the default value of control$init.
  
  if(!is.null(offset.coef)){
      # TODO: Names matching here?
      if(length(control$init[model$etamap$offsettheta])!=length(offset.coef))
          stop("Invalid offset parameter vector offset.coef: ",
               "wrong number of parameters: expected ",
               length(control$init[model$etamap$offsettheta]),
               " got ",length(offset.coef),".")
      control$init[model$etamap$offsettheta]<-offset.coef
  }
  
  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(param_names(model)[is.na(control$init)&model$offsettheta]), ".", sep="")
  
  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, proposal, control$init)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=control$init, target.stats=target.stats, drop=control$drop)
  model <- extremecheck$model; control$init <- extremecheck$init
  
  # MPLE is not supported for valued ERGMs.
  if(estimate=="MPLE"){
    if(is.valued(nw)) stop("Maximum Pseudo-Likelihood (MPLE) estimation for valued ERGM terms is not implemented at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
    if(!is.dyad.independent(proposal$arguments$constraints,
                            proposal.obs$arguments$constraints))
      message("Maximum Pseudo-Likelihood (MPLE) estimation for ERGMs with dyad-dependent constraints is not implemented at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
  }
  
  if (verbose) { message("Fitting initial model.") }
  
  MPLE.is.MLE <- (proposal$reference$name=="Bernoulli"
                  && is.dyad.independent(model)
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
  if(all(model$etamap$offsettheta)){
    # Note that this cannot be overridden with control$force.main.
    message("All terms are either offsets or extreme values. No optimization is performed.")
    return(structure(list(coefficients=control$init,
                          call=ergm_call,
                          iterations=0,
                          loglikelihood=NA,
                          mle.lik=NULL,
                          gradient=rep(NA,length=length(control$init)),
                          failure=TRUE,
                          offset=model$etamap$offsettheta,
                          drop=if(control$drop) extremecheck$extremeval.theta,
                          estimable=constrcheck$estimable,
                          network=nw,
                          reference=reference,
                          newnetwork=nw,
                          formula=formula,
                          constrained=proposal$arguments$constraints,
                          constrained.obs=proposal.obs$arguments$constraints,
                          constraints=constraints,
                          target.stats=target.stats,
                          target.esteq=if(!is.null(target.stats)) ergm.estfun(rbind(target.stats), control$init, model),
                          estimate=estimate,
                          ergm_version=packageVersion("ergm"),
                          control=control
    ),
    class="ergm"))
    
  }
  
  model$nw.stats <- nw.stats
  model$target.stats <- target.stats
  
  if(control$init.method=="CD") if(is.null(names(control$init)))
      names(control$init) <- param_names(model, FALSE)
  
  initialfit <- ergm.initialfit(init=control$init, initial.is.final=!MCMCflag,
                                formula=formula, nw=nw, reference=reference, 
                                m=model, method=control$init.method,
                                MPLEtype=control$MPLE.type, 
                                control=control,
                                proposal=proposal,
                                proposal.obs=proposal.obs,
                                verbose=if(MCMCflag) FALSE else verbose,
                                ...)

  switch(control$init.method,
         MPLE = NVL3(initialfit$xmat.full, check_nonidentifiability(., coef(initialfit), model,
                                                                    tol = control$MPLE.nonident.tol, type="covariates",
                                                                    nonident_action = control$MPLE.nonident,
                                                                    nonvar_action = control$MPLE.nonvar)),
         CD = NVL3(initialfit$sample, check_nonidentifiability(as.matrix(.), coef(initialfit), model,
                                                               tol = control$MPLE.nonident.tol, type="statistics",
                                                               nonident_action = control$MPLE.nonident,
                                                               nonvar_action = control$MPLE.nonvar))
         )

  initialfit$xmat.full <- NULL # No longer needed but takes up space.

  estimate.desc <- switch(estimate,
                          MPLE = if(MPLE.is.MLE) "Maximum Likelihood"
                                 else "Maximum Pseudolikelihood",
                          CD = "Contrastive Divergence",
                          MLE = paste0(if(MCMCflag) # If not, it's just MLE.
                                         switch(control$main.method,
                                                MCMLE = "Monte Carlo ",
                                                `Stochastic-Approximation`="Stochastic Approximation ",
                                                `Robbins-Monro`="Robbins-Monro ",
                                                `Stepping`="Hummel Stepping "),
                                       "Maximum Likelihood"))

  if (!MCMCflag){ # Just return initial (non-MLE) fit and exit.
    message("Stopping at the initial estimate.")
    initialfit$call <- ergm_call
    initialfit$ergm_version <- packageVersion("ergm")
    initialfit$offset <- model$etamap$offsettheta
    initialfit$MPLE_is_MLE <- MPLE.is.MLE
    initialfit$drop <- if(control$drop) extremecheck$extremeval.theta
    initialfit$estimable <- constrcheck$estimable
    initialfit$network <- nw
    initialfit$reference <- reference
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constrained <- proposal$arguments$constraints
    initialfit$constrained.obs <- proposal.obs$arguments$constraints
    initialfit$constraints <- constraints
    initialfit$obs.constraints <- obs.constraints 
    initialfit$target.stats <- suppressWarnings(na.omit(model$target.stats))
    initialfit$nw.stats <- model$nw.stats
      initialfit$etamap <- model$etamap
    initialfit$target.esteq <- suppressWarnings(na.omit(if(!is.null(model$target.stats)) ergm.estfun(rbind(model$target.stats), coef(initialfit), model)))
    initialfit$estimate <- estimate
    initialfit$estimate.desc <- estimate.desc

    initialfit$control<-control

    if(any(!model$etamap$offsettheta) && eval.loglik){
      message("Evaluating log-likelihood at the estimate. ",appendLF=FALSE)
      initialfit<-logLik(initialfit, add=TRUE, control=control$loglik, verbose=verbose)
      message("")
    }
    return(initialfit)
  }
  
  # Otherwise, set up the main phase of estimation:
  
  ergm.getCluster(control, max(verbose-1,0))
  
  # Revise the initial value, if necessary:
  init <- coef(initialfit)
  init[is.na(init)] <- 0
  names(init) <- param_names(model, FALSE)
  
  mainfit <- switch(control$main.method,
                    "Robbins-Monro" = ergm.robmon(init, nw, model, 
                                                  proposal=proposal, verbose=verbose, control=control),
                    "Stochastic-Approximation" = ergm.stocapprox(init, nw, model,
                                                                 control=control, proposal=proposal,
                                                                 verbose=verbose),
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
                                         ...),
                    
                    
                    stop("Method ", control$main.method, " is not implemented.")
  )
  
  initialfit <- NULL
  
  # done with main fit
  
  mainfit$call <- ergm_call
  mainfit$ergm_version <- packageVersion("ergm")
  mainfit$MPLE_is_MLE <- MPLE.is.MLE
  
  mainfit$formula <- formula
  mainfit$target.stats <- suppressWarnings(na.omit(model$target.stats))
  mainfit$nw.stats <- model$nw.stats
  mainfit$target.esteq <- suppressWarnings(na.omit(if(!is.null(model$target.stats)) ergm.estfun(rbind(model$target.stats), coef(mainfit), model)))
  
  mainfit$constrained <- proposal$arguments$constraints
  mainfit$constrained.obs <- proposal.obs$arguments$constraints
  mainfit$constraints <- constraints
  mainfit$obs.constraints <- obs.constraints
  
  # unless the main fitting algorithm passes back a modified control
  NVL(mainfit$control) <- control
  
  mainfit$reference<-reference
  mainfit$estimate <- estimate
  mainfit$estimate.desc <- estimate.desc

  mainfit$offset <- model$etamap$offsettheta
  mainfit$drop <- if(control$drop) extremecheck$extremeval.theta
  mainfit$estimable <- constrcheck$estimable
  mainfit$etamap <- model$etamap

  if (!control$MCMC.return.stats)
    mainfit$sample <- NULL

  if(eval.loglik){
    message("Evaluating log-likelihood at the estimate. ", appendLF=FALSE)
    mainfit<-logLik(mainfit, add=TRUE, control=control$loglik, verbose=verbose)
  }
    
  if (MCMCflag) {
    message(paste(strwrap("This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function."),collapse="\n"))
  }
  
  mainfit
}
