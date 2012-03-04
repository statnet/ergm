################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
#
# --PARAMETERS--
#   formation   : the formation formula, as 'nw ~ term(s)'
#   dissolution : the dissolution formula, as 'nw ~ term(s)'
#   theta.form0 : the intial theta formation coefficients, or optionally if
#                 these are to be estimates, the string "MPLE";
#                 default="MPLE"
#   theta.diss  : the initial theta dissolution coefficients
#   seed        : an integer starting value for the random number generator;
#                 default=NULL
#   MH.burnin   : the number of proposals used in each MCMC step; this is ignored
#                 unless 'control$main.method'="Robbins-Monro"; any other style or
#                 the default style will not recognize this parameter;
#                 default=1000
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     indegreedist
#                      observed  outdegreedist
#                default="~ ."; these may not work currently.
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   control     :  a list of control parameters returned from <control.stergm>;
#                  default=<control.stergm>()
#   verbose     :  whether ergm should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list 
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:
#       <stergm>             = $
#       <stergm.RM>          = @
#       <stergm.NM>          = !
#
#   the components include:
#
#     @   coef.form   : the estimated formation model coefficients
#     @   coef.diss   : the estimated dissolution model coefficients
#      &  eta         : the estimated formation ?? coefficients
#   $     offset      : a logical vector whose ith entry tells whether the
#                        ith curved theta coeffient was offset/fixed
#   $     etamap      :  the list constituting the theta->eta mapping for the
#                        formation model; for details of its components,
#                        see <ergm.etamap>
#   $     MH.burnin   :  the number of proposals made in each MCMC step
#   $     formation   : the formation formula, as 'nw ~ term(s)'
#   $     dissolution : the dissolution formula, as 'nw ~ term(s)'
#   $     constraints : the constraints formula
#     @&  newnetwork  :  the final network sampled
#     @&  network    :  the 'nw' inputted to <ergm> via the 'formula'
#   $     prop.args.form     :  the MHP formation arguments passed to the
#                               InitMHP rountines
#   $     prop.args.diss     :  the MHP dissolution arguments passed to the
#                               InitMHP rountines
#   $     prop.weights.form  :  the method used to allocate probabilities of
#                               being proposed to dyads in the formation stage,
#                               as "TNT", "random", "nonobserved", or "default"
#      &  theta.original     :  the theta values at the start of the MCMC 
#                               sampling
#     @   theta.form.original:  the formation theta values at the start of the
#                               MCMC sampling
#   $     prop.weights.diss  :  as 'prop.weights.form', but for the dissolution
#                               model
#
################################################################################

stergm.CMLE <- function(nw, formation, dissolution, times, offset.coef.form, offset.coef.diss,
                        eval.loglik,
                        estimate,
                        control,
                        verbose) {

  # Translate the "estimate" from the stergm() argument to the ergm() argument.
  estimate <- switch(estimate,
                     CMLE = "MLE",
                     CMPLE = "MPLE")

  if(is.null(times)){
    if(inherits(nw, "network.list") || is.list(nw)){
      times  <- c(1,2)
      warning("Time points not specified for a list. Modeling transition from the first to the second network. This behavior may change in the future.")
    }else if(inherits(nw,"networkDynamic")){
      times  <- c(0,1)
      warning("Time points not specified for a networkDynamic. Modeling transition from time 0 to 1.")
    }
  }
  
  if(length(times)<2) stop("Time points whose transition is to be modeled was not specified.")
  if(length(times)>2) stop("Only two time points (one transition) are supported at this time.")

  if(inherits(nw, "network.list") || is.list(nw)){
    y0 <- nw[[times[1]]]
    y1 <- nw[[times[2]]]
  }else if(inherits(nw,"networkDynamic")){
    require(networkDynamic) # This is needed for the "%t%.network" function
    y0 <- nw %t% times[1]
    y1 <- nw %t% times[2]
  }
  
  
  # Construct the formation and dissolution networks; the
  # network.update is there to copy attributes from y0 to y.form and
  # y.diss.
  
  y.form <- y0 | y1
  y.form <- network.update(y0, as.matrix(y.form, matrix.type="edgelist"), matrix.type="edgelist")
  formation <- ergm.update.formula(formation, y.form~.)

  y.diss <- y0 & y1
  y.diss <- network.update(y0, as.matrix(y.diss, matrix.type="edgelist"), matrix.type="edgelist")
  dissolution <- ergm.update.formula(dissolution, y.diss~.)

  # Apply initial values passed to control.stergm() the separate controls, if necessary.
  if(is.null(control$CMLE.control.form$init)) control$CMLE.control.form$init <- control$init.form
  if(is.null(control$CMLE.control.diss$init)) control$CMLE.control.diss$init <- control$init.diss
  
  # Now, call the ergm()s:
  
  fit.form <- ergm(formation, constraints=~atleast(y0), offset.coef=offset.coef.form, eval.loglik=eval.loglik, estimate=estimate, control=control$CMLE.control.form, verbose=verbose)

  fit.diss <- ergm(dissolution, constraints=~atmost(y0), offset.coef=offset.coef.diss, eval.loglik=eval.loglik, estimate=estimate, control=control$CMLE.control.diss, verbose=verbose)

  # Construct the output list. Conveniently, this is mainly a list consisting of two ergms.
  
  list(network=nw, times=times, formation=formation, dissolution=dissolution, formation.fit=fit.form, dissolution.fit=fit.diss, estimate=estimate)
}
