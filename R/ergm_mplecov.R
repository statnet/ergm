#' Approximate MPLE standard errors in dyad-dependent models
#'
#' Function to approximate the MPLE covariance matrix in a dyad dependence model
#' using the Godambe matrix as described in Schmid and Hunter (2020) or by parametric bootstrap
#' as described by Schmid and Desmarais (2017).
#'
#' @param pl An \code{\link{ergm.pl}} object.
#' @param nw response network.
#' @param m the model, as returned by \code{\link{ergm_model}}
#' @param init a vector a vector of initial theta coefficients
#' @param theta.mple the MPLE of a given model
#' @param invHess the inverse Hessian matrix obtained from glm()
#' @param control a list of MCMC related parameters; recognized
#'   components include: samplesize : the number of networks to sample
#'   Clist.miss : see 'Clist.miss' above; some of the code uses this
#'   Clist.miss,
#' @param verbose whether this and the C routines should be verbose (T
#'   or F); default=FALSE
#'
#' @param constraints {A formula specifying one or more constraints
#' on the support of the distribution of the networks being modeled,
#' using syntax similar to the \code{formula} argument, on the
#' right-hand side.
#'
#' @param family the family to use in the R native routine
#' default="binomial"
#'
#' @return \code{ergm_mplecov} returns a list either
#'   containing a Godambe covariance matrix or a diagonal matrix with bootstrap variances.
#'
#' @references Schmid CS and Desmarais BA (2017). "Exponential random graph
#' models with big networks: Maximum pseudolikelihood estimation and the parametric bootstrap"
#' _IEEE International Conference on Big Data (Big Data)_, pp. 116-121.
#'
#' Schmid CS and Hunter DR (2023).  "Computing Pseudolikelihood Estimators for Exponential-Family Random Graph Models" _Journal of Data Science_.
#' @noRd
#'
#' @examples
#' \donttest{
#' # create initial network with x=50 nodes
#' init.net <- network(x=50, directed = FALSE, density = 0.1)
#'
#' # true parameters for edges, kstar(2), ad triangles
#' truth <- c(-0.25, -0.2, 0.5)
#'
#' # simulate a network from the ERGM defined by the true parameter values
#' sim.net <- simulate(init.net~ edges+kstar(2)+ triangles, nsim=1, coef=truth)
#'
#' # get MPLE with inverse Hessian matrix
#' fisher <- ergm(sim.net ~ edges+kstar(2)+ triangles, estimate = "MPLE",
#' control = control.ergm(MPLE.covariance.method = "invHess"))
#'
#' # get MPLE with Godambe matrix
#' godambe <- ergm(sim.net ~ edges+kstar(2)+ triangles, estimate = "MPLE",
#' control = control.ergm(MPLE.covariance.method = "Godambe"))
#'
#' # get MPLE with bootstrap standard errors
#' bootstrap <- ergm(sim.net ~ edges+kstar(2)+ triangles, estimate="MPLE",
#' control = control.ergm(MPLE.covariance.method="bootstrap"))
#' }
ergm_mplecov <- function(pl,
                         nw,
                         s,
                         init=init,
                         theta.mple,
                         invHess,
                         control=NULL,
                         verbose=FALSE,
                         constraints=NULL,
                         family="binomial",
                         formula = formula){

  m <- s$model
  # get sample size from control.ergm
  R <- control$MPLE.covariance.samplesize
  mple.burnin <- control$MPLE.covariance.sim.burnin
  mple.interval <- control$MPLE.covariance.sim.interval

  # Simulate R networks
  sim.mple <- simulate(m, basis=nw, coef=theta.mple, nsim=R,
                       control=control.simulate.formula(MCMC.burnin=mple.burnin, MCMC.interval=mple.interval))

  num.variables <- ncol(pl$xmat)

  terms.form <- terms(formula)

  new.formula <- reformulate(attr(terms.form, "term.labels"), "sim.mple[[i]]")

  if(control$MPLE.covariance.method == "Godambe"){
    message("Estimating Godambe Matrix using ", R, " simulated networks.")

    # calculation of V(theta) = Var(u(theta,y)) using the sim.num networks
    net.stat <- matrix(0, nrow=length(sim.mple), ncol=num.variables)
    colnames(net.stat) <- colnames(pl$xmat)
    u.data <- matrix(0,nrow=length(sim.mple), ncol=num.variables)

    for(i in 1:length(sim.mple)){

      dat <- ergmMPLE(new.formula)
      net.stat[i,] <- summary(new.formula)

      # write the response, weight and designmatrix into one matrix
      X.dat <- cbind(dat$response, dat$weights, dat$predictor)

      # calculate s(theta)
      u.data[i,] <- sapply(1:num.variables, function(k){
        sum(apply(X.dat, 1, function(x){
          x[2]*(x[k+2]*(x[1] - exp(sum(theta.mple*x[(3:(num.variables+2))]))/(1+exp(sum(theta.mple*x[(3:(num.variables+2))]))) )) } ) ) } )

    } # end for i
    # calculate V.hat by estimating sd
    u.mean <- colMeans(u.data)
    u.sum <- matrix(0,num.variables,num.variables)
    for(i in 1: nrow(u.data)){
      u.diff <- u.data[i,]- u.mean
      u.sum <- u.sum + u.diff%*%t(u.diff)
    }
    u.sum.n <- u.sum/(nrow(u.data)-1)
    G <- invHess%*%u.sum.n%*%invHess
    return(G)

  } # end if Godambe

  if(control$MPLE.covariance.method == "bootstrap"){
    message("Estimating Bootstrap Standard Errors using ", R, " simulated networks.")

    # create empty matrix to store mple of bootstrap samples
    boot.mple.mat <- matrix(0, nrow=length(sim.mple), ncol=num.variables)
    colnames(boot.mple.mat) <- colnames(pl$xmat)

    for(i in 1:length(sim.mple)){

      dat <- ergmMPLE(new.formula)

      # calculate MPLE of simulated network
      glm.sim <- glm(dat$response ~ .-1 , data=data.frame(dat$predictor),
                     weights=dat$weights, family="binomial")
      boot.mple.mat[i,] <- coef(glm.sim)

    }# end for i
    Boot.cov <- matrix(0,num.variables, num.variables)
    diag(Boot.cov) <- apply(boot.mple.mat, 2, var)
    return(Boot.cov)

  } # end if bootstrap
} # end function