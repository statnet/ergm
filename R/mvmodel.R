mvmodel <- function(object, ...)
UseMethod("mvmodel")

mvmodel.default <- function(object,...)
{
  stop("Either a network, an ergm object or a formula argument must be given")
}

mvmodel.formula <- function (formula, ..., theta0, nsim=100,
                             burnin=100, interval=100,
                             constraints=~., prop.weights="default", prop.args=NULL,
                             seed=NULL,  drop=FALSE,
                             statistic=NULL
		      ) {
  trms <- ergm.getterms(formula)
  if(length(trms)>2){
    g <- eval(trms[[2]], sys.parent())
  }else{
     stop("A network object on the RHS of the formula must be given")
  }
  g <- as.network(g)
  if(!is.network(g)){
    stop("A network object on the RHS of the formula must be given")
  }

  m <- ergm.getmodel(formula, g, drop=drop)
  Clist <- ergm.Cprepare(g, m)

  if(missing(theta0)){
      theta0 <- rep(0,Clist$nparam)
      warning("No parameter values given, using 0\n\t")
  }

  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}

  probabilites <- FALSE
  if(!is.function(statistic)){
   if(is.character(statistic) && statistic=="density"){
    statistic <- function(x){summary(x ~ density)}
   }else{
    statistic <- function(x){as.sociomatrix(x)}
    probabilites <- TRUE
   }
  }

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(formula, burnin=burnin, interval=interval,
                                constraints=constraints,
                                control=control.simulate.ergm(prop.args=prop.args,
                                  prop.weights=prop.weights,
                                  drop=drop),
                                theta0=theta0,
                                n=nsim, seed=seed)
  
# cat("\nCollating simulations\n")

  # Set up the output arrays

  simcentrality <- statistic(SimGraphSeriesObj$networks[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj$networks[[i]])
    if(!probabilites){
      sim.mvmodel[i,] <- simcentrality
    }else{
      sim.mvmodel <- sim.mvmodel + simcentrality
    }
  }

  if(!probabilites){
    simcentrality <- apply(sim.mvmodel,2,mean)
    return(list(mean=simcentrality,
                   sim=sim.mvmodel))
  }else{
    return(sim.mvmodel / nsim)
  }
}

mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=100, interval=100,
                          constraints=NULL, prop.weights="default", prop.args=NULL,
                          seed=NULL, drop=FALSE,
                          statistic=NULL) {

# trms <- ergm.getterms(object$formula)
# g <- as.network(eval(trms[[2]], sys.parent()))
  g <- object$network
  g <- as.network(g)
  if(!is.network(g)){
    stop("A network object must be given")
  }
  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  constraints <-
    if(is.null(constraints)) object
    else constraints

  probabilites <- FALSE
  if(!is.function(statistic)){
   if(is.character(statistic) && statistic=="density"){
    statistic <- function(x){summary(x ~ density)}
   }else{
    statistic <- function(x){as.sociomatrix(x)}
    probabilites <- TRUE
   }
  }

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(object, burnin=burnin, interval=interval,
                                constraints=constraints,
                                control=control.simulate.ergm(prop.args=prop.args,
                                  prop.weights=prop.weights,
                                  drop=drop),
                             n=nsim, seed=seed)

  # cat("\nCollating simulations\n")

  # Set up the output arrays

  simcentrality <- statistic(SimGraphSeriesObj$networks[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj$networks[[i]])
    if(!probabilites){
      sim.mvmodel[i,] <- simcentrality
    }else{
      sim.mvmodel <- sim.mvmodel + simcentrality
    }
  }

  if(!probabilites){
    simcentrality <- apply(sim.mvmodel,2,mean)
    return(list(mean=simcentrality,
                   sim=sim.mvmodel))
  }else{
    return(sim.mvmodel / nsim)
  }
}
