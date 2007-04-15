simulatedyn <- function(object, nsim=1, seed=NULL, ...,theta0,
                        burnin=1, interval=1,
                        sequential=TRUE,
                        proposaltype="formationTNT",
                        proposalargs = NULL,
                        dissolve=NULL, gamma=-4.59512,
                        dissolve.order="DissThenForm",
                        algorithm.control=list(),
                        drop=FALSE,
                        verbose=FALSE) {
  formula <- object

  ## Defaults :
  con <- list(boundDeg=NULL, drop=drop,
              dyninterval=1000,
              maxchanges=1000000,
              final=FALSE,
              summarizestats=FALSE
             )

  con[(namc <- names(algorithm.control))] <- algorithm.control
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formula)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula must be given")
  }
  # Resolve conditioning
  # This is laborious to cover the various partial specifications
  #
  BD <- ergm.boundDeg(con$boundDeg, nnodes=network.size(nw))

  model <- ergm.getmodel(formula, nw, drop=con$drop)
  model.dissolve <- ergm.getmodel.dissolve(dissolve, nw, dissolve.order)
#
  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta0)) {
    theta0 <- rep(0,length(model$coef.names))
    warning("No parameter values given, using Bernouli network\n\t")
  }
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
    
  MHproposal <- getMHproposal(proposaltype, proposalargs, nw, model)
  MCMCparams <- c(con,list(samplesize=1, interval=interval,
                         stats=matrix(0,ncol=length(model$coef.names),nrow=1),
                         burnin=nsim,
                         parallel=0,
                         meanstats=theta0-theta0,
                         gamma=gamma))
  
  z <- ergm.getMCMCDynsample(nw, model, model.dissolve, MHproposal,
                             theta0, MCMCparams, verbose, BD)

  if(con$final){
   nw <- network.update(nw,z$newedgelist)
   return(nw)
  }else{
    out.list <- list(formula = formula, networks = nw,
                     changed=z$changed, 
                     dissolved=z$dissolved, 
                     maxchanges=z$maxchanges,
                     stats = NULL, coef=theta0)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
