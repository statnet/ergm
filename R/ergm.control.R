ergm.control<-function(prop.weights="default",prop.args=NULL,
                       prop.weights.diss="default",prop.args.diss=NULL,
                       nr.maxit=100, calc.mcmc.se=TRUE, hessian=FALSE,
                       compress=FALSE,
                       maxNumDyadTypes=10000, 
                       maxedges=20000,
                       maxchanges=1000000,
                       MPLEsamplesize=50000, 
                       trace=0,
                       steplength=0.5,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       mcmc.precision=0.05,
                       metric="Likelihood",
                       method="BFGS",
                       trustregion=20,
                       style="Newton-Raphson",
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       dyninterval=1000,
                       parallel=0,
                       returnMCMCstats=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
      control[[arg]]<-get(arg)
  control
}
