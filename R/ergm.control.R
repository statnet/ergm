ergm.control<-function(prop.weights="default",prop.args=NULL,
                       prop.weights.diss="default",prop.args.diss=NULL,
                       nr.maxit=100, calc.mcmc.se=TRUE, hessian=FALSE,
                       compress=FALSE,
                       maxNumDyadTypes=10000, 
                       maxedges=20000,
                       maxchanges=1000000,
                       MPLEsamplesize=50000, 
                       MPLEtype=c("glm", "penalized"),
                       trace=0,
                       steplength=0.5,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       mcmc.precision=0.05,
                       metric=c("Likelihood","raw"),
                       method=c("BFGS","Nelder-Mead"),
                       trustregion=20,
                       style=c("Newton-Raphson","Robbins-Monro","Stochastic-Approximation"),
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       RobMon.phase1n_base=7,
                       RobMon.phase2n_base=7,
                       RobMon.phase2sub=4,
                       RobMon.init_gain=0.4,
                       RobMon.phase3n=500,
                       dyninterval=1000,
                       parallel=0,
                       returnMCMCstats=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  
  control$MPLEtype<-match.arg(MPLEtype)
  control$metric<-match.arg(metric)
  control$method<-match.arg(method)
  control$style<-match.arg(style)
  control$nsubphases<-match.arg(nsubphases)

  control
}
