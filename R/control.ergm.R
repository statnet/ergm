control.ergm<-function(prop.weights="default",prop.args=NULL,
                       nr.maxit=100, calc.mcmc.se=TRUE, hessian=FALSE,
                       compress=TRUE,
                       maxNumDyadTypes=1e+6, 
                       maxedges=20000,
                       maxchanges=1000000,
                       maxMPLEsamplesize=100000,
                       MPLEtype=c("glm", "penalized"),
                       trace=0,
                       steplength=0.5,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       check.degeneracy=TRUE,
                       mcmc.precision=0.05,
                       metric=c("Likelihood","raw"),
                       method=c("BFGS","Nelder-Mead"),
                       trustregion=20,
                       initial.loglik=NULL,
                       initial.network=NULL,
                       style=c("Newton-Raphson","Robbins-Monro","Stochastic-Approximation"),
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       RobMon.phase1n_base=7,
                       RobMon.phase2n_base=7,
                       RobMon.phase2sub=4,
                       RobMon.init_gain=0.4,
                       RobMon.phase3n=500,
                       dyninterval=1000,
                       packagenames="ergm",
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
  if(missing(trustregion) & control$style=="Stochastic-Approximation"){
   control$trustregion <- 0.5
  }

  control
}
