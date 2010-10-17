control.ergm<-function(prop.weights="default",prop.args=NULL,
                       prop.weights.diss="default",prop.args.diss=NULL,
                       nr.maxit=100,
                       nr.reltol=sqrt(.Machine$double.eps),
                       calc.mcmc.se=TRUE, hessian=TRUE,
                       compress=FALSE,
                       SAN.burnin=NULL,
                       maxNumDyadTypes=1e+6, 
                       maxedges=20000,
                       maxchanges=1000000,
                       maxMPLEsamplesize=100000,
                       miss.MCMCsamplesize=NULL,
                       miss.interval=NULL,
                       miss.burnin=NULL,
                       MPLEtype=c("glm", "penalized"),
                       trace=0,
                       steplength=0.5,
                       adaptive.trustregion=3,
                       adaptive.epsilon=0.01,
                       sequential=TRUE,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       check.degeneracy=FALSE,
                       mcmc.precision=0.05,
                       metric=c("Median.Likelihood", "lognormal", 
                                "EF.Likelihood", "naive"),
                       method=c("BFGS","Nelder-Mead"),
                       trustregion=20,
                       initial.loglik=NULL,
                       initial.network=NULL,
                       style=c("Newton-Raphson","Robbins-Monro",
                               "Stochastic-Approximation","Stepping","PILA"),
                       style.dyn=c("Robbins-Monro","SPSA", "SPSA2"),
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       RobMon.phase1n_base=7,
                       RobMon.phase2n_base=100,
                       RobMon.phase2sub=7,
                       RobMon.init_gain=0.5,
                       RobMon.phase3n=500,
                       PILA.gamma=.99,
                       PILA.steplength=.1,
                       SPSA.a=1,
                       SPSA.alpha=0.602,
                       SPSA.A=100,
                       SPSA.c=1,
                       SPSA.gamma=0.101,
                       SPSA.iterations=1000,
                       stepMCMCsize=100,
                       gridsize=100,
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
  control$style.dyn<-match.arg(style.dyn)
  control$nsubphases<-match.arg(nsubphases)
  if(missing(trustregion) & control$style=="Stochastic-Approximation"){
   control$trustregion <- 0.5
  }

  control
}
