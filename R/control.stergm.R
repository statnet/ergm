###########################################################################
# The <control.stergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
#
# --PARAMETERS--
#   prop.weights.form: the method to allocate probabilities of being proposed
#                      to dyads in the formation stage, as "TNT", "random",
#                      "nonobserved", or "default"
#                      default="default", which is based upon the ergm constraints
#   prop.args.form   : an alternative, direct way of specifying additional
#                      arguments to the formation proposals
#   prop.weights.diss: as 'prop.weights.form', but for the dissolution model
#   prop.args.diss   : as 'prop.args.form', but for the dissoultion model
#   compress         : whether the stats matrix should be compressed to the set
#                      of unique statistics with a column of probability weights
#                      post-pended; default=FALSE
#   SAN.burnin       : the burnin value used to create the SAN-ed network and
#                      formula; default=NULL
#   SAN.interval     :
#   maxNumDyadTypes  : the maximum number of unique psuedolikelihood change stats
#                      to be allowed if 'compress'=TRUE; ignored if 'compress'!=TRUE;
#                      default=1e+6
#   maxedges         : the maximum number of edges to allocate space for; default=20000
#   maxchanges       : the maximum number of changes in dynamic network simulation for
#                      which to allocate space; default=1000000
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the psuedo-
#                      likelihood computation; default=100000
#   maxNumDyadTypes  : the maximum number of unique psuedolikelihood change stats
#                      to be allowed if 'compress'=TRUE; ignored if 'compress'!=TRUE;
#                      default=1e+6
#   maxedges         : the maximum number of edges to allocate space for; default=20000
#   maxchanges       : the maximum number of changes in dynamic network simulation for
#                      which to allocate space; default=1000000
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the psuedo-
#                      likelihood computation; default=100000
#   MPLEtype         : the method for MPL estimation as "penalized", "glm" or
#                      "logitreg"; default="glm"
#   trace            : the number of levels of tracing information to produce during
#                      optimization; see <?optim> for details; default=0
#   sequential       : whether the next iteration of the fit should use the last network
#                      sampled as the starting point; the alternative is to always begin
#                      from the orginial network; default=TRUE
#   drop             : whether degenerate terms should be dropped from the fit (T or F);
#                      default=TRUE
#   style            : the style of ML estimation to use, as one of "Newton-Raphson",
#                      "Robbins-Monro", "Stochastic-Approximation", or "Stepping";
#                      default="Robbins-Monro"
#   RM.init_gain     : this is only used to adjust 'aDdiaginv'in phase1,
#                      in particular:
#                             aDdiaginv = gain/sqrt(aDdiaginv)
#                      default=0.5
#   RM.phase1n_base  : this helps define the 'phase1n' param, which in turn
#                      multiplies 'RM.interval' to control the number of
#                      phase1 iterations; this is the base portion of 'phase1n',
#                      which is added to 3*(the number of formation coefficients)
#                      to form 'phase1n'; default=7
#   RM.phase2sub     : phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
#                      the outer loop counter; default=7
#   RM.phase2n_base  : this helps define the 'phase2n' param, which in turn
#                      limits the phase2 middle loop counter; this is the
#                      base portion of 'phase2n', which is added to 7+(the number
#                      of formation coefficients) to form 'phase2n'; default=100
#   RM.burnin        : the number of MCMC steps to disregard for the burn-in
#                      period; default=1000
#   RM.interval      : like the SPSA.interval, this seems a little more like
#                      a sample size, than an interval, it helps control the
#                      number of MCMCsteps used in phase1 and phase2; in
#                      phase2, this limits the innermost loop counter; default=100
#   SPSA.iterations  : the number of iterations to use in the SPSA sampling;
#                      default=1000
#   SPSA.a           : see the next 2 params; default=1
#   SPSA.alpha       : see the next 2 params; default=.602
#   SPSA.A           : this and the 2 params above help to define the
#                      'gain' paramater as
#                          SPSA.a/(SPSA.A +i +1)^(SPSA.alpha)
#                      where i is indexed from 0 to SPSA.iterations;
#                      default=100
#   SPSA.c           : see the next param; default=1
#   SPSA.gamma       : this and the param above help to define the 'diff'
#                      parameter as
#                           SPSA.c/(i+1)^(SPSA.gamma)
#                      where i is indexed from 0 to SPSA.iterations;
#                      default=.101
#   SPSA.burnin      : the number of MCMC steps to disregard for the burnin
#                      period; default=1000
#   SPSA.interval    : this is eventually received as 'S' and looks like a
#                      a sample size, rather than an interval, since 'S' controls
#                      the number of MCMC steps that contribute to the stats vector
#   NM.abstol        : ??    ; default==0,
#   NM.reltol        : ??    ; default==sqrt(.Machine$double.eps),
#   NM.alpha         : ??    ; default==1,
#   NM.beta          : ??    ; default==.5,
#   NM.gamma         : ??    ; default==2,
#   NM.maxit         : ??    ; default==500,
#   NM.interval      : ??    ; default==1000,
#   NM.burnin        : ??    ; default==1000,
#   packagenames     : the packages in which change statistics are found; default="ergm"
#   parallel         : the number of threads in which to run sampling; default=0
#
# --RETURNED--
#   a list of the above parameters
#
######################################################################################################

control.stergm<-function(prop.weights.form="default",prop.args.form=NULL,
                         prop.weights.diss="default",prop.args.diss=NULL,
                         compress=FALSE,
                         SAN.burnin=10000,
                         SAN.interval=1000,
                         maxNumDyadTypes=1e+6, 
                         maxedges=20000,
                         maxchanges=1000000,
                         maxMPLEsamplesize=100000,
                         MPLEtype=c("glm", "penalized"),
                         trace=0,
                         sequential=TRUE,
                         drop=TRUE,
                         style=c("Robbins-Monro","SPSA", "SPSA2", "Nelder-Mead"),
                         RM.phase1n_base=7,
                         RM.phase2n_base=100,
                         RM.phase2sub=4,
                         RM.init_gain=0.5,
                         RM.interval=100,
                         RM.burnin=1000,
                         SPSA.a=1,
                         SPSA.alpha=0.602,
                         SPSA.A=100,
                         SPSA.c=1,
                         SPSA.gamma=0.101,
                         SPSA.iterations=1000,
                         SPSA.interval=1000,
                         SPSA.burnin=1000,
                         NM.abstol=0,
                         NM.reltol=sqrt(.Machine$double.eps),
                         NM.alpha=1,
                         NM.beta=.5,
                         NM.gamma=2,
                         NM.maxit=500,
                         NM.interval=1000,
                         NM.burnin=1000,
                         packagenames="ergm",
                         parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  
  control$MPLEtype<-match.arg(MPLEtype)
  control$style<-match.arg(style)
  control
}
