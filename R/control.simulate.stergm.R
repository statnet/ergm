########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
#
# --PARAMETERS--
#   prop.weights.form/
#   ptop.weighs.diss: specifies the method used to allocate probabilities
#                     of being proposed to dyads for the formation/dis-
#                     solution processes; options are "TNT",
#                     "random", and "default"; default="default" which
#                      picks a reasonable default considering any constraints
#   final           : whether only the final network of the simulation
#                     process should be returned (T or F); default= FALSE
#                     in which case, models, coefficients, stats matrices,
#                     and the toggle matrix are returned
#   maxchanges      : the maximum number of changes for which to allocate
#                     space; default=1000000
#
# --IGNORED--
#   prop.args.form/
#   prop.args.diss: an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <simulate.stergm> which is the only code using this
#                   control list.
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.simulate.stergm<-control.simulate.network<-function(MCMC.burnin=1000,
                                                            MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                                                            MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,                                  
                                  MCMC.init.maxedges=20000,
                                  MCMC.packagenames="ergm",

                                  MCMC.init.maxchanges=1000000){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))
    control
  }
