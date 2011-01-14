#=================================================================
# This file contains the 2 following functions for controlling
# goodness-of-fit computations
#            <control.gof.ergm>
#            <control.gof.formula>
#=================================================================



########################################################################
# Both of the <control.gof.X> functions return a control list for
# customing the fitting procedure used by the gof code
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                   if being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#   drop          : whether degenerate terms should be dropped from the
#                   fit (T or F); default=TRUE
#   summarizestats: whether to print out a summary of the sufficient
#                   statistics of the generated network (T or F);
#                   default=FALSE
#
# --IGNORED--
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <ergm.gof> which is the only code using this
#                   control list.
#   maxchanges    : ??; default=1000000
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.gof.ergm<-function(prop.weights=NULL,prop.args=NULL,
                           drop=TRUE,
                           summarizestats=FALSE,
                           maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}


control.gof.formula<-function(prop.weights="default",prop.args=NULL,
                              drop=TRUE,
                              summarizestats=FALSE,
                              maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
