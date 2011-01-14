#==========================================================
# This file contains the following 2 functions for
# controlling simulation routines
#        <control.simulate> = <control.simulate.formula>
#        <control.simulate.ergm>
#==========================================================


########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                   if being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#   maxedges      : the maximum number of edges expected in network;
#                   default=20000
#   packagenames  : the names of packages to load on the created
#                   cluster when using parallel threads; currently,
#                   the only recognized package name is "ergm";
#                   default="ergm"
#   network.output: the R class with which to output networks; the
#                   options are "network" and "edgelist.compressed"
#                   (which saves space but only supports networks
#                    without vertex attributes); default="network"

#
# --IGNORED--
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <ergm.san> which is the only code using this
#                   control list.
#   drop          : whether degenerate terms should be dropped from the
#                   fit (T or F); default=TRUE
#   summarizestats: whether to print out a summary of the sufficient
#                   statistics of the generated network (T or F);
#                   default=FALSE.
#   maxchanges    : ??; default=1000000
#   parallel      : number of threads in which to run the sampling;
#                   not implemented in the 2.3 version of <ergm.getMCMCsample>
#
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.simulate<-control.simulate.formula<-function(prop.weights="default",
                                                     prop.args=NULL,
                                                     drop=FALSE,
                                                     summarizestats=FALSE,
                                                     maxedges=20000,
                                                     maxchanges=1000000,
                                                     packagenames="ergm",
                                                     network.output="network",
                                                     parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

control.simulate.ergm<-function(prop.weights="default",
                                prop.args=NULL,
                                drop=FALSE,
                                summarizestats=FALSE,
                                maxchanges=1000000,
                                maxedges=20000,
                                packagenames="ergm",
                                network.output="network",
                                parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
