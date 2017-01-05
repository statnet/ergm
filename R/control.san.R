#  File R/control.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
########################################################################
# The <control.san> function creates a list of paramaters
# for customizing the <ergm.san> routines
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                  of being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#   maxedges      : the maximum number of edges expected in network
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
#   packagenames  : the names of packages in which changestatistics
#                   are found; currently ignored as 'ergm' is presumed;
#                   default="ergm"
#   maxchanges    : ??; default=1000000
#   parallel      : number of threads in which to run the sampling; 
#                   currently unimplemented in <ergm.san> or the
#                   subsequent call to <ergm.mple>
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.san<-function(coef=NULL,

                      SAN.tau=1,
                      SAN.invcov=NULL,
                      SAN.burnin=100000,
                      SAN.interval=10000,
                      SAN.init.maxedges=20000,
                      
                      SAN.prop.weights="default",
                      SAN.prop.args=list(),
                      SAN.packagenames=c(),
                      
                      MPLE.max.dyad.types=1e6,
                      MPLE.samplesize=50000,

                      network.output="network",

                      seed=NULL,
                      parallel=0,
                      parallel.type=NULL,
                      parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))

  set.control.class()
}
