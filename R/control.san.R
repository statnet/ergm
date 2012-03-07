#  File ergm/R/control.san.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# The <control.san> function creates a list of paramaters
# for customizing the <ergm.san> routines
########################################################################
control.san<-function(coef=NULL,

                      SAN.tau=1,
                      SAN.invcov=NULL,
                      SAN.burnin=10000,
                      SAN.interval=10000,
                      SAN.init.maxedges=20000,
                      
                      SAN.prop.weights="default",
                      SAN.prop.args=list(),
                      SAN.packagenames="ergm",

                      network.output="network",

                      seed=NULL,
                      parallel=0,
                      parallel.type=NULL,
                      parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))
  control
}
