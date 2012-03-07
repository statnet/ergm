#  File ergm/R/control.simulate.stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
########################################################################
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
