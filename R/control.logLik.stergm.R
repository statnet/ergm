#  File ergm/R/control.logLik.stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
control.logLik.stergm<-function(control = control.logLik.ergm(),
                                control.form = control,
                                control.diss = control,
                                
                                seed=NULL){
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  control
}
