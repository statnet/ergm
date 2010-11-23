#=================================================================
# This file contains the 2 following functions for controlling
# goodness-of-fit computations
#            <control.gof.ergm>
#            <control.gof.formula>
#=================================================================



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
