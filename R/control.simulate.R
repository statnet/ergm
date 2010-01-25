control.simulate<-control.simulate.formula<-function(prop.weights="default",
                                                     prop.args=NULL,
                                                     drop=FALSE,
                                                     summarizestats=FALSE,
                                                     maxchanges=1000000,
                                                     packagenames="ergm",
                                                     parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

control.simulate.ergm<-function(prop.weights=NULL,
                                prop.args=NULL,
                                drop=FALSE,
                                summarizestats=FALSE,
                                maxchanges=1000000,
                                packagenames="ergm",
                                parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
