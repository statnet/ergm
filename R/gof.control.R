gof.ergm.control<-function(prop.weights=NULL,prop.args=NULL,
                           drop=TRUE,
                           summarizestats=FALSE,
                           maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

gof.formula.control<-function(prop.weights="default",prop.args=NULL,
                              drop=TRUE,
                              summarizestats=FALSE,
                              maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
