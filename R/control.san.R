control.san<-function(prop.weights="default",prop.args=NULL,
                      drop=FALSE,
                      summarizestats=FALSE,
                      maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

control.san.ergm<-function(prop.weights=NULL,prop.args=NULL,
                           drop=FALSE,
                           summarizestats=FALSE,
                           maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
