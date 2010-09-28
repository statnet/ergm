control.san<-function(prop.weights="default", prop.args=NULL,
                      drop=FALSE,
                      network.output="network",
                      maxchanges=1000000){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}


