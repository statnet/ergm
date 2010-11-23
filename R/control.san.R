control.san<-function(prop.weights="default",
                      prop.args=NULL,
                      drop=FALSE,
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
