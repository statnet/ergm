control.simulatedyn<-control.simulatedyn.ergm<-function(prop.weights.form="default",
                                                        prop.args.form=NULL,
                                                        prop.weights.diss="default",
                                                        prop.args.diss=NULL,
                                                        drop=FALSE,
                                                        summarizestats=FALSE,final=FALSE,
                                                        maxchanges=1000000){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[[arg]]<-get(arg)
    control
  }
