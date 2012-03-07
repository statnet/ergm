
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
