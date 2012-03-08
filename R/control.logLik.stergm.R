
control.logLik.stergm<-function(control.form = control.logLik.ergm(),
                                control.diss = control.logLik.ergm(),
                                control = NULL,
                                
                                seed=NULL){
  if(!is.null(control))
    control.form <- control.diss <- control
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class()
}
