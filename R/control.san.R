  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
