InitErgmTerm..edges_times<-function(nw, arglist, ..., times) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="_edges_times", coef.names=paste("edges_times",times), dependence=FALSE,
       minval = min(times*network.dyadcount(nw,FALSE),0), maxval = max(times*network.dyadcount(nw,FALSE),0), conflicts.constraints="edges", inputs=c(times))
}
