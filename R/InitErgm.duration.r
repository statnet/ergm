#  See InitErgm.r for a general explanation 
#  of InitErgm functions

InitErgm.D.off <- function(nw, m, arglist, ...) {
  a=ergm.checkargs("D.off", arglist,
                   varnames = NULL,
                   vartypes = NULL,
                   defaultvalues = list(),
                   required = NULL)
  m$DURATIONflag <- TRUE
  m$coef.names<-c(m$coef.names, "D.off")
  termnumber <- 1 + length(m$terms)
  m$terms[[termnumber]] <- list(name = "D_off", soname="statnet",
                                inputs = c(0, 1, 0))
  m
}
  
InitErgm.D.dyad <- function(nw, m, arglist, ...) {
  a=ergm.checkargs("D.dyad", arglist,
                   varnames = NULL,
                   vartypes = NULL,
                   defaultvalues = list(),
                   required = NULL)
  m$DURATIONflag <- TRUE
  m$coef.names <- c(m$coef.names, "D.dyad")
  termnumber <- 1 + length(m$terms)
  m$terms[[termnumber]] <- list(name = "D_dyad", soname="statnet",
                                inputs = c(0, 1, 0))
  m
}

InitErgm.D.edge <- function(nw, m, arglist, ...) {
  a=ergm.checkargs("D.edge", arglist,
                   varnames = NULL,
                   vartypes = NULL,
                   defaultvalues = list(),
                   required = NULL)
  m$DURATIONflag <- TRUE
  m$coef.names <- c(m$coef.names, "D.edge")
  termnumber <- 1 + length(m$terms)
  m$terms[[termnumber]] <- list(name = "D_edge", soname="statnet",
                                inputs = c(0, 1, 0))
  m
}


