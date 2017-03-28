pack.strasnum <- function(s) c(nchar(s), strtoi(charToRaw(s), 16L))
pack.Clistasnum <- function(Clist){
  fnames <- pack.strasnum(Clist$fnamestring)
  snames <- pack.strasnum(Clist$snamestring)
  c(Clist$nterms, fnames, snames, Clist$inputs)
}


InitErgmTerm.passthrough <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clistasnum(Clist)

  gs <- rep(0, Clist$nstats)
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$terms[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }

  list(name="passthrough_term", coef.names = paste0('meta(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs)
}
