
InitWtErgmTerm.import <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw,...)

  if(!is.dyad.independent(m)) stop("Only dyad-independent binary terms can be imported at this time.")
  
  Clist <- ergm.Cprepare(nw, m)

  fnames <- pack.strtoint(Clist$fnamestring)
  snames <- pack.strtoint(Clist$snamestring)

  inputs <- c(Clist$nterms, fnames, snames, Clist$inputs)

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

  list(name="import_binary_term_sum", coef.names = paste0('sum.',m$coef.names,''), inputs=inputs, dependence=FALSE, emptynwstats = gs)
}
