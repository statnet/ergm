.encode.str <- function(s) c(nchar(s), strtoi(charToRaw(s), 16L))

InitErgmTerm.meta <- function(nw, arglist, response=NULL, ...){
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

  fnames <- .encode.str(Clist$fnamestring)
  snames <- .encode.str(Clist$snamestring)

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

  list(name="meta_term", coef.names = paste0('meta(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs)
}

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.:` <- function(nw, arglist, response=NULL, ...){
  arglist <- substitute(arglist)
  n1 <- arglist[[2]]
  n2 <- arglist[[3]]

  f <- ~nw
  f[[3]] <- call("+",n1,n2)
  
  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  fnames <- .encode.str(Clist$fnamestring)
  snames <- .encode.str(Clist$snamestring)

  inputs <- c(fnames, snames, Clist$inputs)

  cn <- outer(m$terms[[1]]$coef.names,m$terms[[2]]$coef.names,paste,sep=":") 
  
  list(name="interact", coef.names = cn, inputs=inputs, dependence=!is.dyad.independent(m))
}

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.*` <- function(nw, arglist, response=NULL, ...){
  arglist <- substitute(arglist)
  n1 <- arglist[[2]]
  n2 <- arglist[[3]]

  f <- ~nw
  f[[3]] <- call("+",n1,n2)
  
  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  fnames <- .encode.str(Clist$fnamestring)
  snames <- .encode.str(Clist$snamestring)

  inputs <- c(fnames, snames, Clist$inputs)

  cn <- c(m$terms[[1]]$coef.names,m$terms[[2]]$coef.names,outer(m$terms[[1]]$coef.names,m$terms[[2]]$coef.names,paste,sep=":"))
  
  list(name="main_interact", coef.names = cn, inputs=inputs, dependence=!is.dyad.independent(m))
}
