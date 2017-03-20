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

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.:` <- function(nw, arglist, response=NULL, ...){
  arglist <- substitute(arglist)
  e1 <- arglist[[2]]
  e2 <- arglist[[3]]

  e1 <- term.list.formula(e1)
  e2 <- term.list.formula(e2)

  n1 <- length(e1)
  n2 <- length(e2)
  
  f <- ~nw
  f <- append.rhs.formula(f, c(e1,e2))
  
  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  if(!is.dyad.independent(m)) stop("Interactions are not meaningful for dyad-dependent terms.")

  cn1 <- unlist(lapply(m$terms[seq_len(n1)], "[[", "coef.names"))
  cn2 <- unlist(lapply(m$terms[n1+seq_len(n2)], "[[", "coef.names"))

  inputs <- c(length(cn1), length(cn2), pack.Clistasnum(Clist))
  
  cn <- outer(cn1,cn2,paste,sep=":")
  
  list(name="interact", coef.names = cn, inputs=inputs, dependence=FALSE)
}

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.*` <- function(nw, arglist, response=NULL, ...){
  arglist <- substitute(arglist)

  e1 <- term.list.formula(e1)
  e2 <- term.list.formula(e2)

  n1 <- length(e1)
  n2 <- length(e2)
  
  f <- ~nw
  f <- append.rhs.formula(f, c(e1,e2))
  
  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  if(!is.dyad.independent(m)) stop("Interactions are not meaningful for dyad-dependent terms.")
  
  cn1 <- unlist(lapply(m$terms[seq_len(n1)], "[[", "coef.names"))
  cn2 <- unlist(lapply(m$terms[n1+seq_len(n2)], "[[", "coef.names"))

  inputs <- c(length(cn1), length(cn2), pack.Clistasnum(Clist))

  cn <- c(cn1,cn2,outer(cn1,cn2,paste,sep=":"))
  
  list(name="main_interact", coef.names = cn, inputs=inputs, dependence=FALSE)
}
