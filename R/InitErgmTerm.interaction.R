#  File R/InitErgmTerm.interaction.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
check_interact_term <- function(m, dependent_action){
  msg <- paste0("Change statistic interactions are poorly defined for dyad-dependent terms. Use ", sQuote("interact.dependent"), " term option to set the behavior.")
  if(!is.dyad.independent(m))
    switch(dependent_action,
           error = ergm_Init_stop(msg, call.=FALSE),
           warning = ergm_Init_warning(msg, immediate.=TRUE, call.=FALSE), # Warn immediately, so the user gets the warning before the MCMC starts.
           message = message(msg))

  if(is.curved(m)) ergm_Init_stop("Interactions are undefined for curved terms at this time.")
}

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.:` <- function(nw, arglist, ..., env, interact.dependent = c("error", "message", "warning", "silent")){
  arglist <- substitute(arglist)
  e1 <- arglist[[2]]
  e2 <- arglist[[3]]

  e1 <- list_summands.call(e1)
  e2 <- list_summands.call(e2)

  n1 <- length(e1)
  n2 <- length(e2)

  f <- append_rhs.formula(NULL, c(e1, e2), env = env)
  m <- ergm_model(f, nw, ..., offset.decorate=FALSE)

  check_interact_term(m, match.arg(interact.dependent))

  cn1 <- unlist(lapply(m$terms[seq_len(n1)], "[[", "coef.names"))
  cn2 <- unlist(lapply(m$terms[n1+seq_len(n2)], "[[", "coef.names"))

  inputs <- c(length(cn1), length(cn2))
  
  cn <- outer(cn1,cn2,paste,sep=":")

  wm <- wrap.ergm_model(m, nw, NULL)
  if(any(wm$offsettheta) || any(wm$offsetmap)) ergm_Init_warning("The interaction operator does not propagate offset() decorators.")

  c(list(name="interact", coef.names = cn, inputs=inputs, submodel=m, dependence=wm$dependence),
    ergm_propagate_ext.encode(m))

}

## This will always be passed with two arguments in arglist, which
## will cause an error if we actually try to evaluate them. So,
## there's no check.ErgmTerm() but rather an immediate substitute() to
## grab the actual names or calls being passed.
`InitErgmTerm.*` <- function(nw, arglist, ..., env, interact.dependent = c("error", "message", "warning", "silent")){
  arglist <- substitute(arglist)
  e1 <- arglist[[2]]
  e2 <- arglist[[3]]

  e1 <- list_summands.call(e1)
  e2 <- list_summands.call(e2)

  n1 <- length(e1)
  n2 <- length(e2)

  f <- append_rhs.formula(NULL, c(e1, e2), env = env)
  m <- ergm_model(f, nw, ..., offset.decorate=FALSE)

  check_interact_term(m, match.arg(interact.dependent))

  cn1 <- unlist(lapply(m$terms[seq_len(n1)], "[[", "coef.names"))
  cn2 <- unlist(lapply(m$terms[n1+seq_len(n2)], "[[", "coef.names"))

  inputs <- c(length(cn1), length(cn2))

  cn <- c(cn1,cn2,outer(cn1,cn2,paste,sep=":"))
  
  wm <- wrap.ergm_model(m, nw, NULL)
  if(any(wm$offsettheta) || any(wm$offsetmap)) ergm_Init_warning("The interaction operator does not propagate offset() decorators.")

  c(list(name="main_interact", coef.names = cn, inputs=inputs, submodel=m, dependence=wm$dependence),
    ergm_propagate_ext.encode(m))
}
