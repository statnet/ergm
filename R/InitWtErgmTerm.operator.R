#  File R/InitWtErgmTerm.operator.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
InitWtErgmTerm.Passthrough <- function(nw, arglist, ...){
  out <- InitErgmTerm.Passthrough(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

InitWtErgmTerm.B <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "character,formula"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form <- if(is.character(a$form)) match.arg(a$form,c("sum","nonzero"))
          else a$form

  nwb <- nw
  nwb %ergmlhs% "response" <- NULL
  m <- ergm_model(a$formula, nwb, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  if(!is.dyad.independent(m) && form=="sum") stop("Only dyad-independent binary terms can be imported with form 'sum'.")
  
  if(is(form, "formula")){
    form.name <- despace(deparse(ult(form)))
    name <- "import_binary_term_form"
    auxiliaries <- trim_env(~.binary.formula.net(form),"form")
  }else{
    form.name <- form
    name <- paste("import_binary_term",form,sep="_")
    auxiliaries <- if(form=="nonzero") trim_env(~.binary.nonzero.net)
  }

  mw <- wrap.ergm_model(m, nwb, ergm_mk_std_op_namewrap('B', form.name))
  if(form=="sum") mw$emptynwstats <- NULL

  c(list(name=name,
         submodel = m,
         auxiliaries=auxiliaries),
    mw)
}

InitWtErgmTerm..binary.nonzero.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_binary_nonzero_net", depenence=FALSE)
}

InitWtErgmTerm..binary.formula.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  if(!is.dyad.independent(m) || nparam(m)!=1) stop("The binary test formula must be dyad-independent and have exactly one statistc.")

  nw[,] <- FALSE
  gs <- summary(m, nw)
  if(gs!=0) stop("At this time, the binary test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  c(list(name="_binary_formula_net", submodel=m, depenence=FALSE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, NULL))
}

# Arguments and outputs are identical to the binary version, except for the C routine names.
InitWtErgmTerm.Sum <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.Sum
  term <- f(...)
  term$name <- "wtSum"
  term
}

InitWtErgmTerm.Label <- function(nw, arglist, ...){
  out <- InitErgmTerm.Label(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

InitWtErgmTerm.Parametrise <- InitWtErgmTerm.Parametrize <- InitWtErgmTerm.Curve <- function(nw, arglist, ...){
  out <- InitErgmTerm.Curve(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

InitWtErgmTerm..submodel_and_summary <- function(nw, arglist, ...){
  out <- InitErgmTerm..submodel_and_summary(nw, arglist, ...)
  out$name <- "_wtsubmodel_and_summary_term"
  out
}


InitWtErgmTerm.Exp <- function(nw, arglist, ...){
  out <- InitErgmTerm.Exp(nw, arglist, ...)
  out$name <- "wtExp"
  out
}

InitWtErgmTerm.Log <- function(nw, arglist, ...){
  out <- InitErgmTerm.Log(nw, arglist, ...)
  out$name <- "wtLog"
  out
}

InitWtErgmTerm.Prod <- InitErgmTerm.Prod
