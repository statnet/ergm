
InitWtErgmTerm.b <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "character,formula"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form <- if(is.character(a$form)) match.arg(a$form,c("sum","nonzero"))
          else a$form

  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw,...)

  if(!is.dyad.independent(m) && form=="sum") stop("Only dyad-independent binary terms can be imported with form 'sum'.")
  
  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)
  
  gs <- if(form=="nonzero" || is(form, "formula")) ergm.emptynwstats.model(m)
        else rep(0, eta.length.model(m))

  if(is(form, "formula")){
    form.name <- deparse(form[[length(form)]])
    name <- "import_binary_term_form"
    auxiliaries <- ~.binary.formula.net(form)
  }else{
    form.name <- form
    name <- paste("import_binary_term",form,sep="_")
    auxiliaries <- if(form=="nonzero") ~.binary.nonzero.net
  }
  
  list(name=name,
       coef.names = paste0(form.name,'(',m$coef.names,')'),
       inputs=inputs,
       dependence=!is.dyad.independent(m),
       emptynwstats = gs,
       auxiliaries=auxiliaries)
}

InitWtErgmTerm..binary.nonzero.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_binary_nonzero_net", depenence=FALSE)
}

InitWtErgmTerm..binary.formula.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  # Form is a model.
  f<-a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm.getmodel(f, nw, response=response,...)

  if(!is.dyad.independent(m) || coef.length.model(m)!=1) stop("The binary test formula must be dyad-independent and have exactly one statistc.")

  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)
  if(gs!=0) stop("At this time, the binary test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  list(name="_binary_formula_net", inputs=c(inputs), depenence=FALSE)
}
