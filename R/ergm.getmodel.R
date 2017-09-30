#  File R/ergm.getmodel.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#===================================================================================
# This file contains the following 2 functions for creating the 'ergm.model' object
#             <ergm.getmodel>
#             <updatemodel.ErgmTerm>
#===================================================================================





###################################################################################
# The <ergm.getmodel> function parses the given formula, and initiliazes each ergm
# term via the <InitErgmTerm> functions to create a 'model.ergm' object for the
# given network
#
# --PARAMETERS--
#   formula:  a formula of the form 'network ~ model.term(s)'
#   nw     :  the network of interest
#   silent :  whether to print the warning messages from the initialization of each
#             model term (T or F); default=FALSE
#   ...    :  additional parameters for model formulation;
#             recognized parameters include
#               initialfit: whether curved exponential terms have been initially fit
#                           by MPLE (T or F)
#
# --RETURNED--
#   a 'model.ergm' object as a list containing:
#     formula       :  the formula inputted to <ergm.getmodel>
#     coef.names    :  a vector of coefficient names
#     offset        :  a logical vector of whether each term was "offset", i.e. fixed
#     terms         :  a list of terms and 'term components' initialized by the 
#                      appropriate <InitErgmTerm.X> function.  See the <InitErgm> 
#                      function header for details about the 'terms' list
#     network.stats0:  NULL always??
#     etamap        :  the theta -> eta mapping as a list returned from <ergm.etamap> 
#     class         :  the character string "model.ergm" 
#
#####################################################################################

ergm.getmodel <- function (formula, nw, response=NULL, silent=FALSE, role="static",...,extra.aux=list()) {
  if ((dc<-data.class(formula)) != "formula")
    stop (paste("Invalid formula of class ",dc), call.=FALSE)

  if (formula[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.", call.=FALSE)  
  
  if (length(formula) < 3) 
    stop(paste("No model specified for network ", formula[[2]]), call.=FALSE)

  v<-term.list.formula(formula[[3]])
  
  formula.env<-environment(formula)
  
  model <- structure(list(formula=formula, coef.names = NULL,
                      offset = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "model.ergm")
  
  for (i in 1:length(v)) {
    term <- v[[i]]
    
    if (is.call(term) && term[[1]] == "offset"){ # Offset term
      term <- term[[2]]
      model$offset <- c(model$offset,TRUE)
    }else{
      model$offset <- c(model$offset,FALSE)
    }
    ## term is now a call or a name that is not "offset".


    outlist <- call.ErgmTerm(term, formula.env, nw, response=response, role=role, ...)
    

    # If initialization fails without error (e.g., all statistics have been dropped), continue.
    if(is.null(outlist)){
      if(!silent) message("Note: Term ", deparse(v[[i]])," skipped because it contributes no statistics.")
      model$term.skipped <- c(model$term.skipped, TRUE)
      next
    }else model$term.skipped <- c(model$term.skipped, FALSE)
    # Now it is necessary to add the output to the model object
    model <- updatemodel.ErgmTerm(model, outlist)
  }

  model <- ergm.auxstorage(model, nw, response=response, ..., extra.aux=extra.aux)
  
  model$etamap <- ergm.etamap(model)

  # I.e., construct a vector of package names associated with the model terms.
  # Note that soname is not the same, since it's not guaranteed to be a loadable package.
  ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
  ergm.MCMC.packagenames(unlist(sapply(model$model.aux$terms, "[[", "pkgname")))
  
  class(model) <- "ergm.model"
  model
}


call.ErgmTerm <- function(term, env, nw, response=NULL, role="static", ...){
  termroot<-if(is.null(response)) "InitErgm" else "InitWtErgm"
  
  if(is.call(term)) { # This term has some arguments; save them.
    args <- term
    args[[1]] <- as.name("list")
  }else args <- list()
  
  termFun<-locate.InitFunction(term, paste0(termroot,"Term"), "ERGM term")  # check in all namespaces for function found anywhere

  # A kludge so that a term can read its own name.
  term.env <- new.env(parent=env) # Create a temporary new environment within env.
  assign(attr(termFun,"fname"), termFun, pos=term.env) # Assign the function body there.
  term<-as.call(list(as.name(attr(termFun,"fname"))))
  
  term[[2]] <- nw
  names(term)[2] <-  "nw"
  term[[3]] <- args
  names(term)[3] <- "arglist"
  dotdotdot <- c(if(!is.null(response)) list(response=response), list(role=role), list(...))
  for(j in seq_along(dotdotdot)) {
    if(is.null(dotdotdot[[j]])) next
    term[[3+j]] <- dotdotdot[[j]]
    names(term)[3+j] <- names(dotdotdot)[j]
  }
  #Call the InitErgm function in the environment where the formula was created
  # so that it will have access to any parameters of the ergm terms
  out <- eval(term,term.env)
  # If SO package name not specified explicitly, autodetect.
  if(!is.null(out) && is.null(out$pkgname)) out$pkgname <- attr(termFun,"pkgname")
  out
}


#######################################################################
# The <updatemodel.ErgmTerm> function updates an existing model object
# to include an initialized ergm term, X;
#
# --PARAMETERS--
#   model  : the pre-existing model, as created by <ergm.getmodel>
#   outlist: the list describing term X, as returned by <InitErgmTerm.X>
#
# --RETURNED--
#   model: the updated model (with the obvious changes seen below) if
#            'outlist'!=NULL, else
#          the original model; (note that this return is necessary,
#            since terms may be eliminated by giving only 0 statistics,
#            and consequently returning a NULL 'outlist')
#
#######################################################################

updatemodel.ErgmTerm <- function(model, outlist) { 
  if (!is.null(outlist)) { # Allow for no change if outlist==NULL
    model$coef.names <- c(model$coef.names, outlist$coef.names)
    termnumber <- 1+length(model$terms)
    tmp <- attr(outlist$inputs, "ParamsBeforeCov")
    # If the term requests auxiliaries or is an auxiliary itself,
    # reserve space in the input vector. Note that these go before
    # the parameters.
    aux.space <-
      if(!is.null(outlist$auxiliaries)) # requests auxiliaries
        length(term.list.formula(outlist$auxiliaries[[length(outlist$auxiliaries)]]))
      else if(length(outlist$coef.names)==0) 1 # is an auxiliary
      else 0
    outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp)+aux.space,
                        length(outlist$coef.names), 
                        length(outlist$inputs)+aux.space, rep(NA, aux.space), outlist$inputs)
    model$minval <- c(model$minval,
                      rep(NVL(outlist$minval, -Inf),
                          length.out=length(outlist$coef.names)))
    model$maxval <- c(model$maxval,
                      rep(NVL(outlist$maxval, +Inf),
                          length.out=length(outlist$coef.names)))
    model$duration <- c(model$duration,
                      NVL(outlist$duration, FALSE))
    model$terms[[termnumber]] <- outlist
  }
  model
}

