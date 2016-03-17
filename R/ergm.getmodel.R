#  File R/ergm.getmodel.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
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

ergm.getmodel <- function (formula, nw, response=NULL, silent=FALSE, role="static",...) {
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

  termroot<-if(is.null(response)) "InitErgm" else "InitWtErgm"

  
  for (i in 1:length(v)) {
    if (is.call(v[[i]]) && v[[i]][[1]] == "offset"){ # Offset term
      v[[i]] <- v[[i]][[2]]
      model$offset <- c(model$offset,TRUE)
    }else{
      model$offset <- c(model$offset,FALSE)
    }
    ## v[[i]] is now a call or a name that is not "offset".
    
    if(is.call(v[[i]])) { # This term has some arguments; save them.
      args <- v[[i]]
      args[[1]] <- as.name("list")
    }else args <- list()
    
    termFun<-locate.InitFunction(v[[i]], paste0(termroot,"Term"), "ERGM term")  # check in all namespaces for function found anywhere
    
    v[[i]]<-as.call(list(termFun))
    
    v[[i]][[2]] <- nw
    names(v[[i]])[2] <-  ""
    v[[i]][[3]] <- args
    names(v[[i]])[3] <- ""
    dotdotdot <- c(if(!is.null(response)) list(response=response), list(role=role), list(...))
    for(j in seq_along(dotdotdot)) {
      if(is.null(dotdotdot[[j]])) next
      v[[i]][[3+j]] <- dotdotdot[[j]]
      names(v[[i]])[3+j] <- names(dotdotdot)[j]
    }
    #Call the InitErgm function in the environment where the formula was created
    # so that it will have access to any parameters of the ergm terms
    outlist <- eval(v[[i]],formula.env)
    # If initialization fails without error (e.g., all statistics have been dropped), continue.
    if(is.null(outlist)) next
    # If SO package name not specified explicitly, autodetect.
    if(is.null(outlist$pkgname)) outlist$pkgname <- environmentName(environment(termFun))
    # Now it is necessary to add the output to the model object
    model <- updatemodel.ErgmTerm(model, outlist)
  } 
  model$etamap <- ergm.etamap(model)

  # I.e., construct a vector of package names associated with the model terms.
  # Note that soname is not the same, since it's not guaranteed to be a loadable package.
  ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
  
  class(model) <- "ergm.model"
  model
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
    outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp),
                        length(outlist$coef.names), 
                        length(outlist$inputs), outlist$inputs)
    model$minval <- c(model$minval,
                      rep(if(!is.null(outlist$minval)) outlist$minval else -Inf,
                          length.out=length(outlist$coef.names)))
    model$maxval <- c(model$maxval,
                      rep(if(!is.null(outlist$maxval)) outlist$maxval else +Inf,
                          length.out=length(outlist$coef.names)))
    model$duration <- c(model$duration,
                      if(!is.null(outlist$duration)) outlist$duration else FALSE)
    model$terms[[termnumber]] <- outlist
  }
  model
}

