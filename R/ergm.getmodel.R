#  File ergm/R/ergm.getmodel.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###################################################################################
# The <ergm.getmodel> function parses the given formula, and initiliazes each ergm
# term via the <InitErgmTerm> functions to create a 'model.ergm' object for the
# given network
#####################################################################################

ergm.getmodel <- function (formula, nw, silent=FALSE, role="static",...) {
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

  termroot<- "InitErgm"

  
  for (i in 1:length(v)) {
    if (is.call(v[[i]])) { # This term has some arguments
      if(v[[i]][[1]] == "offset"){
        if(length(v[[i]][[2]]) <= 1){
         v[[i]] <- as.call(v[[i]][2])
        }else{
         v[[i]] <- as.call(v[[i]][[2]])
        }
        model$offset <- c(model$offset,TRUE)
      }else{
        model$offset <- c(model$offset,FALSE)
      }
      args=v[[i]]
      args[[1]] = as.name("list")
      fname <- paste(termroot,"Term.", v[[i]][[1]], sep = "")
      newInitErgm <- exists(fname, envir=formula.env, mode="function")
      v[[i]] <- call(ifelse (newInitErgm, fname, 
                             paste(termroot,".", v[[i]][[1]], sep = "")))
    } else { # This term has no arguments
      fname <- paste(termroot,"Term.", v[[i]], sep = "")
      newInitErgm <- exists(fname, envir=formula.env, mode="function")
      v[[i]] <- call(ifelse (newInitErgm, fname, 
                             paste(termroot,".", v[[i]], sep = "")))
      model$offset <- c(model$offset,FALSE)
      args=list()
    }
    if (!newInitErgm) { #Using the old InitErgm style
      v[[i]][[2]] <- nw
      names(v[[i]])[2] <-  ""
      v[[i]][[3]] <- model
      names(v[[i]])[3] <- ""
      v[[i]][[4]] <- args
      dotdotdot <- c(list(role=role), list(...))
      for(j in seq_along(dotdotdot)) {
        if(is.null(dotdotdot[[j]])) next
        v[[i]][[4+j]] <- dotdotdot[[j]]
        names(v[[i]])[4+j] <- names(dotdotdot)[j]
      }    
      # The above steps are preparing the way to make the function call
      # InitErgm.xxxx(g, m, args, ...)
      if(!exists(as.character(v[[i]][[1]]),envir=formula.env, mode="function")){
        stop("The term ", substring(as.character(v[[i]][[1]]),first=nchar(termroot)+2),
             " does not exist for this type of ERGM. Are you sure you have the right name?\n",
             call. = FALSE)
      }
      if(silent){
       silentwarnings <- capture.output(
        model <- eval(v[[i]], formula.env)  #Call the InitErgm function
       )
      }else{
       model <- eval(v[[i]], formula.env)  #Call the InitErgm function
      }
    } else { # New InitErgmTerms style
      v[[i]][[2]] <- nw
      names(v[[i]])[2] <-  ""
      v[[i]][[3]] <- args
      names(v[[i]])[3] <- ""
      dotdotdot <- c(list(role=role), list(...))
      for(j in seq_along(dotdotdot)) {
        if(is.null(dotdotdot[[j]])) next
        v[[i]][[3+j]] <- dotdotdot[[j]]
        names(v[[i]])[3+j] <- names(dotdotdot)[j]
      }
      outlist <- eval(v[[i]], formula.env)  #Call the InitErgm function
      # Now it is necessary to add the output to the model object
      model <- updatemodel.ErgmTerm(model, outlist)
    }
  } 
  model$etamap <- ergm.etamap(model)
  model
}



#######################################################################
# The <updatemodel.ErgmTerm> function updates an existing model object
# to include an initialized ergm term, X;
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
    model$terms[[termnumber]] <- outlist
  }
  model
}

