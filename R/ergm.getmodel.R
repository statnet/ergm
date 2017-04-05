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
#               drop      : whether to drop degenerate terms (T or F)
#               initialfit: whether curved exponential terms have been initially fit
#                           by MPLE (T or F)
#
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

ergm.getmodel <- function (formula, nw, silent=FALSE, ...) {
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
      fname <- paste("InitErgmTerm.", v[[i]][[1]], sep = "")
      newInitErgm <- exists(fname, env=formula.env, mode="function")
      v[[i]][[1]] <- as.name(ifelse (newInitErgm, fname, 
                                     paste("InitErgm.", v[[i]][[1]], sep = "")))
    } else { # This term has no arguments
      fname <- paste("InitErgmTerm.", v[[i]], sep = "")
      newInitErgm <- exists(fname, env=formula.env, mode="function")
      v[[i]] <- call(ifelse (newInitErgm, fname, 
                             paste("InitErgm.", v[[i]], sep = "")))
      model$offset <- c(model$offset,FALSE)
      args=list()
    }
    if (!newInitErgm) { #Using the old InitErgm style
      v[[i]][[2]] <- nw
      names(v[[i]])[2] <-  ""
      v[[i]][[3]] <- model
      names(v[[i]])[3] <- ""
      v[[i]][[4]] <- args
      dotdotdot <- list(...)
      if (length(dotdotdot>0)) {
        for(j in 1:length(dotdotdot)) {
          v[[i]][[4+j]] <- dotdotdot[[j]]
          names(v[[i]])[4+j] <- names(dotdotdot)[j]
        }
      }
      # The above steps are preparing the way to make the function call
      # InitErgm.xxxx(g, m, args, ...)
      if(!exists(as.character(v[[i]][[1]]),env=formula.env, mode="function")){
        stop("The term ", substring(as.character(v[[i]][[1]]),first=10),
             " does not exist. Are you sure you have the right name?\n",
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
      dotdotdot <- list(...)
      if (length(dotdotdot>0)) {
        for(j in 1:length(dotdotdot)) {
          v[[i]][[3+j]] <- dotdotdot[[j]]
          names(v[[i]])[3+j] <- names(dotdotdot)[j]
        }
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
    model$terms[[termnumber]] <- outlist
  }
  model
}

