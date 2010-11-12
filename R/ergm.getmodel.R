ergm.getmodel <- function (formula, nw, response=NULL, silent=FALSE, ...,stergm.order=NULL) {
  # Parse the formula, create an object of class "model.ergm" that contains
  # all relevant information about the model.  As part of this job, call the
  # appropriate InitErgm functions.
  if ((dc<-data.class(formula)) != "formula")
    stop (paste("Invalid formula of class ",dc), call.=FALSE)

  if (formula[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.", call.=FALSE)  
  
  if (length(formula) < 3) 
    stop(paste("No model specified for network ", formula[[2]]), call.=FALSE)

  v<-term.list.formula(formula[[3]])
  
  formula.env<-environment(formula)
  
  model <- structure(list(formula=formula, node.attrib = NULL, coef.names = NULL,
                      offset = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "model.ergm")

  termroot<-if(is.null(response)) "InitErgm" else "InitWtErgm"

  
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
      newInitErgm <- exists(fname, env=formula.env, mode="function")
      v[[i]][[1]] <- as.name(ifelse (newInitErgm, fname, 
                                     paste(termroot,".", v[[i]][[1]], sep = "")))
    } else { # This term has no arguments
      fname <- paste(termroot,"Term.", v[[i]], sep = "")
      newInitErgm <- exists(fname, env=formula.env, mode="function")
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
      dotdotdot <- c(list(response=response), list(...))
      for(j in seq_along(dotdotdot)) {
        if(is.null(dotdotdot[[j]])) next
        v[[i]][[4+j]] <- dotdotdot[[j]]
        names(v[[i]])[4+j] <- names(dotdotdot)[j]
      }    
      # The above steps are preparing the way to make the function call
      # InitErgm.xxxx(g, m, args, ...)
      if(!exists(as.character(v[[i]][[1]]),env=formula.env, mode="function")){
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
      dotdotdot <- c(list(response=response),list(...))
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
  model$stergm.order <- stergm.order
  model
}

# Take the output of an InitErgmTerm.xxx function and add it correctly
# to an existing model object.  If outlist is NULL, then simply return
# original model object.  This is sometimes important, if for example
# a term is to be eliminated because it gives only zero statistics.
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

