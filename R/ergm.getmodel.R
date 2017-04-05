ergm.getmodel <- function (formula, nw, silent=FALSE, ...) {
  # Parse the formula, create an object of class "model.ergm" that contains
  # all relevant information about the model.  As part of this job, call the
  # appropriate InitErgm functions.
  if ((dc<-data.class(formula)) != "formula")
    stop (paste("Invalid formula of class ",dc), call.=FALSE)
  trms<-terms(formula)
  if (trms[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.", call.=FALSE)  
  v <- attr(trms, "variables")
  if (length(v) < 3) 
    stop(paste("No model specified for network ", trms[[2]]), call.=FALSE)
  model <- structure(list(formula=formula, node.attrib = NULL, coef.names = NULL,
                      offset = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "model.ergm")
  for (i in 3:length(v)) {
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
      newInitErgm <- exists(fname, env=.GlobalEnv, mode="function")
      v[[i]][[1]] <- as.name(ifelse (newInitErgm, fname, 
                                     paste("InitErgm.", v[[i]][[1]], sep = "")))
    } else { # This term has no arguments
      fname <- paste("InitErgmTerm.", v[[i]], sep = "")
      newInitErgm <- exists(fname, env=.GlobalEnv, mode="function")
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
      if(!exists(as.character(v[[i]][[1]]),env=.GlobalEnv, mode="function")){
        stop("The term ", substring(as.character(v[[i]][[1]]),first=10),
             " does not exist. Are you sure you have the right name?\n",
             call. = FALSE)
      }
      if(silent){
       silentwarnings <- capture.output(
        model <- eval(v[[i]], .GlobalEnv)  #Call the InitErgm function
       )
      }else{
       model <- eval(v[[i]], .GlobalEnv)  #Call the InitErgm function
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
      outlist <- eval(v[[i]], .GlobalEnv)  #Call the InitErgm function
      # Now it is necessary to add the output to the model object
      model$coef.names <- c(model$coef.names, outlist$coef.names)
      termnumber <- 1+length(model$terms)
      outlist$inputs <- c(0, length(outlist$coef.names), 
                          length(outlist$inputs), outlist$inputs)
      model$terms[[termnumber]] <- outlist
    }
  } 
  model$etamap <- ergm.etamap(model)
  model
}
