ergm.getmodel.dissolve <- function (formula, nw, ...) 
{
  # Parse the formula, create an object of class "model.ergm" that contains
  # all relevant information about the model.  As part of this job, call the
  # appropriate InitErgm functions.
  if(is.null(formula)){
   model <- list(terms=NULL)
   return(model)
  }
  if ((dc<-data.class(formula)) != "formula")
    stop (paste("Invalid formula of class ",dc))
  trms<-terms(formula)
  if (trms[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.")
  
  v <- attr(trms, "variables")
  if (length(v) == 1){
    formula <- update(formula, nw ~ .)
    trms<-terms(formula)
    v <- attr(trms, "variables")
  }
# if (length(v) == 2){
    formula <- update(formula, ~ . + dissolve)
    trms<-terms(formula)
    v <- attr(trms, "variables")
# }
  nw <- ergm.getnetwork(formula)
  model <- structure(list(formula=formula, node.attrib = NULL,
                      coef.names = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "model.ergm")
  for (i in 3:length(v)) {
    if (is.call(v[[i]])) { # This term has some arguments
      args=v[[i]]
      args[[1]] = as.name("list")
      v[[i]][[1]] <- as.name(paste("InitErgm.", v[[i]][[1]],
                                   sep = ""))
    } else { # This term has no arguments
      v[[i]] <- call(paste("InitErgm.", v[[i]], sep = ""))
      args=list()
    }
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
#    m <- eval(v[[i]], sys.parent())  #Call the InitErgm function
    if(!exists(as.character(v[[i]][[1]]),env=.GlobalEnv, mode="function")){
     stop("The term ", substring(as.character(v[[i]][[1]]),first=10),
          " does not exist. Are you sure you have the right name?\n")
    }
    model <- eval(v[[i]], .GlobalEnv)  #Call the InitErgm function
  }
  model$etamap <- ergm.etamap(model)
  model
}
