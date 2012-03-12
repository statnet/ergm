# Patched version that does not use .Internal call, courtesy of Pavel.
# Thanks, Pavel!
# The original version, now commented out, is in formula.utils.R
# FIXME:  If the R developers can be convinced to modify update.formula
# so that it does not automatically call "simply=TRUE" on the output,
# then we can use update.formula.

ergm.update.formula<-function (object, new, ...){
  old.lhs <- if(length(object)==2) NULL else object[[2]]
  old.rhs <- if(length(object)==2) object[[2]] else object[[3]]
 
  new.lhs <- if(length(new)==2) NULL else new[[2]]
  new.rhs <- if(length(new)==2) new[[2]] else new[[3]]
 
  sub.dot <- function(c, dot){
    if(is.null(dot)) c # If nothing to substitute with, just return it.
    else if(is.call(c)) as.call(c(list(c[[1]]), lapply(c[-1], sub.dot, dot))) # If it's a call, construct a call consisting of the call and each of the arguments with the substitution performed, recursively.
    else if(is.name(c) && c==".")  dot # If it's a dot, return substitute.
    else c # If it's anything else, just return it.
  }
 
  # Construct the formula and ensure that the formula's environment gets set to the network's
  # environment.
  out <- if(length(new)==2) call("~", sub.dot(new.rhs, old.rhs)) else call("~", sub.dot(new.lhs, old.lhs), sub.dot(new.rhs, old.rhs))
  as.formula(out, env =  if(new[[2]]==".") environment(object) else environment(new))
}



# Previous version of Pavel's patch:
#ergm.update.formula<-function (object, new, ...){
#  old.lhs <- if(length(object)==2) list() else term.list.formula(object[[2]])
#  old.rhs <- if(length(object)==2) term.list.formula(object[[2]]) else term.list.formula(object[[3]])
# 
#  new.lhs <- if(length(new)==2) list() else term.list.formula(new[[2]])
#  new.rhs <- if(length(new)==2) term.list.formula(new[[2]]) else term.list.formula(new[[3]])
#
#  out.lhs <- list()
#  for(term in new.lhs) out.lhs <- c(out.lhs, if(term==".") old.lhs else list(term))
#  out.rhs <- list()
#  for(term in new.rhs) out.rhs <- c(out.rhs, if(term==".") old.rhs else list(term))
#
#  mk.sum <- function(summands){
#    if(length(summands)==0) return(NULL)
#    if(length(summands)==1) return(summands[[1]])
#
#    out <- call("+",summands[[1]],summands[[2]])
#    for(summand in summands[-(1:2)]) out <- call("+", out, summand)
#    out    
#  }
#
#  # Construct the formula and ensure that the formula's environment gets set to the network's
#  # environment.  
#  as.formula(if(length(out.lhs)==0) call("~", mk.sum(out.rhs)) else call("~", mk.sum(out.lhs), mk.sum(out.rhs)),
#             env =  if(new[[2]]==".") environment(object) else environment(new))
# }
 

