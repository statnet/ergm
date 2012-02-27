###############################################################################
# The <summary.network.list> function provides a summary of the first
# network in the given list
#
# --PARAMETERS--
#   object: a network.list object
#
# --IGNORED PARAMETERS--
#   ...:  additional parameters passed from within
#   wmt:  which matrix type is used to describe the networks;
#         default=which.matrix.type(objects$networks[[1]])) 
#
# --RETURNED--
#   the summary of the first network in the list    
# 
###############################################################################

summary.network.list <-
  function (object, ..., 
  wmt = which.matrix.type(g)) 
{
  if(is.null(object$form)) { # NOTE:  "$form" here could be either "formula" or "formation".
                             # This is very sloppy and should be changed after the
                             # networkDynamic objects are working.
    a <- attributes(object)
    g <- object[[1]]
  } else {
    a <- object
    g <- object$networks
  }
  cat("Number of Networks:",length(object),"\n")
  f.out <- lapply(a$form,as.character)
  if(length(f.out)>2)
    {
      if(length(f.out[[3]])>2)
        f.output <- paste(f.out[[2]],f.out[[1]],paste(f.out[[3]][c(2,1,3)],collapse=" "))
      else
        f.output <- paste(f.out[[2]],f.out[[1]],f.out[[3]])
    }
  else
    {
      if(length(f.out[[2]])>2)
        f.output <- paste(f.out[[1]],paste(f.out[[2]][c(2,1,3)],collapse=" "))
      else
        f.output <- paste(f.out[[1]],f.out[[2]])
    }
  cat("From:",f.output,"\n")
  cat("parameter values:",a$coef,"\n")

  summary(g)
}
