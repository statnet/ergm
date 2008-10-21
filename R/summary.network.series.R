#  File ergm/R/summary.network.series.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
"summary.network.series" <-
  function (object, ..., wmt = which.matrix.type(objects$networks[[1]])) 
{
  g <- object$networks[[1]]
  cat("Number of Networks:",length(object$networks),"\n")
  f.out <- lapply(object$formula,as.character)
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
  cat("parameter values:",object$coef,"\n")

  summary(g)
}
