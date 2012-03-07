#  File ergm/R/wtd.median.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#################################################################################
# The <wtd.median> function determines and returns the weighted median of a
# vector x
################################################################################

wtd.median <- function(x, na.rm = FALSE, weight=FALSE) {
 	if(mode(x) != "numeric")
 		stop("need numeric data")
 	x <- as.vector(x)
 	wnas <- is.na(x)
 	if(sum(wnas)>0) {
 		if(na.rm)
 		 x <- x[!wnas]
 	  	 if(!missing(weight)){weight <- weight[!wnas]}
 		else return(NA)
 	}
 	n <- length(x)
 	half <- (n + 1)/2
 	if(n %% 2 == 1) {
 	  if(!missing(weight)){
 		weight <- weight/sum(weight)
 		sx <- sort.list(x)
 		sweight <- cumsum(weight[sx])
 		min(x[sx][sweight >= 0.5])
 	  }else{
 		x[order(x)[half]]
 	  }
 	}
 	else {
 	  if(!missing(weight)){
 		weight <- weight/sum(weight)
 		sx <- sort.list(x)
 		sweight <- cumsum(weight[sx])
 		min(x[sx][sweight >= 0.5])
 	  }else{
 		half <- floor(half) + 0:1
 		sum(x[order(x)[half]])/2
 	  }
 	}
 }

# Got rid of wtd.mean function because weighted.mean already exists in R
