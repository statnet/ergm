#  File R/wtd.median.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#################################################################################
# The <wtd.median> function determines and returns the weighted median of a
# vector x
#
# --PARAMETERS--
#   x     : a numeric vector
#   na.rm : whether missing values should be removed (T or F); default=FALSE
#   weight: a vector of weights for each value in x
#
# --RETURNED--
#   NA                     if x has missing values and 'na.rm'=FALSE
#   the weighted median    otherwise  
#
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
