#  File R/wtd.median.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

#' Weighted Median
#' 
#' Compute weighted median.
#' 
#' Uses a simple algorithm based on sorting.
#' 
#' @param x Vector of data, same length as \code{weight}
#' @param na.rm Logical: Should NAs be stripped before computation proceeds?
#' @param weight Vector of weights
#' @return Returns an empirical .5 quantile from a weighted sample.
#' @keywords robust
#' @export
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
