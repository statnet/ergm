 "wtd.median"<- function(x, na.rm = FALSE, weight=FALSE) {
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
