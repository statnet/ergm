#  File ergm/R/summary.network.list.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
summary.network.list <- function (object, stats.print=TRUE, 
                       net.print=FALSE, net.summary=FALSE, ...){

  cat("Number of Networks:",length(object),"\n")
  attrmap<-list(formula="Model: ",
                constraints="Constraints: ",
                coef="Parameters:\n",
                stats="Stored network statistics:\n")
  if(!stats.print) attrmap$stats <- NULL
  for(a in names(attrmap)){
    if(a %in% names(attributes(object))){
      s<-attrmap[[a]]
      cat(s)
      if(substr(s,nchar(s),nchar(s))=="\n"){
        print(attr(object,a))
        cat("\n")
      }else{
        cat(format(attr(object,a)),"\n")
      }
    }
  }
  if(net.print){
    cat("Network overviews:\n")
    o <- object
    attributes(o)<-list()
    print(o)
  }
  if(net.summary){
    cat("Network summaries:\n")
    print(lapply(object,summary,...))
  }
}


