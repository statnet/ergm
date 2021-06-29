#  File R/summary.network.list.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################


#' @describeIn network.list A [summary()] method for network lists.
#' 
#' @param stats.print Logical: If TRUE, print network statistics.
#' @param net.print Logical: If TRUE, print network overviews.
#' @param net.summary Logical: If TRUE, print network summaries.
#' @seealso \code{\link{simulate.ergm}}
#' @examples
#' 
#' # Draw from a Bernoulli model with 16 nodes
#' # and tie probability 0.1
#' #
#' g.use <- network(16, density=0.1, directed=FALSE)
#' #
#' # Starting from this network let's draw 3 realizations
#' # of a model with edges and 2-star terms
#' #
#' g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8, 0.03),
#'                basis=g.use, control=control.simulate(
#'                  MCMC.burnin=100000,
#'                  MCMC.interval=1000))
#' print(g.sim)
#' summary(g.sim)
#' 
#' @export
summary.network.list <- function (object, stats.print=TRUE, 
                       net.print=FALSE, net.summary=FALSE, ...){

  cat("Number of Networks:",length(object),"\n")
  attrmap<-list(formula="Model: ",
                reference="Reference: ",
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
  object
}


