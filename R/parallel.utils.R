# Makes a heroic effort to get a cluster of desired size.
ergm.getCluster <- function(control, verbose=FALSE){
  capture.output(require(snow, quietly=TRUE, warn.conflicts = FALSE))
#
# Start PVM if necessary
#
  if(getClusterOption("type")=="PVM"){
    if(verbose){cat("Engaging warp drive using PVM ...\n")}
    capture.output(require(rpvm, quietly=TRUE, warn.conflicts = FALSE))
    PVM.running <- try(.PVM.config(), silent=TRUE)
    if(inherits(PVM.running,"try-error")){
      hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by ergm...\n")
    }
  }else{
    if(verbose){cat("Engaging warp drive using MPI ...\n")}
  }
  #
  #   Start Cluster
  #
  cl<-makeCluster(control$parallel)
  clusterSetupRNG(cl)
  if("ergm" %in% control$packagenames){
    clusterEvalQ(cl,library(ergm))
  }
  cl
}
