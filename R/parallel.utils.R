# Makes a heroic effort to get a cluster of desired size.
ergm.getCluster <- function(control, verbose=FALSE){
  capture.output(require(snow, quietly=TRUE, warn.conflicts = FALSE))

  reqs<-list("PVM"="rpvm", "MPI"="Rmpi")
    
  for(type in unique(c(getClusterOption("type"),"PVM","MPI","SOCK"))){
    if(verbose) cat("Attempting ",type,"... ", sep="")
    #   Start Cluster


    if(type %in% names(reqs)){
      if(! reqs[[type]] %in% rownames(installed.packages())){
        if(verbose) cat("Failed.\n")
        next
      }else{
        if(!require(reqs[[type]], quietly=TRUE, warn.conflicts = FALSE, character.only=TRUE))
          next
      }
    }
    
    if(type=="PVM")
      try({
        PVM.running <- try(.PVM.config(), silent=TRUE)
        if(inherits(PVM.running,"try-error")){
          hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
          .PVM.start.pvmd(hostfile)
          cat("PVM not running. Attempting to start.\n")
        }},silent=TRUE) # OK if fails.
    
    
    capture.output({
      cl<-try(makeCluster(control$parallel,type=type),silent=TRUE)
    })
    
    if(inherits(cl,"try-error")){
      if(verbose) cat("Failed.\n")
    }else{
      if(verbose) cat("Success!\n")
      break
    }
    
  }

  # Set things up. 
  clusterSetupRNG(cl)
  if("ergm" %in% control$packagenames){
    clusterEvalQ(cl,library(ergm))
  }
  cl
}
