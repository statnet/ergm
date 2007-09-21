"update.ergm" <- function(object, ..., 
            contriburl = "http://csde.washington.edu/ergm",
            repos = getOption("repos"), type = getOption("pkgType"))
{
  cran.contriburl <- contrib.url(repos, type)
  cat(paste("Initializing...\n",sep=""))
  cran.Base <- c("network", "coda")
  csde.Base <- c("ergm")
  cran.Recommended <- c("sna", "latentnet")
  csde.Recommended <- c()
# cran.Optional <- c("netdata")
  cran.Optional <- c()
  csde.Optional <- c("netperm", "dnet")
  cran.Base.latentnet <- c("locfit", "mvtnorm", "mclust", 
                           "akima",  "shapes", "KernSmooth")
#
  if(missing(object)){
    object <- c(cran.Base, csde.Base, cran.Recommended,
                csde.Recommended, cran.Optional, csde.Optional)
  }
  inuse <- match(paste("package:","ergm",sep=""), search())
  if(!inherits(try(detach(pos=inuse),silent=TRUE), "try-error")){
    cat(paste("Detaching package '", "ergm","'.\n",sep=""))
  }
#
# local install functions
#
  ergm.install <- function(object, csde, contriburl=cran.contriburl,
                              ask=FALSE, type="recommended",
                              update.pkgs=update.cran.pkgs){
    pma <- match(csde, object)
    csde <- object[pma[!is.na(pma)]]
    for(pkg in csde){
      cat(paste("=========================================\n",sep=""))
      if(pkg %in% update.pkgs){
       cat(paste("Attempting to update package '", pkg,"'.\n",sep=""))
       inuse <- match(paste("package:",pkg,sep=""), search())
       if(!is.na(inuse)){
         if(!inherits(try(detach(pos=inuse),silent=TRUE), "try-error")){
           cat(paste("Detaching package '", pkg,"'.\n",sep=""))}
       }
       if(ask){
        if(type == "recommended"){
         cat(paste("It is recommended you update '",pkg,"'.\n", sep=""))
        }
        if(type == "optional"){
         cat(paste("It is optional to update '",pkg,"'.\n",sep=""))
        }
        answer<-substr(readline(paste("Update ",type," package '",pkg,"' (y/n)?  ", sep="")),1,1)
        if (answer == "y" | answer == "Y"){
         install.packages(pkg, lib=.libPaths()[1], contriburl=contriburl)
        }
       }else{
        install.packages(pkg, lib=.libPaths()[1], contriburl=contriburl)
       }
      }else{
       cat(paste("Package '", pkg,"' is now up-to-date.\n",sep=""))
      }
    }
  invisible()
  }
# CSDE
  new.csde.pkgs <- new.packages(lib.loc=.libPaths()[1],contriburl=contriburl)
  old.csde.pkgs <- old.packages(lib.loc=.libPaths()[1],contriburl=contriburl)[,1]
  update.csde.pkgs <- c(new.csde.pkgs, old.csde.pkgs)
# CRAN
  new.cran.pkgs <- new.packages(lib.loc=.libPaths()[1])
  old.cran.pkgs <- old.packages(lib.loc=.libPaths()[1])[,1]
  update.cran.pkgs <- c(new.cran.pkgs, old.cran.pkgs)
#
  ergm.install(object, cran.Base, ask=FALSE, type="base")
  ergm.install(object, csde.Base, ask=FALSE, type="base",
                  update.pkgs=update.csde.pkgs, contriburl=contriburl)
  ergm.install(object, cran.Recommended, ask=TRUE, type="recommended")
  ergm.install(object, csde.Recommended, ask=TRUE, type="recommended",
                  update.pkgs=update.csde.pkgs, contriburl=contriburl)
  ergm.install(object, cran.Optional, ask=TRUE, type="optional")
  ergm.install(object, csde.Optional, ask=TRUE, type="optional",
                  update.pkgs=update.csde.pkgs, contriburl=contriburl)
#
# check required packages for latentnet
#
  inst.pkgs <- as.vector(installed.packages()[,1])
  if("latentnet" %in% inst.pkgs){
    ergm.install(cran.Base.latentnet, cran.Base.latentnet, ask=FALSE, type="base")
  }
#
  cat(paste("=========================================\n",sep=""))
  cat(paste("'ergm' is now up-to-date.\n",sep=""))
  require(ergm)
  invisible()
}
