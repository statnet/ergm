#  File R/WtReferences.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
# For now, this file contains information about reference
# measures. Eventually, we should create an "InitReference" or similar
# framework.

ergm.init.methods <- local({
  init.methods <- list()
  function(reference, new.methods){
    if(!missing(new.methods)){
      init.methods[[reference]] <<- unique(c(new.methods, init.methods[[reference]]))
    }else{
      init.methods[[reference]]
    }
  }
})
