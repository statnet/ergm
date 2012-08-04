# For now, this file contains information about reference
# meausres. Eventually, we should create an "InitReference" or similar
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
