# This new InitErgmTerm function still needs to be tested:

#################################################################################
InitErgmTerm.transitiveties<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  if (a$diff) stop("diff=TRUE is not currently implemented in transitiveties")
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "transitiveties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to transitiveties() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ transitiveties(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("transitiveties",attrname,sep="."),
                             u[triattr],sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("transitiveties",attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        }
      }
    }
    if (!diff) {
      coef.names <- paste("transitiveties",attrname,sep=".")
      inputs <- c(nodecov)
     } else { 
       coef.names <- paste("transitiveties",attrname, u, sep=".")
       inputs <- c(ui, nodecov)
       attr(inputs, "ParamsBeforeCov") <- length(ui)
     }
  }else{
    coef.names <- "transitiveties"
    inputs <- NULL
  }
  list(name="transitiveties", coef.names=coef.names, inputs=inputs)
}

