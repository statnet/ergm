##### edgecov for ranks is implemented in InitWtErgmTerm.R #####

InitWtErgmTerm.deference<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)

  list(name="deference",
       coef.names="deference",
       inputs=NULL,
       dependence=TRUE)
}

InitWtErgmTerm.inconsistency<-function (nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                     varnames = c("x","attrname"),
                     vartypes = c("matrixnetwork","character"),
                     defaultvalues = list(nw,NULL),
                     required = c(FALSE,FALSE))

    ## Process hamming network ##
  if(is.network(a$x)){ # Arg to hamming is a network
    xm<-as.matrix.network(a$x,matrix.type="adjacency",a$attrname)
  }else if(is.character(a$x)){ # Arg to hamming is the name of an attribute in nw
    xm<-get.network.attribute(nw,a$x)
    xm<-as.matrix.network(xm,matrix.type="adjacency")
  }else{
    xm<-as.matrix(a$x) # Arg to hamming is anything else; attempts to coerce
  }
  ## Process case without dyadcov (i.e. unweighted) ##
  sc03 <- sys.call(0)[[3]]
  coef.names <- "inconsistency"  # This might be modified later
  if (length(sc03)>1) 
    coef.names <- paste("inconsistency", as.character(sc03[[2]]), sep=".")
  
  if(!is.null(a$attrname) && length(sc03)>1){
    coef.names<-paste("inconsistency", as.character(sc03[2]),
                      as.character(a$attrname), sep = ".")
  }else if (length(sc03)>1) {
    coef.names<-paste("inconsistency", as.character(sc03[2]),
                      as.character(sys.call(0)[[3]][3]), sep = ".")
  }

  # A column-major matrix of choices.
  inputs <- c(t(xm))

  list(name="inconsistency_rank", coef.names=coef.names, #name and coef.names: required 
       inputs = inputs, dependence = FALSE)
}

##### nodeicov for ranks is implemented in InitWtErgmTerm.R #####

InitWtErgmTerm.nonconformity<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("form","par"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("all",NULL),
                      required = c(FALSE,FALSE))

  form<-match.arg(a$form,c("all","thresholds","geometric","local"))

  if(form=="all"){
    inputs <- NULL
    coef.names <- "nonconformity"
    name <- "nonconformity"
  }else if(form=="thresholds"){
    inputs <- sort(as.numeric(a$par),decreasing=TRUE)
    coef.names <- paste("nonconformity.over",inputs,sep=".")
    name <- "nonconformity_thresholds"
  }else if(form=="geometric"){
    inputs <- c(a$par,max(nw %e% response))
    coef.names <- paste("nonconformity.gw.",a$par,sep=".")
    name <- "nonconformity_decay"
  }else if(form=="local"){
    inputs <- NULL
    coef.names <- "nonconformity.local"
    name <- "local_nonconformity"
  }
  
  list(name=name,
       coef.names=coef.names,
       inputs=inputs,
       dependence=TRUE)
}
