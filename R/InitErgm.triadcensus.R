#  See InitErgm.r for a general explanation 
#  of InitErgm functions

InitErgm.triadcensus.directedonly<-function (nw, m, arglist, drop=FALSE, ...) {
  ergm.checkdirected("triadcensus", is.directed(nw), requirement=TRUE)
  a=ergm.checkargs("triadcensus", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  d<-a$d
  tcn <- c("003","012", "102", "021D", "021U", "021C", "111D",
           "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  if(is.null(d)){d <- 2:16}
  if(drop){
    mdegree <- paste("c(",paste(d,collapse=","),")",sep="")
    mdegree <- summary(
     as.formula(paste('nw ~ triadcensus(',mdegree,')',sep="")),
     drop=FALSE) == 0
    if(any(mdegree)){
     cat(paste("Warning: There are no triads of type", tcn[d[mdegree]],".\n"))
     dropterms <- paste("triad type", tcn[d[mdegree]],sep="")
     cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
     d <- d[!mdegree]
    }
  }
  lengthd<-length(d)
  if(lengthd==0){return(model)}
  termnumber<-1+length(m$terms)
# No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="triadcensus", soname="ergm",
                                      inputs=c(0, lengthd, lengthd, d),
                                      dependence=TRUE)
  m$coef.names<-c(m$coef.names, 
   paste("triadcensus",tcn,sep=".")[d])
  m
}
InitErgm.triadcensus<-function (nw, m, arglist, drop=FALSE, ...) {
  a=ergm.checkargs("triadcensus", arglist,
    varnames = c("d","drop"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  d<-a$d
  drop<-a$drop
  detach(a)
  if(is.directed(nw)){
   tcn <- c("003","012", "102", "021D", "021U", "021C", "111D",
            "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
   if(is.null(d)){d <- 1:15}
   if(is.character(d)){d <- match(d, tcn)-1}
  }else{
#  Undirected
   tcn <- c("0", "1", "2", "3")
   if(is.null(d)){d <- 1:3}
   if(is.character(d)){d <- match(d, tcn)-1}
  }
  if(drop){
    mdegree <- paste("c(",paste(d,collapse=","),")",sep="")
    mdegree <- summary(
     as.formula(paste('nw ~ triadcensus(',mdegree,')',sep="")),
     drop=FALSE) == 0
    if(any(mdegree)){
     cat(paste("Warning: There are no triads of type", tcn[d[mdegree]],".\n"))
     dropterms <- paste("triad type", tcn[d[mdegree]],sep="")
     cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
     d <- d[!mdegree]
    }
  }
  d <- d + 1
  lengthd<-length(d)
  if(lengthd==0){return(model)}
  termnumber<-1+length(m$terms)
# No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="triadcensus", soname="ergm",
                                      inputs=c(0, lengthd, lengthd, d),
                                      dependence=TRUE)
  m$coef.names<-c(m$coef.names, paste("triadcensus",tcn,sep=".")[d])
  m
}
InitErgm.balance<-function (nw, m, arglist, drop=TRUE, ...) {
  a=ergm.checkargs("balance", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    mdegree <- summary(
     as.formula('nw ~ balance'),
     drop=FALSE) == 0
    if(mdegree){
     cat(paste("Warning: There are no balanced triads.\n"))
     cat(paste("To avoid degeneracy the balance term has been dropped.\n"))
     return(m)
    }
  }
  termnumber<-1+length(m$terms)
# No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="balance", soname="ergm",
                                      inputs=c(0, 1, 0),
                                      dependence=TRUE)
  m$coef.names<-c(m$coef.names, "balance")
  m
}
