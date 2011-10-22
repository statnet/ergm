#  See InitErgm.R for a general explanation 
#  of InitErgm functions
InitErgm.dissolve<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("dissolve", arglist=arglist,
    varnames = c("x"),
    vartypes = c("matrix,network"),
    defaultvalues = list(nw),
    required = c(FALSE))
  attach(a)
  x<-a$x
#
# Extract reference network as an edgelist
#
  nwm<-as.matrix.network(nw,matrix.type="edgelist")
  if (is.null(nwm) || ncol(nwm)!=2){
    stop("dissolve() requires a proper network as its reference")
  }
  termnumber<-1+length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "dissolve", soname="ergm",
                                inputs = c(1, 1, 1+length(nwm),
                                   nrow(nwm), as.integer(nwm)),
                                 dependence=FALSE)
  cn<-paste("dissolve", as.character(sys.call(0)[[4]][2]), sep = ".")
  m$coef.names <- c(m$coef.names, cn)
  m
}

InitErgm.formation<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("formation", arglist=arglist,
    varnames = c("x"),
    vartypes = c("matrix,network"),
    defaultvalues = list(nw),
    required = c(FALSE))
  attach(a)
  x<-a$x
#
# Extract reference network as an edgelist
#
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist")
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else if(is.null(x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    stop("formation() requires a proper network as its reference")
  }
  termnumber<-1+length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "formation", soname="ergm",
                                 inputs = c(1, 1,
                                   1+2*nrow(xm),
                                   nrow(xm), as.integer(xm)),
                                 dependence=TRUE)
   cn<-paste("formation", as.character(sys.call(0)[[4]][2]), sep = ".")
   m$coef.names <- c(m$coef.names, cn)
   m
}
