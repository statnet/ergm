
InitErgmTerm..spcache.net<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("type"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  type <- match.arg(tolower(a$type), c("otp","osp","isp","utp")) # ITP not included, because it's just OTP with direction reversed.

  if(is.directed(nw)==(type=="utp")) stop("Type UTP may only be used with undirected networks, the others only with directed.")
  
  list(name=paste0("_",type,"_wtnet"),
       coef.names=c(), dependence=TRUE)
}
