
InitErgmTerm..spcache.netL<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("type", "Ls.path", "L.in_order"),
                      vartypes = c("character", "formula,list", "logical"),
                      defaultvalues = list(NULL,NULL,NULL),
                      required = c(TRUE, TRUE, TRUE))

  type <- match.arg(tolower(a$type), c("otp","osp","isp","utp")) # ITP not included, because it's just OTP with direction reversed.

  if(is.directed(nw)==(type=="utp")) stop("Type UTP may only be used with undirected networks, the others only with directed.")

  dname <- paste0("_",type,"_wtnet")
  
  linfo <- .sp.handle_layers(nw, a, type, FALSE)

  
  list(name=paste0(dname,linfo$name_suffix), auxiliaries=linfo$auxiliaries, inputs=c(linfo$any_order),
       coef.names=c(), dependence=TRUE)
}
