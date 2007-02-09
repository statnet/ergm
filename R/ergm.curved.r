"ergm.curved" <- function(theta0,m,m.expanded,nw,theta1=NULL){
  geodf <- rep(FALSE, length(m.expanded$coef.names))
  geod <- grep("geodegree#",m.expanded$coef.names)     
  if(length(geod)>0){
    geodf[geod] <- TRUE
  }
  geosdf <- rep(FALSE, length(m.expanded$coef.names))
  geosd <- grep("geospartner#",m.expanded$coef.names)     
  if(length(geosd)>0){
    geosdf[geosd] <- TRUE
  }
  wdegf <- rep(FALSE, length(m.expanded$coef.names))
  wdeg <- grep("gwdegree#",m.expanded$coef.names)     
  if(length(wdeg)>0){
    wdegf[wdeg] <- TRUE
  }
  gwespf <- rep(FALSE, length(m.expanded$coef.names))
  gwesp <- grep("gwesp#",m.expanded$coef.names)     
  if(length(gwesp)>0){
    gwespf[gwesp] <- TRUE
  }
#
#     Extend for curved parameters
#
  if(any(geodf) | any(geosdf) | any(wdegf) | any(gwespf) ){
    nongeonames <- is.na(match(names(theta0),
                               c("gd.alpha", "gsd.alpha", "geodegree",
                                 "geospartner", "gwdegree.alpha",
                                 "gwesp.alpha","gwdegree", "gwesp")))
    eta0 <- rep(0,length=length(m.expanded$coef.names))
    names(eta0) <- m.expanded$coef.names
  }else{
    eta0 <- theta0
  }
  parms.curved <- vector(mode="list")
  if(any(geodf)){
    loc.geodegree <- match("geodegree",names(theta0))
    for(i in 1:length(m$terms)){
      if(m$terms[[i]]$name=="geodegree"){
        parms.curved$geodegree.alpha <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$geodegree.alpha) <- "gd.alpha"
    theta0 <- c(theta0, parms.curved$geodegree.alpha)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-1))
    eta0[geodf] <- theta0["geodegree"]*(exp(-theta0["gd.alpha"]*iseq)-1)
  }
  if(any(geosdf)){
    loc.geospartner <- match("geospartner",names(theta0))
    for(i in 1:length(m$terms)){
      if(m$terms[[i]]$name=="geospartner"){
        parms.curved$geospartner.alpha <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$geospartner.alpha) <- "gsd.alpha"
    theta0 <- c(theta0, parms.curved$geospartner.alpha)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-2))
    eta0[geosdf] <- theta0["geospartner"]*(exp(-theta0["gsd.alpha"]*iseq)-1)
  }
  if(any(wdegf)){
    loc.gwdegree <- match("gwdegree",names(theta0))
    for(i in 1:length(m$terms)){
      if(m$terms[[i]]$name=="gwdegree"){
        parms.curved$gwdegree.alpha <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$gwdegree.alpha) <- "gwdegree.alpha"
    theta0 <- c(theta0, parms.curved$gwdegree.alpha)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-1))
    eta0[wdegf] <- theta0["gwdegree"]*
      (iseq*exp(theta0["gwdegree.alpha"]) +
       exp(2*theta0["gwdegree.alpha"])*
       (((1-exp(-theta0["gwdegree.alpha"]))^iseq)-1))
  }
  if(any(gwespf)){
    loc.gwesp <- match("gwesp",names(theta0))
    for(i in 1:length(m$terms)){
      if(m$terms[[i]]$name=="gwesp"){
        parms.curved$gwesp.alpha <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$gwesp.alpha) <- "gwesp.alpha"
    theta0 <- c(theta0, parms.curved$gwesp.alpha)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-2))
    eta0[gwespf] <- theta0["gwesp"]*
      exp(theta0["gwesp.alpha"])*(1-(1-exp(-theta0["gwesp.alpha"]))^iseq)
  }
  if(any(geodf) | any(geosdf) | any(wdegf) | any(gwespf) ){
    nongeonames <- is.na(match(names(theta0),
                               c("gd.alpha", "gsd.alpha", "geodegree",
                                 "geospartner", "gwdegree.alpha",
                                 "gwesp.alpha","gwdegree", "gwesp")))
    eta0[!geodf & !geosdf & !wdegf & !gwespf] <- theta0[nongeonames]
    parms.curved$eta0 <-  eta0
  }else{
    parms.curved$eta0 <-  theta0
  }
  list(parms.curved=parms.curved,eta0=eta0,theta0=theta0,theta1=theta1,
       geodf=geodf,geosdf=geosdf,wdegf=wdegf,gwespf=gwespf)
}
