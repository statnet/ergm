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
  wodegf <- rep(FALSE, length(m.expanded$coef.names))
  wodeg <- grep("gwodegree#",m.expanded$coef.names)     
  if(length(wodeg)>0){
    wodegf[wodeg] <- TRUE
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
                                 "geospartner", "gwdegree.decay",
                                 "gwodegree.decay",
                                 "gwesp.alpha","gwdegree","gwodegree",
                                 "gwesp")))
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
        parms.curved$gwdegree.decay <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$gwdegree.decay) <- "gwdegree.decay"
    theta0 <- c(theta0, parms.curved$gwdegree.decay)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-1))
    eta0[wdegf] <- theta0["gwdegree"]*
      (iseq*exp(theta0["gwdegree.decay"]) +
       exp(2*theta0["gwdegree.decay"])*
       (((1-exp(-theta0["gwdegree.decay"]))^iseq)-1))
  }
  if(any(wodegf)){
    loc.gwodegree <- match("gwodegree",names(theta0))
    for(i in 1:length(m$terms)){
      if(m$terms[[i]]$name=="gwodegree"){
        parms.curved$gwodegree.decay <-  m$terms[[i]]$inputs[4] 
      }
    }
    names(parms.curved$gwodegree.decay) <- "gwodegree.decay"
    theta0 <- c(theta0, parms.curved$gwodegree.decay)
    theta1 <- c(theta1, FALSE)
    iseq <- (1:(network.size(nw)-1))
    eta0[wodegf] <- theta0["gwodegree"]*
      (iseq*exp(theta0["gwodegree.decay"]) +
       exp(2*theta0["gwodegree.decay"])*
       (((1-exp(-theta0["gwodegree.decay"]))^iseq)-1))
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
  if(any(geodf) | any(geosdf) | any(wdegf) | any(wodegf) | any(gwespf) ){
    nongeonames <- is.na(match(names(theta0),
                               c("gd.alpha", "gsd.alpha", "geodegree",
                                 "geospartner", "gwdegree.decay",
                                 "gwodegree.decay",
                                 "gwesp.alpha","gwdegree", "gwodegree",
                                 "gwesp")))
    eta0[!geodf & !geosdf & !wdegf & !wodegf & !gwespf] <- theta0[nongeonames]
    parms.curved$eta0 <-  eta0
  }else{
    parms.curved$eta0 <-  theta0
  }
  list(parms.curved=parms.curved,eta0=eta0,theta0=theta0,theta1=theta1,
       geodf=geodf,geosdf=geosdf,wdegf=wdegf,wodegf=wodegf,gwespf=gwespf)
}
