"ergm.curved.update" <- function(ecurved,theta0,m.expanded,g){
#
        geodf <- ecurved$geodf
        geosdf <- ecurved$geosdf
        wdegf <- ecurved$wdegf
        wodegf <- ecurved$wodegf
        gwespf <- ecurved$gwespf
        parms.curved <- ecurved$parms.curved
#
        if(any(geodf) | any(geosdf) | any(wdegf) | any(wodegf) | any(gwespf) ){
          eta0 <- rep(0,length=length(m.expanded$coef.names))
          names(eta0) <- m.expanded$coef.names
        }else{
          eta0 <- theta0
        }
        if(any(geodf)){
         iseq <- (1:(network.size(g)-1))
         parms.curved$geodegree.alpha <- theta0["gd.alpha"]
         eta0[geodf] <- theta0["geodegree"]*(exp(-theta0["gd.alpha"]*iseq)-1)
        }
        if(any(geosdf)){
         parms.curved$geospartner.alpha <- theta0["gsd.alpha"]
         iseq <- (1:(network.size(g)-2))
         eta0[geosdf] <- theta0["geospartner"]*(exp(-theta0["gsd.alpha"]*iseq)-1)
        }
        if(any(wdegf)){
         parms.curved$gwdegree.alpha <- theta0["gwdegree.decay"]
         iseq <- (1:(network.size(g)-1))
         eta0[wdegf] <- theta0["gwdegree"]*(iseq*exp(theta0["gwdegree.decay"]) +
          exp(2*theta0["gwdegree.decay"])*(((1-exp(-theta0["gwdegree.decay"]))^iseq)-1))
        }
        if(any(wodegf)){
         parms.curved$gwodegree.decay <- theta0["gwodegree.decay"]
         iseq <- (1:(network.size(g)-1))
         eta0[wodegf] <- theta0["gwodegree"]*(iseq*exp(theta0["gwodegree.decay"]) +
          exp(2*theta0["gwodegree.decay"])*(((1-exp(-theta0["gwodegree.decay"]))^iseq)-1))
        }
        if(any(gwespf)){
         parms.curved$gwesp.alpha <- theta0["gwesp.alpha"]
         iseq <- (1:(network.size(g)-2))
         eta0[gwespf] <- theta0["gwesp"]*
          exp(theta0["gwesp.alpha"])*(1-(1-exp(-theta0["gwesp.alpha"]))^iseq)
        }
        if(any(geodf) | any(geosdf) | any(wdegf) | any(wodegf) | any(gwespf) ){
         nongeonames <- is.na(match(names(theta0),
              c("gd.alpha", "gsd.alpha", "geodegree","geospartner",
                "gwdegree.decay","gwodegree.decay","gwesp.alpha",
                "gwdegree","gwodegree", "gwesp")))
         eta0[!geodf & !geosdf & !wdegf & !wodegf & !gwespf] <- theta0[nongeonames]
         parms.curved$eta0 <-  eta0
        }else{
         parms.curved$eta0 <-  theta0
        }
        list(parms.curved=parms.curved,eta0=eta0,theta0=theta0)
}
