#  File ergm/R/ergm.pen.glm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.pen.glm <- function(formula = attr(data, "formula"), 
  data = sys.parent(), alpha = 0.05, 
  maxit = 25, maxhs = 5, epsilon = 0.0001, maxstep = 10, 
  start=NULL,
  weights=NULL)
{
  y <- as.vector(model.extract(model.frame(formula, data = data), "response"))
  n <- length(y)
  x <- model.matrix(formula, data = data)  # Model-Matrix

#  aaa <- rep(seq(along=weights),weights)
#  x <- as.matrix(x[aaa,])
#  y <- y[aaa]
#  n <- length(y)
#  weights <- y-y+1

  k <- ncol(x)  # Anzahl Effekte

  int <- 0
  coltotest <-1:k

  if(missing(weights)){weights <- rep(1,length=n)}
  beta <- c(log((sum(y*weights)/sum((1-y)*weights))),
            rep(0, k - 1))
  if(!missing(start) && !is.null(start) && ncol(x)==length(start)){
    beta[1] <- beta[1] - sum((x %*% start)*weights)
  }
  iter <- 0
  pi <- as.vector(1/(1 + exp( - x %*% beta)))
  loglik <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
  XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")  #### X' (W ^ 1/2)
  Fisher <- crossprod(XW2)  #### X' W  X
  loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
  repeat {
   iter <- iter + 1
   XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
   Fisher <- crossprod(XW2)  #### X' W  X
   covs <- robust.inverse(Fisher)  ### (X' W  X) ^ -1
#  H <- crossprod(XW2, covs) %*% XW2
#  H <- XW2 %*% covs %*% t(XW2)
   diagH <- pi
   for(i in seq(along=diagH)){
    diagH[i] <- XW2[i,] %*% covs %*% XW2[i,]
   }
#  U.star <- crossprod(x, y - pi)
#  U.star <- crossprod(x, (y - pi)*weights)
#  U.star <- crossprod(x, (y - pi + diagH * (0.5 - pi))*weights)
   U.star <- crossprod(x, (y - pi)*weights + diagH * (0.5 - pi))
   delta <- as.vector(covs %*% U.star)
   mx <- max(abs(delta))/maxstep
   if(mx > 1){
    delta <- delta/mx
   }
   beta <- beta + delta
   loglik.old <- loglik
   for(halfs in 1:maxhs) {
## Half-Steps
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
    loglik <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
    XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
    Fisher <- crossprod(XW2) #### X' W  X
    loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    if(loglik > loglik.old) break
    beta <- beta - delta * 2^( - halfs) ##beta-A enderung verkleinern
   }
   if(iter == maxit | sum(abs(delta)) <= epsilon){break}
  }

  fit <- list(coefficients = beta, alpha = alpha, var = covs, 
              df = (k-int), 
              loglik = loglik,
#             loglik = logistftest(formula=formula, data=data, 
#                                  test=coltotest,
#                                  weights=weights),
              iter = iter, n = n, 
              terms = colnames(x), y = y, 
              formula = as.formula(formula), call=match.call(),
              data = data, 
              model.matrix = x, 
              method = "pen.glm.fit",
              linear.predictors=as.vector(x %*% beta))
  class(fit) <- c("pen.glm")
# vars <- diag(covs)
  fit
}
logistftest <- function(formula = attr(data, "formula"),
  data = sys.parent(), test, values, maxit = 25, maxhs = 5, 
  epsilon = 0.0001, maxstep = 10, start, weights=NULL)
{
  #n <- nrow(data)
  y <- model.extract(model.frame(formula, data = data), "response")
  n <- length(y)
  x <- model.matrix(formula, data = data)  ## Model-Matrix
  cov.name <- labels(x)[[2]]
  k <- ncol(x)
  int <- 0
  coltotest <-1:k

  if(missing(weights)){weights <- rep(1,length=n)}
  beta <- c(log((sum(y*weights)/sum((1-y)*weights))),
            rep(0, k - 1))  
  if(!missing(start) && !is.null(start) && ncol(x)==length(start)){
    beta[1] <- beta[1] - sum((x %*% start)*weights)
  }
##berechne Startwerte
  iter <- 0
  loglik <- rep(0, 2)
  pi <- as.vector(1/(1 + exp( - x %*% beta)))
  if(missing(start)) {
################## coxphfplot braucht dies nicht! ###
  loglik[2] <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
  XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
  Fisher <- crossprod(XW2)  #### X' W  X
# loglik[2] <- loglik[2] + 0.5 * determinant(Fisher)$modulus[1]
  repeat {
    iter <- iter + 1
    XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
    Fisher <- crossprod(XW2)  #### X' W  X
    covs <- robust.inverse(Fisher)  ### (X' W  X) ^ -1
    diagH <- pi
    for(i in seq(along=diagH)){
       diagH[i] <- XW2[i,] %*% covs %*% XW2[i,]
    }
#   U.star <- crossprod(x, y - pi + diagH * (0.5 - pi))
    U.star <- crossprod(x, y - pi + diagH * (0.5 - pi))
    U.star <- crossprod(x, (y - pi)*weights)
    delta <- as.vector(covs %*% U.star)
    mx <- max(abs(delta))/maxstep
    if(mx > 1) delta <- delta/mx
    beta <- beta + delta
    loglik.old <- loglik[2]
    for(halfs in 1:maxhs) {
## 5 Half-Steps
      pi <- as.vector(1/(1 + exp( - x %*% beta)))
      loglik[2] <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
      XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")
      Fisher <- crossprod(XW2)  #### X' W  X
#     loglik[2] <- loglik[2] + 0.5 * determinant(Fisher)$modulus[1]
      if(loglik[2] > loglik.old) break
      beta <- beta - delta * 2^( - halfs)  
    }
    if(iter == maxit | sum(abs(delta)) <= epsilon) break
  }
  }
########################################################
## Labels der Test-Fakt.
  if(missing(test))
    test <- coltotest
  if(is.vector(test))
    cov.name2 <- cov.name[test]
  else cov.name2 <- labels(model.matrix(test, data = data))[[2]]
  pos <- match(cov.name2, cov.name)  ## Position der Testfakt.
  OK <- !is.na(pos)
  pos <- pos[OK]
  cov.name2 <- cov.name2[OK]
  k2 <- length(cov.name2)  ## Anzahl Faktoren
  if(!missing(start))
    offset <- start
  else offset <- rep(0, k)  ## Vektor der fixierten Werte
  if(!missing(values))
    offset[pos] <- values
  beta <- offset  ########################################
  iter <- 0
  pi <- as.vector(1/(1 + exp( - x %*% beta)))
  loglik[2] <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
  XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")  #### X' (W ^ 1/2)
  Fisher <- crossprod(XW2)  #### X' W  X
# loglik[1] <- loglik[1] + 0.5 * determinant(Fisher)$modulus[1]
  repeat {
   if(k2 == k) break  ## -> Overall Test
   iter <- iter + 1
   XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")  #### X' (W ^ 1/2)
   Fisher <- crossprod(XW2)  #### X' W  X
   covs <- robust.inverse(Fisher)  ### (X' W  X) ^ -1
   diagH <- pi
   for(i in seq(along=diagH)){
     diagH[i] <- XW2[i,] %*% covs %*% XW2[i,]
   }
#  U.star <- crossprod(x, y - pi + diagH * (0.5 - pi))
   U.star <- crossprod(x, (y - pi)*weights)
   XX.XW2 <- sweep(x[,  - pos, drop = FALSE], 1, (weights*(pi * (1 - pi)))^0.5, "*")
   #### Teil von X' (W ^ 1/2)
   XX.Fisher <- crossprod(XX.XW2)  #### Teil von  X' W  X
   XX.covs <- matrix(0, k, k)
   XX.covs[ - pos,  - pos] <- robust.inverse(XX.Fisher)  
   ### aufblasen der Cov-Matrix
   delta <- as.vector(XX.covs %*% U.star)
   mx <- max(abs(delta))/maxstep
   if(mx > 1)
     delta <- delta/mx
   beta <- beta + delta
   loglik.old <- loglik[1]
   for(halfs in 1:maxhs) {
## Half-Steps
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
    loglik[1] <- sum((y * log(pi) + (1 - y) * log(1 - pi))*weights)
    XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")
    Fisher <- crossprod(XW2)  #### X' W  X
#   loglik[1] <- loglik[1] + 0.5 * determinant(Fisher)$modulus[1]
    if(loglik[1] > loglik.old) break
    beta <- beta - delta * 2^( - halfs)  
   ##beta-Aenderung verkleinern
   }
   if(iter == maxit | sum(abs(delta)) <= epsilon) break
  }
#######################
  loglik
}
model.matrix.pen.glm <- function(object, ...)
{
	object$model.matrix
}
