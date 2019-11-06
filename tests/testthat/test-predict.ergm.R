context("Test predict.formula() and predict.ergm()")

library(ergm)
library(magrittr)

test_that("works for edges model on small digraph", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    r.f <- predict(net ~ edges, log(1/5)) # predict.formula()
  )
  fit <- ergm(net ~ edges)
  expect_silent(
    r.e <- predict(fit) # predict.ergm()
  )
  expect_true( all.equal(unique(r.f$p), 1/6) )
  expect_identical(
    names(r.f),
    c("tail", "head", "p")
  )
  expect_identical(
    names(r.e),
    c("tail", "head", "p")
  )
  expect_true( all.equal(unique(r.e$p), 1/6) )
})



test_that("works for edges model on small graph", {
  net <- network.initialize(3, directed=FALSE)
  net[1,2] <- 1
  expect_silent(
    r.f <- predict(net ~ edges, log(1/2)) # predict.formula()
  )
  fit <- ergm(net ~ edges)
  expect_silent(
    r.e <- predict(fit) # predict.ergm()
  )
  expect_identical(
    names(r.f),
    c("tail", "head", "p")
  )
  expect_identical(
    names(r.e),
    c("tail", "head", "p")
  )
  expect_true( all.equal(unique(r.f$p), 1/3) )
  expect_true( all.equal(unique(r.f$p), 1/3) )
})





















# Scrapbook ---------------------------------------------------------------

# Some ideas for tests not [yet] implemented.

if(FALSE) {
  
  logit <- function(p, ...) log(p / (1-p), ...)
  expit <- function(lgit, ...) 1 / (1 + exp(-lgit))

  r1 <- predict_ergm(fit, flomarriage)
  
  data("flo")
  fit <- ergm(flomarriage ~ edges + gwesp(0.25, fixed = TRUE))
  # term: indicies
  z <- ergmMPLE( 
    update(fit$formula, . ~ . + indices), 
    output="matrix"
  )
  str(z)
  coef(fit)
  v <- ergm.eta(fit$coef, fit$etamap)
  z$predictor[,1:2] %*% v -> m
  
  p <- predict_formula(flomarriage ~ edges + gwesp(0.25, fixed = TRUE), c(-1, 0.5))
  
  
  # OTPs --------------------------------------------------------------------
  
  net <-
    igraph::make_graph(~A --+ C:D --+ B +-- A) %>%
    intergraph::asNetwork() 
  plot(net, displaylabels=TRUE)
  
  # Formula
  fit <- ergm(
    net ~ edges + dgwesp(decay=0.5, fixed=TRUE, type="OTP"),
    control = control.ergm(MCMLE.maxit = 2)
  )
  
  predict(
    net ~ edges + dgwesp(decay=0.5, fixed=TRUE, type="OTP"),
    net,
    c(0.4, 1)
  )
  
  predict(fit, net)
  ff <- fix.curved(fit)
  coef(fit)
  summary(ff$formula)
  
  # Sampson data ------------------------------------------------------------
  
  # Curved ERGM
  
  data(sampson)
  gest <- ergm(
    samplike ~ edges + gwesp(decay=0.5, fixed=FALSE),
    control=control.ergm(MCMLE.maxit=1)
  )
  coef(gest)
  summary(samplike ~ edges + gwesp(decay=0.5, fixed=FALSE))
  enf <- enformulate.curved(gest)
  ff <- fix.curved(gest)
  summary(enf$formula)
  summary(ff$formula)
  
  predict(gest, samplike)
  
  # IBE121 - play --------------------------------------------------------
  
  
  netw <- isnar::IBE121 %>%
    igraph::delete_edges(igraph::E(.)[question != "play"]) %>%
    intergraph::asNetwork()
  
  fit <- ergm(
    netw ~ edges + nodematch("female") + gwesp(decay=0.2, fixed=FALSE),
    control = control.ergm(MCMLE.maxit=1)
  )
  
  ff <- enformulate.curved(fit)
  
  
  summary(netw ~ edges + nodematch("female") + gwesp(decay=0.2, fixed=TRUE) )
  p <- cond_probs(fit, netw)  
  

  # Subsets -----------------------------------------------------------------
  
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  r <- predict(net ~ edges, -1.609, constraints= ~ fixallbut(net))
  theta <- -1.609
  predmat <- ergmMPLE(net ~ edges + indices, constraints = ~ fixallbut(net))$predictor
}
