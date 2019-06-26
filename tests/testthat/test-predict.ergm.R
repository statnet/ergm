context("Test predict.formula()")

library(ergm)
library(magrittr)

test_that("works on a small digraph",{
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    r <- predict(net ~ edges, -1.609)
  )
})




context("Test predict.ergm()")

test_that("works on a small digraph",{
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  fit <- ergm(net ~ edges)
  expect_silent(
    r <- predict(fit)
  )
})


# TODO --------------------------------------------------------------------

if(FALSE) {
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
}
