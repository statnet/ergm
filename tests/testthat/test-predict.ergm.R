context("Test predict.formula() and predict.ergm()")

library(ergm)

test_that("predict.formula(type=) give correct results", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(net ~ edges, theta = log(1/5), type="response") # predict.formula()
  )
  expect_silent(
    p.link <- predict(net ~ edges, theta = log(1/5), type="link") # predict.formula()
  )
  expect_true(
    all.equal(p.link$p, log(p.prob$p / (1 - p.prob$p)))
  )
})


test_that("predict.formula(conditional=FALSE)", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(
      net ~ edges, 
      theta = log(1/5),
      nsim = 5,
      type="response", 
      conditional=FALSE
    )
  )
  
})


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






test_that("predict.formula(output='matrix') works correctly", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(net ~ edges, theta = log(1/5), type="response", output="matrix")
  )
  
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
  summary(fit)
  
  simulate(fit, nsim=1)

  predict(fit)
  predict(fit, conditional=FALSE, nsim=10)
  
  predict(
    net ~ edges + dgwesp(decay=0.5, fixed=TRUE, type="OTP"),
    theta = c(-191, 66)
  )

  simulate_formula(
    net ~ edges + dgwesp(decay=0.5, fixed=TRUE, type="OTP"),
    coef = c(-191, 66)
  )

  # Sampson data ------------------------------------------------------------
  
  # Curved ERGM
  
  data(sampson)
  gest <- ergm(
    samplike ~ edges + gwesp(decay=0.5, fixed=TRUE)
  )
  coef(gest)
  summary(samplike ~ edges + gwesp(decay=0.5, fixed=TRUE))

  predict(gest)
  predict(gest, conditional=FALSE, output="matrix")
  
  # predict on net different than the one the model was fit to --------------
  
  predict(
    net ~ edges + gwesp(decay=0.5, fixed=TRUE),
    theta = coef(gest)
  )  
  
  
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
  
  
  
  data(sampson)
  gest2 <- ergm(
    samplike ~ edges + nodematch("group")
  )
  summary(gest2)
  g <- gof(gest2)
  plot(g)
  
  # Probs for all dyads
  predict(gest2)
  
  # Probs for dyads in `net`

  
}

