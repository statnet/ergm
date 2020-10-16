
options(ergm.eval.loglik=FALSE)

data(florentine)

test_that("Stochastic Approximation produces similar results to MCMLE (linear ERGM)",{
  set.seed(2)

  mod.sa <- ergm(flomarriage~edges+triangle,control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+triangle)

  expect_equivalent(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1)
})

test_that("Stochastic Approximation produces similar results to MCMLE (curved ERGM)",{
  set.seed(2)

  mod.sa <- ergm(flomarriage~edges+gwdegree(),control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+gwdegree())

  expect_equivalent(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1)
})


test_that("Stochastic Approximation produces similar results to MCMLE (valued ERGM)",{
  set.seed(2)

  flomarriage %e% "w" <- 1

  mod.sa <- ergm(flomarriage~edges+transitiveweights, response="w", reference=~DiscUnif(0,1), control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+transitiveweights, response="w", reference=~DiscUnif(0,1))

  expect_equivalent(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1)
})