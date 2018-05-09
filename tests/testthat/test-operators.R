context("test-operators.R")

data(florentine)

test_that("Simulation for passthrough() and .submodel() and .summary()", {
  text <- capture.output(out <- simulate(flomarriage~edges+degree(0)+absdiff("wealth")+passthrough(~edges+degree(0)+absdiff("wealth"))+submodel.test(~edges+degree(0)+absdiff("wealth"))+summary.test(~edges+degree(0)+absdiff("wealth")), statsonly=TRUE, nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1), coef=numeric(10)))
  text.out <- matrix(scan(textConnection(paste(text, collapse="")),quiet=TRUE),byrow=TRUE,ncol=3)
  text.out <- text.out[nrow(text.out)-nrow(out)+seq_len(nrow(out)),]
  
  expect_equivalent(out[,1:3],out[,4:6])
  expect_equivalent(out[,1:3],out[,7:9])
  expect_equivalent(out[,1:3],text.out)
})

data(sampson)
g <- samplike%v%"group"
sameg <- outer(g,g,"==")

test_that("Simulation for NodematchFilter() and F()", {
  out <- simulate(samplike~nodematch("group")+odegree(0:5, by="group", homophily=TRUE)+idegree(0:5, by="group", homophily=TRUE)+localtriangle(sameg)+
                    NodematchFilter(~edges+odegree(0:5)+idegree(0:5)+triangle,"group")+
                    F(~edges+odegree(0:5)+idegree(0:5)+triangle,~nodematch("group")), statsonly=TRUE, nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1), coef=numeric(42))

  expect_equivalent(out[,1:14], out[,15:28])
  expect_equivalent(out[,1:14],out[,29:42])
})

library(ergm.count)
data(zach)
test_that("Summary for the B() operator with nonzero criteria",{
  summ <- summary(zach~B(~edges+triangles, "nonzero") + B(~edges+triangles, ~nonzero), response="contexts")
  expect_equivalent(summ, rep(summary(zach~edges+triangles),2))
})

test_that("Summary for the B() operator with interval criteria",{
  summ <- summary(zach~B(~edges+triangles, ~ininterval(3,5,c(FALSE,FALSE))), response="contexts")
  expect_equivalent(summ, summary(network_view(zach, "contexts", function(x) x>=3 & x<=5)~edges+triangles))
})
