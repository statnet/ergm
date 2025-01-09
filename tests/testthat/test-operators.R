#  File tests/testthat/test-operators.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(florentine)

test_that("Simulation for Passthrough() and .submodel() and .summary()", {
  text <- capture.output(
    out <- simulate(
      flomarriage ~ edges+degree(0)+absdiff("wealth")+
        Passthrough(~edges+degree(0)+absdiff("wealth"))+
        Passthrough(~edges+degree(0)+absdiff("wealth"), submodel=FALSE)+
        submodel.test(~edges+degree(0)+absdiff("wealth"))+
        summary.test(~edges+degree(0)+absdiff("wealth")),
      output="stats", nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1), coef=numeric(13)))
  text.out <- matrix(scan(textConnection(paste(text, collapse="")),quiet=TRUE),byrow=TRUE,ncol=3)
  text.out <- text.out[nrow(text.out)-nrow(out)+seq_len(nrow(out)),]
  
  expect_equal(out[,1:3],out[,4:6], ignore_attr=TRUE)
  expect_equal(out[,1:3],out[,7:9], ignore_attr=TRUE)
  expect_equal(out[,1:3],out[,10:12], ignore_attr=TRUE)
  expect_equal(out[,1:3],text.out, ignore_attr=TRUE)
})

data(sampson)
g <- samplike%v%"group"
sameg <- outer(g,g,"==")

test_that("Simulation for NodematchFilter() and F()", {
  out <- simulate(samplike~nodematch("group")+odegree(0:5, by="group", homophily=TRUE)+idegree(0:5, by="group", homophily=TRUE)+localtriangle(sameg)+
                    NodematchFilter(~edges+odegree(0:5)+idegree(0:5)+triangle,"group")+
                    F(~edges+odegree(0:5)+idegree(0:5)+triangle,~nodematch("group"))+
                    edges+
                    F(~edges, ~!nodematch("group")),
                  output="stats", nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1), coef=numeric(44))

  expect_equal(out[,1:14],out[,15:28], ignore_attr=TRUE)
  expect_equal(out[,1:14],out[,29:42], ignore_attr=TRUE)
  expect_equal(out[,1]+out[,44],out[,43], ignore_attr=TRUE)
})

test_that("Summary for F() with complex form", {
  m <- abs(outer(w <- flomarriage %v% "wealth", w, FUN="-"))[c(as.matrix(flomarriage))!=0]
  out <- summary(flomarriage ~
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") == 93) +
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") < 5) +
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") != 5) +
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") <= 5) +
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") > 5) +
                   F(~edges + absdiff("wealth"), ~absdiff("wealth") >= 5))

  expect_equal(out,
               sapply(list(m==93, m[m==93],
                           m<5, m[m<5],
                           m!=5, m[m!=5],
                           m<=5, m[m<=5],
                           m>5, m[m>5],
                           m>=5, m[m>=5]),
                      sum)/2, ignore_attr=TRUE)
})

test_that("Symmetrize() summary", {
  m <- as.matrix(samplike)
  expect_equal(
    c(sum(m*t(m))/2, sum(m+t(m)>0)/2, sum(m[lower.tri(m)]), sum(m[upper.tri(m)])),
    summary(samplike ~ Symmetrize(~edges,"strong") + Symmetrize(~edges,"weak") + Symmetrize(~edges,"lower") + Symmetrize(~edges,"upper")),
    ignore_attr=TRUE
  )
})

test_that("S() summary directed->bipartite", {
  m <- as.matrix(samplike)
  b1 <- sample.int(network.size(samplike), 5)
  b2 <- sample(setdiff(seq_len(network.size(samplike)), b1), 4)

  expect_equal(
    c(sum(m[b1,b2])),
    summary(samplike ~ S(~edges,I(b1)~I(b2))), ignore_attr=TRUE
  )
})

test_that("S() summary undirected->bipartite", {
  m <- as.matrix(flomarriage)
  b1 <- sample.int(network.size(flomarriage), 5)
  b2 <- sample(setdiff(seq_len(network.size(flomarriage)), b1), 4)

  expect_equal(
    c(sum(m[b1,b2])),
    summary(flomarriage ~ S(~edges,I(b1)~I(b2))), ignore_attr=TRUE
  )
})

test_that("S() summary directed->directed", {
  m <- as.matrix(samplike)
  i <- sample.int(network.size(samplike), 5)

  expect_equal(
    c(sum(m[i,i])),
    summary(samplike ~ S(~edges,~i)), ignore_attr=TRUE
  )
})


test_that("S() summary undirected->undirected", {
  m <- as.matrix(flomarriage)
  i <- sample.int(network.size(flomarriage), 5)

  expect_equal(
    c(sum(m[i,i])/2),
    summary(flomarriage ~ S(~edges,~i)), ignore_attr=TRUE
  )
})


test_that("Binary Label() summary", {
  expect_equal(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), "abc")),
    summary(flomarriage ~ edges+absdiff("wealth")), ignore_attr=TRUE
  )

  expect_named(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), "abc")),
    c("abc(edges)","abc(absdiff.wealth)")
  )

  expect_named(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), "abc", "prepend")),
    c("abcedges","abcabsdiff.wealth")
  )

  expect_named(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), c("abc","def"), "append")),
    c("edgesabc","absdiff.wealthdef")
  )

  expect_named(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), c("abc","def"), "replace")),
    c("abc","def")
  )

  expect_named(
    summary(flomarriage ~ Label(~edges+absdiff("wealth"), ~gsub(".","!",.,fixed=TRUE))),
    c("edges","absdiff!wealth")
  )
})

test_that("Binary Label() estimation, offsets, and curved terms", {
  expect_equal(
    coef(ergm(flomarriage ~ Label(~edges+offset(absdiff("wealth")), "abc"), offset.coef=-.5)),
    coef(ergm(flomarriage ~ edges+offset(absdiff("wealth")), offset.coef=-.5)), ignore_attr=TRUE
  )

  expect_equal(
    coef(ergm(flomarriage ~ Label(~edges+offset(gwesp), "abc"), offset.coef=c(-.5,1), estimate="MPLE")),
    coef(ergm(flomarriage ~ edges+offset(gwesp), offset.coef=c(-.5,1), estimate="MPLE")), ignore_attr=TRUE
  )

  ## list label
  ca <- c("abc", paste0("def", 1:14))
  cu <- c("abc", "ijk", "lmn")
  f <- flomarriage ~ Label(~edges+gwesp, list(cu, ca), "replace")
  expect_named(coef(ergm(f, estimate="MPLE")), cu)
  expect_named(summary(f), ca)

  ## list label with omitted vector
  f <- flomarriage ~ Label(~edges+gwesp, list(cu, NULL), "replace")
  f0 <- flomarriage ~edges+gwesp
  expect_named(coef(ergm(f, estimate="MPLE")), cu)
  expect_named(summary(f), c("abc", names(summary(f0))[-1]))

  f <- flomarriage ~ Label(~edges+gwesp, list(NA, ca), "replace")
  f0 <- flomarriage ~edges+gwesp
  expect_named(coef(ergm(f, estimate="MPLE")), param_names(ergm(f, estimate="MPLE")))
  expect_named(summary(f), ca)

  ## mapper label
  f <- flomarriage ~ Label(~edges+gwesp, ~gsub("[.#]","!",.))
  expect_named(coef(ergm(f, estimate="MPLE")),
               c("edges", "gwesp", "gwesp!decay"))
  expect_named(summary(f), c("edges", paste0("esp!", 1:14)))

})


library(ergm.count)
data(zach)
test_that("Summary for the B() operator with nonzero criteria",{
  summ <- summary(zach~B(~edges+triangles+degree(0:5), "nonzero") + B(~edges+triangles+degree(0:5), ~nonzero), response="contexts")
  expect_equal(summ, rep(summary(zach~edges+triangles+degree(0:5)),2), ignore_attr=TRUE)
})

test_that("Summary for the B() operator with interval criteria",{
  summ <- summary(zach~B(~edges+triangles+degree(0:5), ~ininterval(3,5,c(FALSE,FALSE))), response="contexts")
  expect_equal(summ, summary(zach~edges+triangles+degree(0:5), response= ~ contexts>=3 & contexts<=5), ignore_attr=TRUE)
})

test_that("Valued Label() summary", {
  expect_equal(
    summary(zach ~ Label(~edges+absdiff("faction.id"), "abc")),
    summary(zach ~ edges+absdiff("faction.id")), ignore_attr=TRUE
  )

  expect_named(
    summary(zach ~ Label(~edges+absdiff("faction.id"), "abc")),
    c("abc(edges)","abc(absdiff.faction.id)")
  )

  expect_named(
    summary(zach ~ Label(~edges+absdiff("faction.id"), "abc", "prepend")),
    c("abcedges","abcabsdiff.faction.id")
  )

  expect_named(
    summary(zach ~ Label(~edges+absdiff("faction.id"), c("abc","def"), "append")),
    c("edgesabc","absdiff.faction.iddef")
  )

  expect_named(
    summary(zach ~ Label(~edges+absdiff("faction.id"), c("abc","def"), "replace")),
    c("abc","def")
  )

  expect_named(
    summary(zach ~ Label(~edges+absdiff("faction.id"), ~gsub(".","!",.,fixed=TRUE))),
    c("edges","absdiff!faction!id")
  )
})

test_that("Interaction terms", {
  # TODO: Need better tests.
  expect_equal(summary(flomarriage~edges:absdiff("wealth") + absdiff("wealth"):edges), summary(flomarriage~absdiff("wealth") + absdiff("wealth")), ignore_attr = TRUE)
  expect_equal(summary(flomarriage~edges*absdiff("wealth") + absdiff("wealth")*edges), summary(flomarriage~edges + absdiff("wealth")+ absdiff("wealth") + absdiff("wealth") + edges + absdiff("wealth")), ignore_attr = TRUE)
})

test_that("Interaction terms handling of interact.dependent", {
  expect_error(summary(flomarriage~triangles:absdiff("wealth")), ".*poorly defined.*")
  expect_warning(summary(flomarriage~triangles:absdiff("wealth"), interact.dependent = "warning"), ".*poorly defined.*")
  expect_message(summary(flomarriage~triangles:absdiff("wealth"), interact.dependent = "message"), ".*poorly defined.*")
})
