#  File tests/testthat/test-term-For.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
test_that("For() operator (binary and valued) with list input", {
  data(florentine)
  add_eattr(flomarriage, "w")
  expect_equal(summary(flomarriage ~ For(~absdiff("wealth", pow=.), . = 1:3)),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3)))
  expect_equal(summary(flomarriage ~ For(~absdiff("wealth", pow=.), . = 1:3), response = "w"),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3), response = "w"))
})

test_that("For() operator (binary and valued) with formula input", {
  data(sampson)
  add_eattr(samplike, "w")
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = ~sort(unique(.%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = ~sort(unique(.%v%"group"))), response = "w"),
               summary(samplike ~ nodematch("group", diff=TRUE), response = "w"))
})

test_that("For() operator (binary and valued) with function input", {
  data(sampson)
  add_eattr(samplike, "w")
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = function(nw) sort(unique(nw%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = function(nw) sort(unique(nw%v%"group"))), response = "w"),
               summary(samplike ~ nodematch("group", diff=TRUE), response = "w"))
})

test_that("For() operator (binary and valued) with multiple list inputs", {
  data(florentine)
  add_eattr(flomarriage, "w")
  expect_equal(summary(flomarriage ~ For(~absdiff(a, pow=.), a = c("wealth", "priorates"), . = 1:3)),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3) +
                         absdiff("priorates", pow=1) + absdiff("priorates", pow=2) + absdiff("priorates", pow=3)))

  expect_equal(summary(flomarriage ~ For(. = 1:3, a = c("wealth", "priorates"), ~absdiff(a, pow=.))),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("priorates", pow=1) +
                         absdiff("wealth", pow=2) + absdiff("priorates", pow=2) +
                         absdiff("wealth", pow=3) + absdiff("priorates", pow=3)))

  expect_equal(summary(flomarriage ~ For(~absdiff(a, pow=.), a = c("wealth", "priorates"), . = 1:3), response = "w"),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3) +
                         absdiff("priorates", pow=1) + absdiff("priorates", pow=2) + absdiff("priorates", pow=3), response = "w"))

  expect_equal(summary(flomarriage ~ For(. = 1:3, a = c("wealth", "priorates"), ~absdiff(a, pow=.)), response = "w"),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("priorates", pow=1) +
                         absdiff("wealth", pow=2) + absdiff("priorates", pow=2) +
                         absdiff("wealth", pow=3) + absdiff("priorates", pow=3), response = "w"))
})
