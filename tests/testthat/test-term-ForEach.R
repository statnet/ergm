test_that("ForEach() operator with list input", {
  data(florentine)
  expect_equal(summary(flomarriage ~ ForEach(~absdiff("wealth", pow=.), ".", 1:3)),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3)))
})

test_that("ForEach() operator with formula input", {
  data(sampson)
  expect_equal(summary(samplike ~ ForEach(~nodematch("group", levels=., diff=TRUE), ".", ~sort(unique(.%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
})

test_that("ForEach() operator with function input", {
  data(sampson)
  expect_equal(summary(samplike ~ ForEach(~nodematch("group", levels=., diff=TRUE), ".", function(nw) sort(unique(nw%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
})
