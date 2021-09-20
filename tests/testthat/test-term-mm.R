#  File tests/testthat/test-term-mm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# a bipartite nw
set.seed(143)
b1 <- floor(runif(60, 1,100))
b2 <- floor(runif(60, 101, 130))
exbip.el <- cbind(b1,b2)
bipnw <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw %v% "Letter" <- letters[1:3]
bipnw %v% "Cost" <- c(3,2,1)
                          

# another bipartite nw with more ties and 2 attributes
set.seed(258)
b1 <- floor(runif(150, 1,200))
b2 <- floor(runif(150, 201, 400))
exbip.el <- cbind(b1,b2)
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color


# a directed nw
data(sampson)
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")


# an undirected nw
data(faux.mesa.high)
fmh <- faux.mesa.high
set.seed(7)
set.edge.attribute(fmh, "GradeMet", rbinom(203, 6, .5))


# a small undirected nw w/ lots o' triangles
set.seed(20)
t<-trunc(runif(160, 1, 20))
set.seed(21)
h<-trunc(runif(160, 1, 20))
el <- cbind(t,h)
bad <- which(el[,2]==el[,1])
el[bad,2] = el[bad,2]+1
unnw <- network(el, directed=FALSE)
unnw %v% "Pet" <- c("dog", "cat")

test_that("Undirected mm() summary", {
  s.a <- summary(fmh ~ mm("Grade"))
  expect_equal(s.a, c(`mm[Grade=7,Grade=8]` = 0, `mm[Grade=8,Grade=8]` = 33,
                      `mm[Grade=7,Grade=9]` = 0, `mm[Grade=8,Grade=9]` = 2,
                      `mm[Grade=9,Grade=9]` = 23, `mm[Grade=7,Grade=10]` = 1,
                      `mm[Grade=8,Grade=10]` = 4, `mm[Grade=9,Grade=10]` = 7,
                      `mm[Grade=10,Grade=10]` = 9, `mm[Grade=7,Grade=11]` = 1,
                      `mm[Grade=8,Grade=11]` = 2, `mm[Grade=9,Grade=11]` = 6,
                      `mm[Grade=10,Grade=11]` = 1, `mm[Grade=11,Grade=11]` = 17,
                      `mm[Grade=7,Grade=12]` = 1, `mm[Grade=8,Grade=12]` = 1,
                      `mm[Grade=9,Grade=12]` = 4, `mm[Grade=10,Grade=12]` = 5,
                      `mm[Grade=11,Grade=12]` = 5, `mm[Grade=12,Grade=12]` = 6))
})


test_that("Undirected mm() summary with expression, levels, and levels2 selectors", {
  s.a <- summary(fmh ~ mm(~Grade >= 10, levels=-1, levels2=NULL))
  expect_equal(s.a, c(`mm[Grade>=10=TRUE,Grade>=10=TRUE]` = 43))
})

test_that("Undirected mm() summary with expression, levels selector, no terms", {
  s.a <- summary(fmh ~ mm(~Grade >= 10, levels=-1))
  expect_equal(s.a, numeric(0))
})


test_that("Directed mm() ERGM", {
  e.a <- ergm(samplike ~ mm("group", levels2=TRUE))
  expect_equivalent(coef(e.a),
                    c(0.191055236762708, -3.29583686600359,
                      -2.17475172140943, -2.5649493574587,
                      1.6094379124341, -3.29583686600359,
                      -1.49165487677766, -1.09861228866811,
                      0.916290731874155))
})

test_that("Bipartite mm() summary", {
  s.ab <- summary(bipnw ~ mm("Letter"))
  expect_equal(s.ab, c(`mm[Letter=b,Letter=a]` = 8, `mm[Letter=c,Letter=a]` = 8,
                       `mm[Letter=a,Letter=b]` = 7, `mm[Letter=b,Letter=b]` = 7,
                       `mm[Letter=c,Letter=b]` = 5, `mm[Letter=a,Letter=c]` = 4,
                       `mm[Letter=b,Letter=c]` = 6, `mm[Letter=c,Letter=c]` = 6))
})


test_that("Bipartite mm() ERGM with level2 filter", {
  e.ab <- ergm(bipnw ~ mm("Letter", levels2=-(2:6)))
  expect_equivalent(coef(e.ab),
                    c(-3.49650756126405, -4.43081679479004,
                      -3.98898404450688, -3.98898404450688 ))
})

fmh.mm.Race <-
  c(`mm[Race=Black,Race=Black]` = 0, `mm[Race=Black,Race=Hisp]` = 8,
    `mm[Race=Hisp,Race=Hisp]` = 53, `mm[Race=Black,Race=NatAm]` = 13,
    `mm[Race=Hisp,Race=NatAm]` = 41, `mm[Race=NatAm,Race=NatAm]` = 46,
    `mm[Race=Black,Race=Other]` = 0, `mm[Race=Hisp,Race=Other]` = 1,
    `mm[Race=NatAm,Race=Other]` = 0, `mm[Race=Other,Race=Other]` = 0,
    `mm[Race=Black,Race=White]` = 5, `mm[Race=Hisp,Race=White]` = 22,
    `mm[Race=NatAm,Race=White]` = 10, `mm[Race=Other,Race=White]` = 0,
    `mm[Race=White,Race=White]` = 4)


test_that("Undirected mm() summary with level2 filter", {
  s.ab2 <- summary(fmh ~ mm("Race", levels2=-1))
  expect_equivalent(s.ab2, fmh.mm.Race[-1])
})

test_that("Undirected mm() summary with level2 filter by logical matrix", {
  M <- matrix(TRUE, 5, 5)
  M[1,1] <- M[3,2] <- M[2,3] <- FALSE
  s.ab2 <- summary(fmh ~ mm("Race", levels2=M))
  expect_equivalent(s.ab2, fmh.mm.Race[M[upper.tri(M, TRUE)]])
})

test_that("Undirected mm() summary with level2 filter by numeric matrix", {
  M <- cbind(2,3)
  s.ab2 <- summary(fmh ~ mm("Race", levels2=M))
  expect_equivalent(s.ab2, fmh.mm.Race[5])
})

test_that("Directed mm() ERGM with level2 filter", {
  e.ab2 <- ergm(samplike ~ mm("Trinity", levels2=-(3:9)))
  expect_equivalent(coef(e.ab2),
                    c(-1.01160091056776, -0.693147180549145))
})

test_that("Directed mm() ERGM with level2 filter", {
  e.ab2 <- ergm(samplike ~ mm("Trinity", levels2=-(3:9)))
  expect_equivalent(coef(e.ab2),
                    c(-1.01160091056776, -0.693147180549145))
})

samp.mm.Trinity <-
  c(`mm[Trinity=F,Trinity=F]` = 8, `mm[Trinity=H,Trinity=F]` = 12,
    `mm[Trinity=S,Trinity=F]` = 12, `mm[Trinity=F,Trinity=H]` = 8,
    `mm[Trinity=H,Trinity=H]` = 5, `mm[Trinity=S,Trinity=H]` = 9,
    `mm[Trinity=F,Trinity=S]` = 13, `mm[Trinity=H,Trinity=S]` = 12,
    `mm[Trinity=S,Trinity=S]` = 9)


test_that("Directed mm() summary with level2 matrix filter", {
  M <- matrix(FALSE, 3, 3)
  M[1,2] <- M[1,3] <- TRUE
  s.ab2 <- summary(samplike ~ mm("Trinity", levels2=M))
  expect_equivalent(s.ab2, samp.mm.Trinity[c(M)])
})

test_that("Undirected mm() marginal summary", {
  s.a <- summary(fmh ~ mm(.~Grade))
  expect_equal(s.a, c(`mm[.,Grade=8]` = 75, `mm[.,Grade=9]` = 65,
                      `mm[.,Grade=10]` = 36, `mm[.,Grade=11]` = 49,
                      `mm[.,Grade=12]` = 28))
})

test_that("Undirected mm() asymmetric two-sided summary", {
  s.a.tab <- as.matrix(fmh, matrix.type="edgelist") %>%
    rbind(.,.[,2:1,drop=FALSE]) %>%
    (function(x) data.frame(Race=factor(fmh%v%"Race")[x[,1]], Grade=factor(fmh%v%"Grade")[x[,2]])) %>%
    table()
  s.a <- summary(fmh ~ mm(Race~Grade, levels2=TRUE))
  expect_equivalent(s.a, c(s.a.tab))
})

test_that("Undirected mm() asymmetric two-sided summary and pipe operator on one side", {
  Grade <- fmh%v%"Grade" %>% factor()
  Race <- fmh%v%"Race" %>% replace(., . %in% c("Black","White","Other"), "BWO") %>% factor()

  s.a.tab <- as.matrix(fmh, matrix.type="edgelist") %>%
    rbind(.,.[,2:1,drop=FALSE]) %>%
    (function(x) data.frame(Grade=Grade[x[,2]], Race=Race[x[,1]])) %>%
    table()
  s.a <- summary(fmh ~ mm(Grade~(~Race) %>% COLLAPSE_SMALLEST(3,"BWO"), levels2=TRUE))
  expect_equivalent(s.a, c(s.a.tab))
})

test_that("Undirected mm() marginal summary with fixed levels set", {
  s.a <- summary(fmh ~ mm(.~Grade, levels=I(c(7,6,9,8)), levels2=TRUE))
  expect_equivalent(s.a, c(153, 0, 65, 75))
})

test_that("Undirected mm() summary with fixed levels set", {
  s.a <- summary(fmh ~ mm("Grade", levels=I(c(7,6,9,8)), levels2=TRUE))
  expect_equal(s.a, c(`mm[Grade=7,Grade=7]` = 75, `mm[Grade=7,Grade=6]` = 0,
                      `mm[Grade=6,Grade=6]` = 0, `mm[Grade=7,Grade=9]` = 0,
                      `mm[Grade=6,Grade=9]` = 0, `mm[Grade=9,Grade=9]` = 23,
                      `mm[Grade=7,Grade=8]` = 0, `mm[Grade=6,Grade=8]` = 0,
                      `mm[Grade=9,Grade=8]` = 2, `mm[Grade=8,Grade=8]` = 33))
})


test_that("Undirected valued mm() sum summary", {
  s.a <- summary(fmh ~ mm("Grade"), response="GradeMet")
  expect_equal(s.a,
               c(`mm.sum[Grade=7,Grade=8]` = 0, `mm.sum[Grade=8,Grade=8]` = 96,
                 `mm.sum[Grade=7,Grade=9]` = 0, `mm.sum[Grade=8,Grade=9]` = 6, 
                 `mm.sum[Grade=9,Grade=9]` = 75, `mm.sum[Grade=7,Grade=10]` = 3, 
                 `mm.sum[Grade=8,Grade=10]` = 15, `mm.sum[Grade=9,Grade=10]` = 21, 
                 `mm.sum[Grade=10,Grade=10]` = 24, `mm.sum[Grade=7,Grade=11]` = 3, 
                 `mm.sum[Grade=8,Grade=11]` = 6, `mm.sum[Grade=9,Grade=11]` = 19, 
                 `mm.sum[Grade=10,Grade=11]` = 4, `mm.sum[Grade=11,Grade=11]` = 54, 
                 `mm.sum[Grade=7,Grade=12]` = 2, `mm.sum[Grade=8,Grade=12]` = 5, 
                 `mm.sum[Grade=9,Grade=12]` = 11, `mm.sum[Grade=10,Grade=12]` = 17, 
                 `mm.sum[Grade=11,Grade=12]` = 18, `mm.sum[Grade=12,Grade=12]` = 16
                 ))
})

test_that("Undirected valued mm() nonzero summary", {
  fmh %ergmlhs% "response" <- "GradeMet"
  s.a <- summary(fmh ~ mm("Grade", form="nonzero"))

  expect_equal(s.a, c(`mm.nonzero[Grade=7,Grade=8]` = 0,
                      `mm.nonzero[Grade=8,Grade=8]` = 33,
                      `mm.nonzero[Grade=7,Grade=9]` = 0,
                      `mm.nonzero[Grade=8,Grade=9]` = 2,
                      `mm.nonzero[Grade=9,Grade=9]` = 22,
                      `mm.nonzero[Grade=7,Grade=10]` = 1,
                      `mm.nonzero[Grade=8,Grade=10]` = 4,
                      `mm.nonzero[Grade=9,Grade=10]` = 7,
                      `mm.nonzero[Grade=10,Grade=10]` = 8,
                      `mm.nonzero[Grade=7,Grade=11]` = 1,
                      `mm.nonzero[Grade=8,Grade=11]` = 2,
                      `mm.nonzero[Grade=9,Grade=11]` = 6,
                      `mm.nonzero[Grade=10,Grade=11]` = 1,
                      `mm.nonzero[Grade=11,Grade=11]` = 17,
                      `mm.nonzero[Grade=7,Grade=12]` = 1,
                      `mm.nonzero[Grade=8,Grade=12]` = 1,
                      `mm.nonzero[Grade=9,Grade=12]` = 4,
                      `mm.nonzero[Grade=10,Grade=12]` = 5,
                      `mm.nonzero[Grade=11,Grade=12]` = 5,
                      `mm.nonzero[Grade=12,Grade=12]` = 6))
})


test_that("Undirected mm() summary with just one unique attribute level and levels imposed exogenously", {
  nw2 <- network.initialize(2, directed = FALSE)
  nw2[1, 2] <- TRUE
  nw2 %v% "x" <- c("b", "b")

  s.a <- summary(nw2 ~ mm("x", levels = letters[1:3], levels2 = TRUE))
  expect_equal(s.a, c(`mm[x=a,x=a]` = 0, `mm[x=a,x=b]` = 0, `mm[x=b,x=b]` = 1,
                      `mm[x=a,x=c]` = 0, `mm[x=b,x=c]` = 0, `mm[x=c,x=c]` = 0))
})
