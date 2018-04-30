context("test-multinet.R")

data(samplk)

logit <- function(p) log(p/(1-p))
times <- 1:3

samplk1%n%"t" <- times[1]
samplk2%n%"t" <- times[2]
samplk3%n%"t" <- times[3]
samplkl <- list(samplk1, samplk2, samplk3)
samplks <- Networks(samplk1, samplk2, samplk3)

test_that("N() summary with an LM and noncompacted stats", {
  summ.N <- summary(samplks~N(~edges+nodematch("cloisterville"), ~1+t), term.options=list(N.compact_stats=FALSE))
  summ.l <- unlist(lapply(samplkl, function(nw) summary(statnet.common::nonsimp_update.formula(~edges+nodematch("cloisterville"), nw~., from.new="nw"))))
  expect_equivalent(summ.l, summ.N)
})


test_that("N() summary with offset and compacted stats", {
  summ.N <- summary(samplks~N(~edges, offset=~t))
  summ.l <- sapply(samplkl, function(nw) summary(statnet.common::nonsimp_update.formula(~edges, nw~., from.new="nw")))
  summ.l <- c(sum(summ.l), sum(summ.l*1:3))
  expect_equivalent(summ.l, summ.N)
})

pl <- lapply(samplkl, function(nw) ergmMPLE(statnet.common::nonsimp_update.formula(~edges+nodematch("cloisterville"), nw~., from.new="nw")))

for(N.compact_stats in c(FALSE,TRUE)){
  testlab <- if(N.compact_stats) "compacted stats" else "noncompacted stats"
  
  test_that(paste("N() estimation with an LM and", testlab),{
    # Three networks, jointly.
    ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+t), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats)))
  
    nr <- sapply(lapply(pl, `[[`, "response"),length)
  
    y <- unlist(lapply(pl, `[[`, "response"))
    x <- do.call(rbind,lapply(pl, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,t=rep(times, nr)))
    w <- unlist(lapply(pl, `[[`, "weights"))
    glmfit <- glm(y~t*nodematch.cloisterville,data=x,weights=w,family="binomial")
    expect_equivalent(coef(glmfit),coef(ergmfit))
  })
  
  test_that(paste("N() estimation with an LM, subsetting, and", testlab),{
  # Ignore second (test subset).
    ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+t, subset=quote(t!=2)), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats)))
    
    pl2 <- pl[-2]
    times2 <- times[-2]
    
    nr <- sapply(lapply(pl2, `[[`, "response"),length)
    
    y <- unlist(lapply(pl2, `[[`, "response"))
    x <- do.call(rbind,lapply(pl2, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,t=rep(times2, nr)))
    w <- unlist(lapply(pl2, `[[`, "weights"))
    glmfit <- glm(y~t*nodematch.cloisterville,data=x,weights=w,family="binomial")
    expect_equivalent(coef(glmfit),coef(ergmfit))
  })
    
  test_that(paste("N() estimation with an LM, subsetting down to only 1 network, and", testlab),{
    # Ignore first and third (test subset).
    ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+t, subset=quote(t==2)), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats)))
    
    pl2 <- pl[2]
    times2 <- times[2]
    
    nr <- sapply(lapply(pl2, `[[`, "response"),length)
    
    y <- unlist(lapply(pl2, `[[`, "response"))
    x <- do.call(rbind,lapply(pl2, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,t=rep(times2, nr)))
    w <- unlist(lapply(pl2, `[[`, "weights"))
    glmfit <- glm(y~t*nodematch.cloisterville,data=x,weights=w,family="binomial")
    # Note that unlike glmift, ergm cannot detect nonidentifiability at
    # this time.
    if(N.compact_stats){
      expect_equivalent(na.omit(coef(glmfit)),na.omit(coef(ergmfit)))
    }else{
      expect_equivalent(na.omit(coef(glmfit)),
                        c(matrix(c(1,2,0,0,
                                   0,0,1,2),2,4,byrow=TRUE)%*%coef(ergmfit)))
    }
  })
}
