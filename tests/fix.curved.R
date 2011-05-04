library(ergm)
data(sampson)
out<-fix.curved(samplike~edges+gwnsp(alpha=.5,fixed=TRUE)+gwesp(alpha=.5,fixed=FALSE)+gwodegree(decay=.5,fixed=FALSE)+edges,c(1:7))
stopifnot(out$formula==(samplike ~ edges + gwnsp(alpha = 0.5, fixed = TRUE) + gwesp(alpha = 4L, fixed = TRUE) + gwodegree(decay = 6L, fixed = TRUE) + edges),
          all(out$theta==c(1,2,3,5,7)))
