library(ergm)

data(florentine)

tryCatch(ergm(flomarriage~edges+degree(1)+offset(degree(0)),target.stats=summary(flomarriage~edges+degree(1)),
              offset.coef=-0.25,control=control.ergm(init=c(-1.47,0.462,-0.25))) ,
         error=function(e) stop('error in target + offset test'), 
         warning=function(w) stop('unexpected warning in target + offset test'))

