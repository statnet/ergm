library(ergm)

data(florentine)
text <- capture.output(out <- simulate(flomarriage~edges+degree(0)+absdiff("wealth")+passthrough(~edges+degree(0)+absdiff("wealth"))+submodel.test(~edges+degree(0)+absdiff("wealth"))+summary.test(~edges+degree(0)+absdiff("wealth")), statsonly=TRUE, nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1)))
text.out <- matrix(scan(textConnection(paste(text, collapse=""))),byrow=TRUE,ncol=3)
text.out <- text.out[nrow(text.out)-nrow(out)+seq_len(nrow(out)),]

stopifnot(all(out[,1:3]==out[,4:6]),all(out[,1:3]==out[,7:9]),all(out[,1:3]==text.out))

data(sampson)
g <- samplike%v%"group"
sameg <- outer(g,g,"==")

out <- simulate(samplike~nodematch("group")+odegree(0:5, by="group", homophily=TRUE)+idegree(0:5, by="group", homophily=TRUE)+localtriangle(sameg)+
                  NodematchFilter(~edges+odegree(0:5)+idegree(0:5)+triangle,"group")+
                  F(~edges+odegree(0:5)+idegree(0:5)+triangle,~nodematch("group")), statsonly=TRUE, nsim=20, control=control.simulate.formula(MCMC.burnin=0, MCMC.interval=1))

stopifnot(all(out[,1:14]==out[,15:28]), all(out[,1:14]==out[,29:42]))
