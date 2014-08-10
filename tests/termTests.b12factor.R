library(ergm)



bipnet<-network.initialize(4,bipartite=2,directed=FALSE)
add.edges(bipnet,tail=c(1),head=c(4))
# set an attribute on the the vertices of the second mode
set.vertex.attribute(bipnet,'felines',c('cat','tiger'),v=3:4)

tryCatch(summary(bipnet~b1factor('felines')),error=function(e)cat("expected error message \n"))


if(summary(bipnet~b2factor('felines')) !=1)
	stop("b2factor error A")


