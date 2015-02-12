#  File tests/constrain_degrees.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(statnet.common)

# fixedas

opttest({
			
	library(ergm)
	
	net1 <- network(10,directed=FALSE,density=0.2)
	
	el1 <- as.edgelist(net1)
	
# both present and absent
	present <- as.edgelist(el1[c(1,2),],n=10,directed=FALSE)
	
	absent <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)
	
	while (any(as.data.frame(t(absent)) %in% as.data.frame(t(el1)))){
		{absent <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)}
	}
	
	t1 <- ergm(net1~edges,constraint=~fixedas(present=present,absent=absent))
	
	s1 <-simulate(t1,1000)
	
	# check if all the simulated network have 'present' edges
	stopifnot(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))
	
	# check if all the simulated network do not have 'absent' edges
	stopifnot(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))
	
	
# only present
			present <- as.edgelist(el1[c(1,2),],n=10,directed=FALSE)
			
			t1 <- ergm(net1~edges,constraint=~fixedas(present=present))
			
			s1 <-simulate(t1,1000)
			
			stopifnot(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

			
# only absent
			
			absent <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)
			
			while (any(as.data.frame(t(absent)) %in% as.data.frame(t(el1)))){
				{absent <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)}
			}
			
			t1 <- ergm(net1~edges,constraint=~fixedas(absent=absent))
			
			s1 <-simulate(t1,1000)
			
			stopifnot(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))
			
# input is network
			
			present <- network.initialize(10,directed=FALSE)
			add.edges(x=present,tail=el1[c(1,2),1],head=el1[c(1,2),2])
			
			absent.el <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)
			
			while (any(as.data.frame(t(absent.el)) %in% as.data.frame(t(el1)))){
				{absent.el <- as.edgelist(matrix(sample(1:10,4,replace=F),2,2),n=10,directed=FALSE)}
			}
			
			absent <- network.initialize(10,directed=FALSE)
			
			add.edges(x=absent,tail=absent.el[,1],head=absent.el[,2])
			
			
			t1 <- ergm(net1~edges,constraint=~fixedas(present=present, absent=absent))
			
			s1 <-simulate(t1,1000)
			
			stopifnot(all(sapply(s1,function(x)as.data.frame(t(as.edgelist(present))) %in% as.data.frame(t(as.edgelist(x))))))
			
			stopifnot(all(!sapply(s1,function(x)as.data.frame(t(as.edgelist(absent))) %in% as.data.frame(t(as.edgelist(x))))))
			
			
		})
		
		
# fixallbut

opttest({
			
	library(ergm)
	
	net1 <- network(10,directed=FALSE,density=0.5)

	free.dyads <- as.edgelist(matrix(sample(1:10,8,replace=F),4,2),n=10,directed=FALSE)
	
	t1 <- ergm(net1~edges,constraint=~fixallbut(free.dyads=free.dyads))
	
	s1 <-simulate(t1,1000)
	
	fixed.dyads.state <- net1[as.edgelist(invert.network(network.update(net1,free.dyads,matrix.type="edgelist")))]
	
	stopifnot(all(sapply(s1,function(x) all.equal(x[as.edgelist(invert.network(network.update(x,free.dyads,matrix.type="edgelist")))],fixed.dyads.state))))
	
	
})