#  File tests/utils_tests.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################


# ------------- tests for the update.network function -----

require(ergm)
aaa <- network.initialize(10,directed=TRUE,loops=TRUE)
aaa %v% 'race' <- rep(c('B','W'),10)
aaa %v% 'letters'<-rep(LETTERS,10)
aaa %v% 'numbers'<-runif(10)
aaa %v% 'astructure'<-rep(list(first="A",second="B"),10)
aaa %v% 'inf.status' <-rbinom(10,1,0.2)
aaa %n% 'verymeta' <- 'so meta'
aaa %n% 'verybeta' <- 'so beta'
aaa[2,1]<-1  # add a single edge

# a random matrix
amat<-structure(c(1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0,  0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1,  0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0,  1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,  0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0), .Dim = c(10L,  10L))

# and edgelist version of the same matrix
ael <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L,  4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L,  7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 10L, 10L,  10L, 10L, 10L, 10L, 1L, 2L, 4L, 5L, 6L, 8L, 3L, 8L, 9L, 6L, 9L,  10L, 1L, 3L, 4L, 6L, 10L, 2L, 3L, 5L, 7L, 8L, 9L, 10L, 3L, 4L,  5L, 7L, 8L, 3L, 4L, 5L, 7L, 8L, 10L, 1L, 4L, 5L, 8L, 9L, 1L,  2L, 1L, 2L, 4L, 5L, 6L, 7L), .Dim = c(48L, 2L))

bb<-update(aaa,new = amat,matrix.type = 'adjacency')

# correct number of edges created
if(!all(as.matrix(bb)==amat)){
  stop("update.network did not create correct edges corresponding to input adjacency matrix")
}

# edges removed
if(bb[2,1]!=0){
  stop("update.network did not correctly remove edges in input network")
}

# original unmodified
if(network.edgecount(aaa) >1){
  stop("update.network modified its input argument")
}

# attributes preserved
if (!all(c(
aaa%v%'race'==bb%v%'race',
aaa%v%'letters'==bb%v%'letters',
aaa%v%'numbers'==bb%v%'numbers',
aaa%v%'astructure'==bb%v%'astructure',
aaa %n% 'verymeta' == bb %n% 'verymeta',
aaa %n% 'verybeta' == bb %n% 'verybeta'
))){
  stop("update.network did not copy network and vertex attributes correctly")
}

# flags preserved
if(!has.loops(bb)){
  stop("network flags were not copied correctly by update.network")
}

# try creating from an edgeslist

ccc<-update(aaa,ael,matrix.type='edgelist')
if(!all(as.matrix(ccc)==amat)){
  stop('update.network did not create edges corresponding to input edgelist matrix')
}

