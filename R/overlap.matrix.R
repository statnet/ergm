# overlap.matrix is a wrapper for the OverlapDuration function in C.
# Its job is to look at a network.series object (of which the only important
# parts are the original network and a 3-column matrix giving the times
# and nodes of each edge-toggle starting from the original) and return
# a summary of all instances of concurrent relationships -- i.e., nodes
# having two or more edges.  It's assumed that the network is bipartite
# and that the smaller-numbered nodes are called "females" and the 
# larger-numbered ones are called "males"
#
# I do not know what would happen if there were instances of a node having
# three edges --DH
#
# The summary is in the form of a 10-column data frame.  Each row describes
# one instance of overlapping, or concurrent, relationships.
#   Columns 1, 2, 5:  The female and male node numbers and the start time
#                     of the edge between them
#   Columns 3, 4, 6:  Same as 1, 2, 5 for the second edge
#   Column 7:  The time when the concurrency ended, i.e., the time when one 
#              of the two edges ended.  If both edges continue past the end of
#              the simulation, then this value will equal N, where N is
#              1 more than the largest time observed for any network change.
#   Column 8:  Which of the two edges (1 or 2) was the first to end.  If they
#              ended at the same time, this value will be 3.  If the concurrency
#              continues past the end of the simulation, this value will be 0.
#   Column 9:  The duration of the concurrency.  Easy to calculate:  Just column
#              7 minus the larger of columns 5 and 6.                                                
#   Column 10: The type of the overlap.  There are three types:
#      transitional: The first relationship to start is also the first to end 
#      embedded:  The second relationship starts and ends while the first continues.
#                 Note that any overlaps involving ties are considered embedded
#                 as long as both the start and the end of the overlap are 
#                 observed.  If either the start of the end is unobserved, 
#                 the concurrency is considered truncated even if the observed 
#                 end is a tie.  Should this be changed?  Should any ties be 
#                 called embedded, even if the start is 0 or the end is N?
#      truncated: We can't identify the overlap type due to censoring at either
#                 the beginning or the end

overlap.matrix <- function(gsim, maxoverlaps=100000) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.series") {
    stop("This function requires that the argument be a network.series")
  }
  if(!is.bipartite(gsim$networks)) stop("This function only works for bipartite networks at the moment.")
  cha <- gsim$changed
  cha <- cha[order(cha[,1]),]
  N <- max(cha[,1])+1
  Nfem <- g0%n%"bipartite"
  Ntot <- g0%n%"n"
  edges <- edgelist.ergm(g0)
  nedge <- nrow(edges)
  nchange <- nrow(cha)

  overlap <- .C("OverlapDurations", as.integer(Ntot), as.integer(nedge), as.integer(edges), 
                as.integer(N), as.integer(Nfem), as.integer(Ntot), 
                as.integer(nchange), as.integer(cha),
                as.integer(maxoverlaps),
                omatrix = as.integer(rep(0,8*maxoverlaps)),
                PACKAGE = "ergm")$omatrix
  overlap <- matrix(overlap, ncol=8)
  colnames(overlap) <- c("Fem1", "Male1", "Fem2", "Male2", "start1", "start2", 
                         "endtime", "firsttoend")
  overlap <- overlap[overlap[,1]>0,] # Get rid of unused rows
  overlap[overlap==-1] <- N # Overlaps that didn't end will be censored at N
  duration <- overlap[,7] - apply(overlap[,5:6], 1, max) 
  fts <- apply(overlap[,c("start1","start2")], 1, 
                        function(x) (x[1]<=x[2]) + 2*(x[2]<=x[1]))
  fte <- overlap[,"firsttoend"]
  type <- (((fts == 1 & fte == 2) | (fts == 2 & fte == 1)) 
           + 2*((fts == 3 & fte > 0) | fte == 3 | fte == fts) 
           + 3*(fte==0))
  type[overlap[,"start1"]==0 & overlap[,"start2"]==0] <- 3
  type <- c("transitional","embedded","truncated")[type]
  data.frame(overlap, duration = duration, type=type)
}


