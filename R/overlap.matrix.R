#############################################################################
# The <overlap.matrix> function calculates a summary of concurrent
# partnerships in a bipartite network list ; this includes the overlap
# durations and is computed via <OverlapsDurations.C> function; it's
# assumed that the smaller-numbered nodes are "females" and the larger-
# numbered ones are "males"
#
# NB: I do not know what would happen if there were instances of a node having
# three edges --DH
#
# --PARAMETERS--
#   gsim       : a network list, which is assumed to be bipartite; the elements
#                of interest are the original network and the 3-column matrix
#                giving the times and nodes of each edge-toggle starting from
#                the original
#   maxoverlaps: the maximum number of overlaps to provide space for; if more
#                overlaps are observed than space is made for, an error will
#                be printed and the first 'maxoverlaps' overlaps are returned
#                
#
# --RETURNED--
#   a 10-column data frame where each row describes one instance of an
#   overlapping, or concurrent, partner-pair; the columns are:  
#     1: the female node number of the first partnership
#     2: the male node number of the first partnership
#     3: the female nocde number of the second partnership
#     4: the male node number of the second partnership
#     5: the start time of the first partnership
#     6: the start time of the second partnership
#     7: the time when the concurrency ended, i.e., the time when one 
#        of the two edges ended.  If both edges continue past the end of
#        the simulation, then this value will equal N, where N is
#        1 more than the largest time observed for any network change.
#     8: an inidicator of which edge ended first, where
#            0 implies neither, because of truncation
#            1 implies the first
#            2 implies the second
#            3 implies both, because the ended at the same time
#     9: the duration of the concurrency, as  (column 7 - max(col 5, col 6))
#    10: the type of the overlap, as one of:
#            "transitional":  the first relationship to start is also the
#                             first to end 
#             "embedded"   :  the second relationship starts and ends while
#                             the first continues.
#             "truncated"   : we can't identify the overlap type due to
#                             censoring at either the beginning or the end
#         Note that any overlaps involving ties are considered embedded as
#         long as both the start and the end of the overlap are observed.
#         If either the start of the end is unobserved, the concurrency is
#         considered truncated even if the observed end is a tie.
#         Should this be changed?  
#         Should any ties be called embedded, even if the start is 0 or the end is N?
#
#############################################################################

overlap.matrix <- function(gsim, maxoverlaps=1000000) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.list") {
    stop("This function requires that the argument be a network.list")
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


