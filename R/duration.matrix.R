#############################################################################
# The <duration.matrix> function computes and returns a duration matrix, via
# <DurationMatrix.C>, for a given network list.
#
# --PARAMETERS--
#   gsim: a network list, as returned by <simulate.stergm>
#
# --RETURNED --
#   allties: the duration matrix as a dataframe, with rows as the ties  and
#   with the 6 following columns:
#     "Ego"        : the ego node in the tie
#     "Alter"      : the alter node in the tie
#     "Start"      : the start of the tie
#     "End"        : the end of the tie, or N if the end has been censored,
#                    where N is the maximum time step + 1
#     "Noncensored": an indicator for whether the end of the tie was observed;
#                    0 indicates 'no, it was censored', 1 indicates 'yes,
#                    it was observed'
#     "duration"   : the duration of the tie as "End"-"Start"; note that this
#                    will include the duration of censored ties
#
#############################################################################

## FIXME:  gsim object is currently an obsolete version of network.list but
## it should be a networkDynamic object
duration.matrix <- function(gsim) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.list") {
    stop("This function requires that the argument be a network.list")
  }
  cha <- gsim$changed[,c(1,2,3)]
  N <- max(cha[,1])+1
  Nbip <- g0%n%"bipartite"
  if(is.null(Nbip)) Nbip<-0
  Ntot <- g0%n%"n"
  edges <- as.edgelist(g0)
  nedge <- nrow(edges)
  # Workaround --- if a network is undirected (or bipartite), force tails<heads.
  if(!is.directed(g0) && nedge) edges<-t(apply(edges,1,sort))
  nchange <- nrow(cha)

  allties <- .C("DurationMatrix", as.integer(nedge), as.integer(edges), 
                 as.integer(N), as.integer(Ntot), 
                 as.integer(nchange), as.integer(cha),
                 duration = as.integer(rep(0,5*(nedge+nchange))),
                 PACKAGE = "ergm")$duration
  allties <- matrix(allties, ncol=5)
  colnames(allties) <- c("Ego", "Alter", "Start", "End", "Noncensored")
  allties <- allties[allties[,1]!=0,] # Get rid of unused rows
  allties[allties==-1] <- N # Edges that didn't end will be censored at N
  
  as.data.frame(allties <- cbind(allties, duration=allties[,4]-allties[,3]))
}
