check.ERGMterm.situation <- function(nw, directed=NULL, bipartite=NULL) {
  # Check to see whether there are any conditions that preclude
  # the use of this ERGM term. (Called by InitERGMterm.xxx functions)
  message <- NULL
  if (!is.null(directed) && directed != is.directed(nw)) {
    message <- paste("networks with directed==",is.directed(nw),sep="")
  }
  if (!is.null(bipartite) && bipartite != (nw %v% "bipartite" > 0)) {
    message <- paste("networks with bipartite", 
                     ifelse(nw %v% "bipartite">0,"!= FALSE", " > 0"))
  }
  if (!is.null(message)) {
    sc <- sys.calls()
    fname <- ifelse(length(sc)>1, as.character(sc[[length(sc)-1]]), NULL)
    fname <- substring(fname,14) # get rid of leading "InitERGMterm."
    stop(paste("The ERGM term",fnames,"may not be used with",message))
  }
}

