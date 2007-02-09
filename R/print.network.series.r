"print.network.series" <-
  function (x, ..., wmt = which.matrix.type(x$networks[[1]])) 
{
  g <- x$networks[[1]]
  cat("Number of Networks:",length(x$networks),"\n")
  f.out <- lapply(x$formula,as.character)
  if(length(f.out)>2)
    {
      if(length(f.out[[3]])>2)
        f.output <- paste(f.out[[2]],f.out[[1]],paste(f.out[[3]][c(2,1,3)],collapse=" "))
      else
        f.output <- paste(f.out[[2]],f.out[[1]],f.out[[3]])
    }
  else
    {
      if(length(f.out[[2]])>2)
        f.output <- paste(f.out[[1]],paste(f.out[[2]][c(2,1,3)],collapse=" "))
      else
        f.output <- paste(f.out[[1]],f.out[[2]])
    }
  cat("From:",f.output,"\n")
  cat("parameter values:",x$coef,"\n")
  cat("Networks attributes:\n")
  for (i in 1:length(g$gal)) {
    if (names(g$gal)[i] == "n") {
      attributeName <- "nodes"
      attributeValue <- g$gal[[i]]
    }
    else if (names(g$gal)[i] == "mnext") {
      if (!is.directed(g)) 
        attributeName <- "edges"
      else attributeName <- "arcs"
      attributeValue <- network.edgecount(g)
    }
    else {
      attributeName <- names(g$gal)[i]
      attributeValue <- g$gal[[i]]
    }
    cat("  ", attributeName, "=", attributeValue, "\n")
  }
  cat("\n", wmt, " matrix:\n")
  print(as.matrix.network(g, matrix.type=wmt))
  invisible(g)
}
