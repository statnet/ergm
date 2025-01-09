#  File R/build_term_index.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

.termTable<-function(terms,merge.names=FALSE){
  cat("<table border=1 cellpadding='8'>\n")
  cat("<tr><th>Description</th><th>Categories</th></tr>\n")
  for (term in terms) {
    u <- sapply(term$usages, '[[', 'usage')
    if (merge.names) {
      usages <- c(paste(u, collapse=' '))
    } else {
      usages <- c()
      for (alias in unique(gsub(' *\\(.*', '', u))) {
        usages <- c(usages, paste(u[which(startsWith(u, alias))], collapse=' '))
      }
    }

    for (usage in usages) {
      cat(sprintf('<tr><td><a id="%s">%s</a><br /><em>%s</em>: %s</p></td><td>%s</td></tr>\n',
        term$link, usage, term$title, term$description, paste(term$concepts, collapse=', ')))
    }
  }
  cat("</table>")
}
  
# go through the structure parsed from term documentation to 
# ensure that it meets our expectations
.checkTermDocs <-function(terms){
  for (term in terms){
    
    # every term must include at least one of 'directed', 'undirected', or 'operator'
    if (!any(c('directed','undirected','operator')%in%term$concepts)){
      stop('the term ',term$name,' must be marked as directed and/or undirected in the documentation')
    }
    # every term must include either 'valued' or 'binary'
    # check that there is a visable init function defined for the term
    # some terms have both valued an binary forms
    if ('valued'%in%term$concepts){
      if(!is.function(eval(locate_prefixed_function(term$name, 'InitWtErgmTerm')))){
        stop('unable to locate an InitWtErgmTerm function defined for weighted term ',term$name,' in documentation')
      }
    } 
    if ('binary'%in%term$concepts){
      if(!is.function(eval(locate_prefixed_function(term$name, 'InitErgmTerm')))){
        stop('unable to locate an InitErgmTerm function defined for term ',term$name,' in documentation')
      }
    }
    if (!any(c('binary','valued')%in%term$concepts)) {
      stop('the term ',term$name,' is not marked as binary or valued in the documentation')
    }
  }
}
