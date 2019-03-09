#  File tests/ergm_term_doc_test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
# tests to veryify the structure of the ergm-terms documentation
# this executes functions similar to what would be run when building the ergm=term=crossRef.Rmd file
# and verifys that the structure of ergm-terms.Rd matches what the parser expects
library(ergm)
termBlock<-ergm:::.extractTermBlock()
items<-ergm:::.extractTags(termBlock,"\\item")
terms<-lapply(items,ergm:::.extractTerms)
terms<-unlist(terms,recursive=FALSE)

# this function is defined in build_term_index.R
# it checks certain known assumptions about term tags, like each term must be either binary or valued, etc
ergm:::.checkTermDocs(terms)

# crude checks for search.ergmTerms are in the search.ergmTerms man page

# expect to find at least eight terms mentioning triangles
found<-search.ergmTerms('triangle')
if(length(found)<8){
  stop(' search.ergmTerms unexpectly found less than 8 terms mentioning triangles')
}

found<-search.ergmTerms(categories = 'bipartite')
if(length(found)<20){
  stop(' search.ergmTerms unexpectly found less than 20 terms with the category "bipartite"')
}