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
