# tests to veryify the structure of the ergm-terms documentation
library(ergm)
termBlock<-ergm:::.extractTermBlock()
items<-ergm:::.extractTags(termBlock,"\\item")
terms<-lapply(items,ergm:::.extractTerms)
terms<-unlist(terms,recursive=FALSE)

# this function is defined in build_term_index.R
ergm:::.checkTermDocs(terms)