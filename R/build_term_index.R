#  File R/build_term_index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

# output listings of terms, grouped by category
.termToc<-function(terms){
  cats<-unique(unlist(sapply(terms,'[[','concepts')))
  # print out a category header
  cat("Jump to category:")
  for(cat in cats){
    cat("<a href='#cat_",cat,"'>",cat,"</a> ",sep='')
  }
  for (cat in cats){
    cat("<h3> <a id='cat_",cat,"'>",cat,"</a></h3>",sep='')
    #find ones with matching cats
    matchingTerms<-sapply(terms,function(term){
      any(term$concepts==cat)
    })
    cat("<ul>\n")
    lapply(terms[matchingTerms],function(term){
      cat("<li><a href='#",term$link,"'>",term$name,"</a> : ",term$title,"</li>\n",sep='')
    })
    cat("</ul>\n")
  }
}

.termTable<-function(terms,merge.names=FALSE){
  cat("<table border=1 cellpadding='8'>\n")
  cat("<tr>><th>Description</th><th>Categories</th></tr>\n")
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
      cat(sprintf('<tr><td><a id="%s">%s</a><br /><em>%s</em>: %s</p></td><td>%s</td></tr>',
        term$link, usage, term$title, term$description, paste(term$concepts, collapse=', ')))
    }
  }
  cat("</table>")
}

# terms : a list structure of the documentation data
# categores : an optional vector of column names to print and include
# only.include : an optional vector of categories, only terms that match the category will be printed 
.termMatrix<-function(terms,categories=NULL,only.include=NULL){
  
  # if list of categories not supplied, generate it
  # otherwise, use the categories (and column order) provided
  cats<-unique(unlist(sapply(terms,'[[','concepts')))
  if(is.null(categories)){
    categories<-cats
  }##  else { 
  ##   # check that not requesting something that doesn't exist
  ##   if (any(!categories%in%cats)){
  ##     stop("requested column name does not appear in documentation category tags")
  ##   }
  ## }
  
  # figure out which terms are members of each cat
  membership<-lapply(categories,function(cat){
    # return checkmark for terms that match cat, otherwise blank
    sapply(terms,function(term){
      if(any(term$concepts==cat)){
        return('&#10004;')
      } else {
        return('')
      }
    })
  })
  
  # figure out which terms should be included
  if(!is.null(only.include)){
    included<-sapply(terms,function(term){
      if(any(term$concepts%in%only.include)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    terms<-terms[included]
	# fix a bug , the membership didn't update for only.include
	membership <- lapply(membership,"[",included)
  }
  
  
  
  # generate the html table
  cat("<table>\n")
  cat("<tr><th>Term name</th><th>",paste(categories,collapse='</th><th>'),"</th></tr>\n",sep='')
  for (t in seq_along(terms)){
    term<-terms[[t]]
    cat("<tr><td align='center'><a href='#",term$link,"'>",term$name,"</a></td>",sep="")
    for(c in seq_along(categories)){
      cat("<td align='center'>",membership[[c]][[t]],"</td>")
    }
    cat("</tr>",sep='')
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
