
# parse the ergm-terms doc file to create structured data about each set of terms
# so that we can generate indexed documentation



#currently, the content we are interested is in index 491, but need more robust way to discover this
# should be the 3rd section?
.extractTermBlock<-function(rawdoc){
  
  return<-rawdoc[[491]][[2]][[3]]
}

# takes a chunk of text like "   (tag1) (tag2)" and returns just the tags
.extractCats<-function(text){
  openParens<-gregexpr("(",text,fixed=TRUE)[[1]]
  closeParens<-gregexpr(")",text,fixed=TRUE)[[1]]
  lapply(seq_len(length(openParens)),function(i){
    substr(text,openParens[i]+1,closeParens[i]-1)
  })
}

# function to extract parts of an Rd object with specific tags
.extractTags<-function(block,tagvalue){
  indices<-sapply(block,attr,'Rd_tag')==tagvalue
  return(block[indices])
}

# function takes an index item, and extracts the term info from it
.extractTerms<-function(item){
  terms<-list()
  #first part contains the term name(s) and args and tags
  # assume that the parts in \code blocks are the term names, and deliniate multiple terms
  termIndices<-which(sapply(item[[1]],attr,'Rd_tag')=="\\code")
  if(length(termIndices)==0){
    warning('no terms found for item')
  }
  for (t in seq_len(length(termIndices))){
    rawterm<-as.character(item[[1]][[termIndices[t]]])
    # debug
    #print(rawterm)
    # the part before the ( is the term, part after is args
    funparts<-strsplit(rawterm,'(',fixed=TRUE)[[1]]
    term.name<-funparts[1]
    term.id<-paste('term',term.name,t,sep="_") # for document cross refs
    term.args<-substr(funparts[2],start=0,stop=nchar(funparts[2])-1)
    # now try to find the category tags
    if(((termIndices[t]+1) > length(item[[1]])) | (t<length(termIndices) & (termIndices[t]+1==termIndices[t+1]))){
      warning('could not locate category tags for term ',term.name)
    } else {
       cats<-.extractCats(as.character(item[[1]][[termIndices[t]+1]]))
    }
    # try to extract short definition from the first sentance
    short.desc<-paste(item[[2]][[1]],sep='',collapse=' ')
    short.desc<-sub("\n"," ",short.desc)
    short.desc<-substring(short.desc,0,nchar(short.desc)-1)
    
    terms[[t]]<-list(term.id=term.id,term.name=term.name,term.args=term.args,categories=cats,short.desc=short.desc, description.rd=item[[2]])
  }
  
    
  # second part has the definition text and example 
  return(terms)
}







# output listings of terms, grouped by category
.termToc<-function(terms){
  cats<-unique(unlist(sapply(terms,'[[','categories')))
  # print out a category header
  cat("Jump to category:")
  for(cat in cats){
    cat("<a href='#cat_",cat,"'>",cat,"</a> ",sep='')
  }
  for (cat in cats){
    cat("<h3> <a id='cat_",cat,"'>",cat,"</a></h3>",sep='')
    #find ones with matching cats
    matchingTerms<-sapply(terms,function(term){
      any(term$categories==cat)
    })
    cat("<ul>\n")
    lapply(terms[matchingTerms],function(term){
      cat("<li><a href='#",term$term.id,"'>",term$term.name,"</a> : ",term$short.desc,"</li>\n",sep='')
    })
    cat("</ul>\n")
  }
}

.termTable<-function(terms){
  cat("<table border=1 cellpadding='8'>\n")
  cat("<tr>><th>Description</th><th>Categories</th></tr>\n")
  for (term in terms){
    cat("<tr><td><a id='",term$term.id,"'>",term$term.name,"(",ifelse(is.na(term$term.args),'',term$term.args),")</a><br>",sep='')
    cat(capture.output(tools::Rd2HTML(term$description.rd,fragment=TRUE)),"</p></td><td>",paste(term$categories,collapse=", "),"</td></tr>")
  }
  cat("</table>")
}

.termMatrix<-function(terms){
  cats<-unique(unlist(sapply(terms,'[[','categories')))
  # figure out which terms are members of each cat
  membership<-lapply(cats,function(cat){
    # return checkmark for terms that match cat, otherwise blank
    sapply(terms,function(term){
      if(any(term$categories==cat)){
        return('&#10004;')
      } else {
        return('')
      }
    })
  })
  # generate the html table
  cat("<table border=1 cellpadding='8'>\n")
  cat("<tr><th>Term name</th><th>",paste(cats,collapse='</th><th>'),"</th></tr>\n",sep='')
  for (t in seq_along(terms)){
    term<-terms[[t]]
    cat("<tr><td><a href='#",term$term.id,"'>",term$term.name,"</a></td>",sep="")
    for(c in seq_along(cats)){
      cat("<td>",membership[[c]][[t]],"</td>")
    }
    cat("</tr>",sep='')
  }
  cat("</table>")
}

