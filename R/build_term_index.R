
# parse the ergm-terms doc file to create structured data about each set of terms
# so that we can generate indexed documentation




# grab a the relevent section of .Rd document structure
# based on the location of a specific comment tag
.extractTermBlock<-function(){
  # query the install documentation
  rawdoc<-tools::Rd_db('ergm')$'ergm-terms.Rd'
  # find the tag indicating where the term definitions begin
  
  for (i in rawdoc) {
    if (length(i) > 1) {
      defIndex<-grep('beginTerms',i[[2]])
      if (length(defIndex)>0) return(i[[2]][[defIndex]])
    }
  }
  
  stop("Error extracting term definitions")
}

# takes a chunk of text like "   (tag1) (tag2)" and returns just the tags
.extractCats<-function(text){
  openParens<-gregexpr("(",text,fixed=TRUE)[[1]]
  closeParens<-gregexpr(")",text,fixed=TRUE)[[1]]
  lapply(seq_len(length(openParens)),function(i){
    substr(text,openParens[i]+1,closeParens[i]-1)
  })
}

# utility function to extract parts of an Rd object with specific tags
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
       if(length(cats)==1 && cats[[1]]==''){
         warning('failed to extract category tags for term ',term.name,'. Are all of the tags enclosed in () and listed on the same line as the term name?')
       }
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

.termTable<-function(terms,merge.names=FALSE){
  cat("<table border=1 cellpadding='8'>\n")
  cat("<tr>><th>Description</th><th>Categories</th></tr>\n")
  if (merge.names){
    # merge all the terms with the same names
    termNames<-sapply(terms,'[[','term.name')
    for (term.name in unique(termNames)){
      matchedTerms<-terms[which(termNames==term.name)]
      cat("<tr><td><a id='",matchedTerms[[1]]$term.id,"'>")
      for(term in matchedTerms){
        cat(' ',term$term.name,"(",ifelse(is.na(term$term.args),'',term$term.args),")",sep='')
      }
      cat('</a><br>')
      cat(capture.output(tools::Rd2HTML(matchedTerms[[1]]$description.rd,fragment=TRUE)),"</p></td>")
      cat("<td>",paste(unique(unlist(lapply(matchedTerms,'[[','categories'))),collapse=", "),"</td></tr>")
    }
  } else {
    for (term in terms){
      cat("<tr><td><a id='",term$term.id,"'>",term$term.name,"(",ifelse(is.na(term$term.args),'',term$term.args),")</a><br>",sep='')
      cat(capture.output(tools::Rd2HTML(term$description.rd,fragment=TRUE)),"</p></td><td>",paste(term$categories,collapse=", "),"</td></tr>")
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
  cats<-unique(unlist(sapply(terms,'[[','categories')))
  if(is.null(categories)){
    categories<-cats
  } else { 
    # check that not requesting something that doesn't exist
    if (any(!categories%in%cats)){
      stop("requested column name does not appear in documentation category tags")
    }
  }
  
  # figure out which terms are members of each cat
  membership<-lapply(categories,function(cat){
    # return checkmark for terms that match cat, otherwise blank
    sapply(terms,function(term){
      if(any(term$categories==cat)){
        return('&#10004;')
      } else {
        return('')
      }
    })
  })
  
  # figure out which terms should be included
  if(!is.null(only.include)){
    included<-sapply(terms,function(term){
      if(any(term$categories%in%only.include)){
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
    cat("<tr><td align='center'><a href='#",term$term.id,"'>",term$term.name,"</a></td>",sep="")
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
    
    # every term must include at least one of 'directed' and 'undirected'
    if (!any(c('directed','undirected')%in%term$categories)){
      stop('the term ',term$term.name,' must be marked as directed and/or undirected in the documentation')
    }
    # every term must include either 'valued' or 'binary'
    # check that there is a visable init function defined for the term
    # some terms have both valued an binary forms
    if ('valued'%in%term$categories){
      if(!is.function(get(paste('InitWtErgmTerm',term$term.name,sep='.')))){
        stop('unable to locate an InitWtErgmTerm function defined for weighted term ',term$term.name,' in documentation')
      }
    } 
    if ('binary'%in%term$categories){
      if(!is.function(get(paste('InitErgmTerm',term$term.name,sep='.')))){
        stop('unable to locate an InitErgmTerm function defined for term ',term$term.name,' in documentation')
      }
    }
    if (!any(c('binary','valued')%in%term$categories)) {
      stop('the term ',term$term.name,' is not marked as binary or valued in the documentation')
    }
  }
}

# function to look up the set of terms applicable for a specific network

search.ergmTerms<-function(keyword,net,categories,name){
  
  if (!missing(net)){
    if(!is.network(net)){
      stop("the 'net' argument must be the network argument that applicable terms are to be searched for")
    }
  }
  
  termBlock<-.extractTermBlock()
  items<-.extractTags(termBlock,"\\item")
  terms<-lapply(items,.extractTerms)
  terms<-unlist(terms,recursive=FALSE)
  
  
  
  if(missing(categories)){
    categories<-character(0)
  }
  if (!missing(net)){
    if(is.directed(net)){
      categories<-c(categories,'directed')
    } else {
      categories<-c(categories,'undirected')
    }
    if(is.bipartite(net)){
      categories<-c(categories,'bipartite')
    } 
  }
  found<-rep(TRUE,length(terms))
  
  # if name is specified, restrict to terms with that name
  if(!missing(name)){
    for (t in seq_along(terms)){
      term<-terms[[t]]
      found[t]<-any(term$term.name==name)
    }
  }
  
  # restrict by categories
  for (t in which(found)){
    term<-terms[[t]]
    if(!all(categories%in%term$categories)){
      found[t]<-FALSE
    }
  }
  
  # next (optionally) restrict by keyword matches
  if (!missing(keyword)){
    for (t in which(found)){
      term<-terms[[t]]
      # if we don't find the keyword in the text grep, mark it as false
      descText<-capture.output(tools::Rd2txt(term$description.rd,fragment=TRUE))
      if(length(grep(keyword,descText,ignore.case=TRUE))==0){
        found[t]<-FALSE
      } 
    }
  }
  
  
  # if term name was specified, print out all the matching terms
  # otherwise,  loop over the remaining terms to format output as condensed
  output<-list()
  if(!missing(name)){
    if(sum(found)==0){
      cat("No terms named '",name,"' were found. Try searching with keyword='",name,"'instead.",sep='')
    } else {
      cat("Definitions for term(s) ",name,":\n")
      for (t in which(found)){
        term<-terms[[t]]
        output<-c(output,term)
        cat(paste(term$term.name,"(",ifelse(is.na(term$term.args),'',term$term.args),")\n    ",sep=''), capture.output(tools::Rd2txt(term$description.rd,fragment=TRUE)),"\n    Categories:",paste(term$categories,collapse=', '),"\n\n")
      }
    }
  }else{
    for (t in which(found)){
      term<-terms[[t]]
      outText<-paste(term$term.name,"(",ifelse(is.na(term$term.args),'',term$term.args),")\n    ",term$short.desc,"\n",sep='')
      outText<-sub("\t",'',outText)
      outText<-gsub(" +",' ',outText)
      output<-c(output,outText)
    }
    cat("Found ",length(output)," matching ergm terms:\n")
    cat(paste(output,collapse='\n'))
  }
  invisible(output)
}

