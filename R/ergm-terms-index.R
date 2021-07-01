#  File R/ergm-terms-index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons
#######################################################################

SUPPORTED_TERM_TYPES <- c('Term', 'Constraint', 'Reference')
SUPPORTED_TERM_TYPE_REGEX <- sprintf('-ergm(%s)(.Rd)?', paste(SUPPORTED_TERM_TYPES, collapse='|'))
CONCEPT_ABBREVIATIONS <- list('directed'='dir', 'dyad-independent'='dyad-indep', 'quantitative nodal attribute'='quant nodal attr',
    'undirected'='undir', 'binary'='bin', 'valued'='val', 'categorical nodal attribute'='cat nodal attr',
    'curved'='curved', 'triad-related'='triad rel', 'operator'='op', 'bipartite'='bipartite',
    'frequently-used'='freq', 'non-negative'='non-neg')

DISPLAY_TEXT_INDEX_MAX_WIDTHS <- list('Term'=25, 'Pkg'=5, 'Description'=33, 'Concepts'=12)
DISPLAY_TEXT_MAX_WIDTH <- sum(unlist(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) + length(DISPLAY_TEXT_INDEX_MAX_WIDTHS) - 1
DISPLAY_LATEX_INDEX_PCT_WIDTHS <- c(0.3, 0.05, 0.5, 0.1)
DISPLAY_LATEX_TOC_PCT_WIDTHS <- function(n_concepts) c(2.4, rep(.7, n_concepts))

# Return the index entry for a single term in the new format
.parseTerm <- function(name, pkg, pkg_name) {
  doc <- pkg[[name]]
  name <- substr(name, 1, nchar(name) - 3)
  tags <- sapply(doc, attr, 'Rd_tag')
  comps <- strsplit(name, '-')[[1]]

  concepts <- doc[tags == '\\concept'] %>% unlist
  type <- comps[length(comps)]
  raw_usage <- doc[tags == '\\usage'] %>% 
    unlist %>% 
    paste(collapse="") %>%
    gsub("\n *# *"," ", .) %>%
    strsplit("\n") %>%
    .[[1]] %>%
    trimws()

  if (type %in% c('ergmTerm', 'ergmProposal')) {
    usages <- list()
    for (usage_line in regmatches(raw_usage, regexec("^(binary|valued): *(.+)$", raw_usage))) {
      usages[[length(usages) + 1]] <- list('type'=usage_line[2], 'usage'=usage_line[3])
      concepts <- c(concepts, usage_line[2])
    }
  } else {
    usages <- lapply(raw_usage, function(u) list(type=NULL, usage=u))
  }

  ret <- list(
    link=name,
    name=comps[1],
    type=type,
    alias=doc[tags == '\\alias'] %>% unlist,
    package=pkg_name,
    usages=usages,
    title=doc[tags == '\\title'] %>% unlist %>% paste(collapse='') %>% trimws(),
    description=doc[tags == '\\details'] %>% unlist %>% paste(collapse='') %>% trimws(),
    concepts=unique(concepts),
    keywords=c())
}

#' A simple dictionary to cache loaded terms
#'
#' @param the name of the term
#' @param env the environment name for the function; if `NULL`, look
#'   up in cache, otherwise insert or overwrite.
#'
#' @return A character string giving the name of the environment
#'   containing the function, or `NULL` if not in cache.
#' @noRd
ergmTermCache <- local({
  cache <- lapply(SUPPORTED_TERM_TYPES, function(x) list())
  pkglist <- character(0) # Current list of packages.

  # Reset the cache and update the list of watched packages.
  unload <- function(pkg_name) {
    pkglist <<-.packages(TRUE)

    for (term_type in SUPPORTED_TERM_TYPES) {
      for (term in names(cache[[term_type]])) {
        if (cache[[term_type]][[term]]$package == pkg_name) {
          cache[[term_type]][[term]] <<- NULL
        }
      }
    }
  }

  # Crawl all loaded packages for terms, parse them once and store in the singleton store
  load <- function(pkg_name) {
    pkg <- tools::Rd_db(pkg_name)
    all_doco <- attributes(pkg)$names
    converted <- all_doco[grep(SUPPORTED_TERM_TYPE_REGEX, all_doco)]

    for (term in lapply(converted, .parseTerm, pkg, pkg_name)) {
      cache[[term$type]][[term$link]] <<- term
    }
  }

  # Check if new namespaces have been added.
  checknew <- function() {
    loaded_packages <- .packages(TRUE)
    db <- utils::hsearch_db()$Base
    term_packages <- unique(db$Package[grep(SUPPORTED_TERM_TYPE_REGEX, db$Topic)])
    for (pkg_name in intersect(loaded_packages, term_packages)) {
      if (!pkg_name %in% pkglist) {
        load(pkg_name)

        setHook(packageEvent(pkg_name, "detach"), unload)
        setHook(packageEvent(pkg_name, "onUnload"), unload)
      }
    }

    cache <<- lapply(cache, function(terms) terms[sort(names(terms))])

    pkglist <<- loaded_packages
  }

  function (term_type) {
    checknew()

    cache[[term_type]]
  }
})

.buildTermsDataframe <- function(terms) {
  df <- c()
  for (term in terms) {
    if (!is.null(term$usages[[1]]$type)) {
      usage <- paste(sprintf('`%s` (%s)',
        sapply(term$usages, "[[", 'usage') %>% gsub('\\$', '\\\\$', .) %>% gsub('`', '', .) %>% gsub(' *=[^,)]*(,|\\)) *', '\\1 ', .) %>% trimws,
        sapply(term$usages, "[[", 'type')), collapse='\n')
    } else {
      usage <- paste(sprintf('`%s`', sapply(term$usages, '[[', 'usage')), collapse='\n')
    }

    df <- rbind(df, c(usage, term$package, term$title, if (length(term$concepts) > 0) paste(term$concepts, collapse='\n') else '', term$link))
  }

  df <- data.frame(df, stringsAsFactors=FALSE)
  colnames(df) <- c('Term', 'Package', 'Description', 'Concepts', 'Link')

  df
}

# terms : a list structure of the documentation data
# categores : an optional vector of column names to print and include
# only.include : an optional vector of categories, only terms that match the category will be printed 
.termMatrix <- function(terms, categories=NULL, only.include=NULL) {
  # if list of categories not supplied, generate it
  # otherwise, use the categories (and column order) provided
  if (is.null(categories)) {
    categories <- unique(unlist(sapply(terms,'[[','concepts')))
  }

  # figure out which terms should be included
  if (!is.null(only.include)) {
    included <- sapply(terms, function(term) any(term$concepts %in% only.include))
    terms <- terms[included]
  }

  # figure out which terms are members of each cat
  membership <- lapply(categories, function(cat) {
    # return true for terms that match cat
    sapply(terms, function(term) cat %in% term$concepts)
  })

  df <- data.frame(membership)

  categories <- unlist(CONCEPT_ABBREVIATIONS[categories])
  colnames(df) <- categories
  rownames(df) <- NULL
  df$Link <- sapply(terms,'[[','link')
  df$Term <- sapply(terms,'[[','name')

  df[c('Term', categories, 'Link')]
}

# output listings of terms, grouped by category
.termToc <- function(terms) {
  cats <- unique(unlist(sapply(terms, '[[', 'concepts')))

  ret <- list()
  for (cat in cats) {
    #find ones with matching cats
    matchingTerms <- terms[sapply(terms, function(term) any(term$concepts==cat))]
    ret[[cat]] <- list(
      link=sapply(matchingTerms, '[[', 'link'),
      name=sapply(matchingTerms, '[[', 'name')
    )
  }

  ret
}

# Generate the dynamic index text
.generateDynamicIndex <- function(term_type, formatter) {
  terms <- ergmTermCache(term_type)

  index <- .buildTermsDataframe(terms)

  if (term_type == 'ergmTerm') {
    m_freq <- .termMatrix(terms, categories=c('binary', 'valued', 'directed', 'undirected', 'bipartite', 'dyad-independent','operator','layer-aware'), only.include='frequently-used')
    m_op <- .termMatrix(terms, categories=c('binary', 'valued', 'directed', 'undirected', 'bipartite', 'dyad-independent', 'layer-aware'), only.include='operator')
  } else {
    m_freq <- NULL
    m_op <- NULL
  }

  m_all <- .termMatrix(terms)
  toc <- .termToc(terms)

  formatter(index, m_freq, m_op, m_all, toc)
}

.formatIndexText <- function(df) {
  df$Term <- gsub('`', '', gsub('valued', 'val', gsub('binary', 'bin', df$Term)))

  line_wrap <- function(lines, max_width) {
    lines <- unlist(strsplit(sapply(strsplit(lines, '\n'), stringr::str_wrap, max_width), '\n'))

    out <- c()
    for (line in lines) {
      while (nchar(line) > max_width) {
        out <- c(out, substr(line, 1, max_width))
        line <- substr(line, max_width + 1, nchar(line))
      }
      out <- c(out, line)
    }
    out
  }

  pad_lines <- function(lines, max_lines) {
    c(lines, rep('', max_lines - length(lines)))
  }

  colnames(df)[[2]] <- 'Pkg'
  out <- sprintf('|%s|\n', paste(stringr::str_pad(names(DISPLAY_TEXT_INDEX_MAX_WIDTHS), DISPLAY_TEXT_INDEX_MAX_WIDTHS, side='right', pad='-'), collapse='|'))
  empty_row <- sprintf('|%s|\n', paste(stringr::str_pad(rep('', length(DISPLAY_TEXT_INDEX_MAX_WIDTHS)), DISPLAY_TEXT_INDEX_MAX_WIDTHS), collapse='|'))

  r <- list()
  for (i in 1:dim(df)[1]) {
    for (c in names(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) {
      r[[c]] <- if (df[i, c] != '') line_wrap(df[i, c], DISPLAY_TEXT_INDEX_MAX_WIDTHS[[c]]) else c()
    }

    max_lines <- max(sapply(r, length))
    for (c in names(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) {
      r[[c]] <- pad_lines(r[[c]], max_lines)
    }

    for (j in 1:max_lines) {
      out <- sprintf('%s|%s|\n', out, paste(stringr::str_pad(sapply(r, "[[", j), DISPLAY_TEXT_INDEX_MAX_WIDTHS, side='right'), collapse='|'))
    }
    out <- paste(out, empty_row, sep='')
  }

  out
}

.formatMatrixText <- function(df) {
  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  paste(knitr::kable(df, 'pipe'), collapse='\n')
}

.formatTocText <- function(toc) {
  out <- '**** Term index by category ****'
  for (cat in names(toc)) {
    out <- sprintf('%s\n\n%s:\n%s', out, cat,
      paste(paste('   ', stringr::str_wrap(paste(toc[[cat]]$name, collapse=', '), DISPLAY_TEXT_MAX_WIDTH)), collapse='\n'))
  }
  out
}

.formatText <- function(index, m_freq=NULL, m_op=NULL, m_all=NULL, toc=NULL) {
  sprintf('\\preformatted{%s\n%s\n\n%s\n\n%s\n\n%s}',
    .formatIndexText(index),
    if (!is.null(m_freq)) .formatMatrixText(m_freq) else '',
    if (!is.null(m_op)) .formatMatrixText(m_op) else '',
    if (!is.null(m_all)) .formatMatrixText(m_all) else '',
    if (!is.null(toc)) .formatTocText(toc) else '')
}

# Format the table for text output
.formatIndexLatex <- function(df) {
  df$Term <- gsub('valued', 'val', gsub('binary', 'bin', df$Term))
  df$Term <- gsub('\n', ' \\\\newline ', df$Term) %>% gsub('`([^`]*)`', '\\1', .) %>%
    strsplit(' ') %>% sapply(., function(x) paste(sprintf('\\code{%s}', x), collapse=' ')) %>%
    gsub('\\\\code\\{\\(([^(]*)\\)\\}', '(\\1)', .) %>% gsub('\\\\code\\{\\\\newline\\}', '\\\\newline', .) %>%
    paste('\\\\raggedright \\\\allowbreak', .)
  df$Link <- NULL
  knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1f\\textwidth}', DISPLAY_LATEX_INDEX_PCT_WIDTHS), vline="") %>%
    gsub(' *\n *', ' ', .) %>%
    gsub('\\\\ ', '\\\\\\\\ ', .)
}

.formatMatrixLatex <- function(df) {
  df$Link <- NULL

  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1fcm}', DISPLAY_LATEX_TOC_PCT_WIDTHS(dim(df)[2] - 1)), vline="") %>%
    gsub(' *\n *', ' ', .) %>%
    gsub('\\\\ ', '\\\\\\\\ ', .)
}

.formatTocLatex <- function(toc) {
  out <- '\\noindent\\textbf{Terms by concepts}\n\n\\begin{description}'
  for (cat in names(toc)) {
    out <- sprintf('%s\n\\item[%s] %s', out, cat, paste(toc[[cat]]$name, collapse=', '))
  }
  paste(out, "\n\\end{description}")
}

.formatLatex <- function(index, m_freq=NULL, m_op=NULL, m_all=NULL, toc=NULL) {
  sprintf('\\out{%s%s%s%s%s}',
    .formatIndexLatex(index),
    if (!is.null(m_freq)) .formatMatrixLatex(m_freq) else '',
    if (!is.null(m_op)) .formatMatrixLatex(m_op) else '',
    if (!is.null(m_all)) .formatMatrixLatex(m_all) else '',
    if (!is.null(toc)) .formatTocLatex(toc) else '')
}

# Format the table for text output
.formatIndexHtml <- function(df) {
  df$Term <- gsub('valued', 'val', gsub('binary', 'bin', df$Term))

  # Hack! because HTML code has to be wrapped \out{} which prevents \link{} from being parsed, the link has
  # to be constructed with the actual address. To find the address, generate the correct link with 
  # \link[=absdiff-ergmTerm]{test} and check that it works with \out{<a href="../help/absdiff-ergmTerm">test</a>}.
  # This address may change from an upstream R-studio change

  df$Term <- sprintf(gsub('`([^`(]*)([^`]*)`', '<span class="code"><a href="../help/%1$s" id="%1$s">\\1\\2</a></span>', gsub('\n', '<br />', df$Term)), df$Link)
  df$Link <- NULL

  knitr::kable(df, 'html', escape=FALSE, table.attr='class="striped"')
}

.formatMatrixHtml <- function(df) {
  df$Term <- sprintf('<a href="#%s">%s</a>', df$Link, df$Term)
  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], '&#10004;', '')
  }

  knitr::kable(df, 'html', escape=FALSE, table.attr='class="matrix"')
}

.formatTocHtml <- function(toc) {
  out <- paste('Jump to category:', paste(sprintf('<a href="#cat_%s">%s</a>', names(toc), names(toc)), collapse=' '))
  for (cat in names(toc)) {
    out <- sprintf('%s<h3><a id="%s">%s</a></h3>%s', out, cat, cat,
      paste(sprintf('<a href="#%s">%s</a>', toc[[cat]]$link, toc[[cat]]$name), collapse=' '))
  }
  out
}

.formatHtml <- function(index, m_freq=NULL, m_op=NULL, m_all=NULL, toc=NULL) {
  css <- '<style>.striped th,.striped td {padding:3px 10px} .striped tbody tr:nth-child(odd) {background: #eee} .striped .code {font-family: monospace; font-size:75\\%} .matrix td {align: center} .matrix th,.matrix td {padding-right:5px; width: 75px}</style>'
  sprintf('\\out{%s%s%s%s%s%s}', css,
    .formatIndexHtml(index),
    if (!is.null(m_freq)) .formatMatrixHtml(m_freq) else '',
    if (!is.null(m_op)) .formatMatrixHtml(m_op) else '',
    if (!is.null(m_all)) .formatMatrixHtml(m_all) else '',
    if (!is.null(toc)) .formatTocHtml(toc) else '')
}

# function to look up the set of terms applicable for a specific network

#' Search the ergm-terms documentation for appropriate terms
#' 
#' Searches through the \code{\link{ergm.terms}} help page and prints out a
#' list of terms appropriate for the specified network's structural
#' constraints, optionally restricting by additional categories and keyword
#' matches.
#' 
#' Uses \code{\link{grep}} internally to match keywords against the term
#' description, so \code{keywords} is currently matched as a single phrase.
#' Category tags will only return a match if all of the specified tags are
#' included in the term.
#' 
#' @param keyword optional character keyword to search for in the text of the
#' term descriptions. Only matching terms will be returned. Matching is case
#' insensitive.
#' @param net a network object that the term would be applied to, used as
#' template to determine directedness, bipartite, etc
#' @param categories optional character vector of category tags to use to
#' restrict the results (i.e. 'curved', 'triad-related')
#' @param name optional character name of a specific term to return
#' @return prints out the name and short description of matching terms, and
#' invisibly returns them as a list.  If \code{name} is specified, prints out
#' the full definition for the named term.
#' @author skyebend@uw.edu
#' @seealso See also \code{\link{ergm.terms}} for the complete documentation
#' @examples
#' 
#' # find all of the terms that mention triangles
#' search.ergmTerms('triangle')
#' 
#' # two ways to search for bipartite terms:
#' 
#' # search using a bipartite net as a template
#' myNet<-network.initialize(5,bipartite=3)
#' search.ergmTerms(net=myNet)
#' 
#' # or request the bipartite category
#' search.ergmTerms(categories='bipartite')
#' 
#' # search on multiple categories
#' search.ergmTerms(categories=c('bipartite','dyad-independent'))
#' 
#' # print out the content for a specific term
#' search.ergmTerms(name='b2factor')
#' 
#' @importFrom utils capture.output
#' @export search.ergmTerms
search.ergmTerms<-function(keyword,net,categories,name){
  
  if (!missing(net)){
    if(!is.network(net)){
      stop("the 'net' argument must be the network argument that applicable terms are to be searched for")
    }
  }
  
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

  terms <- ergmTermCache('ergmTerm')
  
  found<-rep(TRUE,length(terms))
  
  # if name is specified, restrict to terms with that name
  if(!missing(name)){
    for (t in seq_along(terms)){
      term<-terms[[t]]
      found[t]<-any(term$name==name || term$link==name)
    }
  }
  
  # restrict by categories
  for (t in which(found)){
    term<-terms[[t]]
    if(!all(categories%in%term$concepts)){
      found[t]<-FALSE }
  }
  
  # next (optionally) restrict by keyword matches
  if (!missing(keyword)){
    for (t in which(found)){
      term<-terms[[t]]
      # if we don't find the keyword in the text grep, mark it as false
      if(length(grep(keyword,c(term$description, term$title),ignore.case=TRUE))==0){
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
        output<-c(output, term)
        cat(sprintf('%s\n    %s: %s\n    Categories: %s\n\n',
          term$usages[[1]]$usage,
          term$title,
          term$description,
          paste(term$concepts, collapse=', ')))
      }
    }
  }else{
    for (t in which(found)){
      term<-terms[[t]]
      for (usage in term$usages) {
        outText <- sprintf('%s\n    %s\n', usage$usage, term$title)
        output<-c(output,outText)
      }
    }
    cat("Found ",length(output)," matching ergm terms:\n")
    cat(paste(output,collapse='\n'))
  }
  invisible(output)
}

#' An index of Ergm terms
#'
#' @name ergmTerm
#' @docType package
#' @description An index of Ergm terms
#'
#' @section Term index:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmTerm", ergm:::.formatLatex)}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmTerm", ergm:::.formatText)}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmTerm", ergm:::.formatHtml)}}
NULL

#' An index of Ergm constraints
#'
#' @name ergmConstraint
#' @docType package
#' @description An index of Ergm constraints
#'
#' @section Constraint index:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmConstraint", ergm:::.formatLatex)}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmConstraint", ergm:::.formatText)}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmConstraint", ergm:::.formatHtml)}}
NULL

#' An index of Ergm reference
#'
#' @name ergmReference
#' @docType package
#' @description An index of Ergm references
#'
#' @section Reference index:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmReference", ergm:::.formatLatex)}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmReference", ergm:::.formatText)}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex("ergmReference", ergm:::.formatHtml)}}
NULL
