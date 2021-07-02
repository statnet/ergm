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

FREQUENTLY_USED_TERM_CATEGORIES <- c('binary', 'valued', 'directed', 'undirected', 'bipartite', 'dyad-independent','operator','layer-aware')
OPERATOR_CATEGORIES <- c('binary', 'valued', 'directed', 'undirected', 'bipartite', 'dyad-independent', 'layer-aware')

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
    concepts=if (!is.null(concepts)) unique(concepts) else c(),
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
  names(cache) <- paste('ergm', SUPPORTED_TERM_TYPES, sep='')
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

.buildTermsDataframe <- function(term_type) {
  terms <- ergmTermCache(term_type)

  if (length(terms) == 0) return(NULL)

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
.termMatrix <- function(term_type, categories=NULL, only.include=NULL) {
  terms <- ergmTermCache(term_type)

  if (length(terms) == 0) return(NULL)

  # if list of categories not supplied, generate it
  # otherwise, use the categories (and column order) provided
  if (is.null(categories)) {
    categories <- unique(unlist(sapply(terms,'[[','concepts')))
  }

  if(length(categories)==0) return(NULL)

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
.termToc <- function(term_type) {
  terms <- ergmTermCache(term_type)

  if (length(terms) == 0) return(NULL)

  cats <- unique(unlist(sapply(terms, '[[', 'concepts')))
  if(length(cats)==0) return(NULL)

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

.formatIndexText <- function(df) {
  if(is.null(df)) return(NULL)

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

  sprintf('\\preformatted{%s}', out)
}

.formatMatrixText <- function(df) {
  if(is.null(df)) return(NULL)

  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  sprintf('\\preformatted{%s}', paste(knitr::kable(df, 'pipe'), collapse='\n'))
}

.formatTocText <- function(toc) {
  if(is.null(toc)) return(NULL)

  out <- ''
  for (cat in names(toc)) {
    out <- sprintf('%s\n\n%s:\n%s', out, cat,
      paste(paste('   ', stringr::str_wrap(paste(toc[[cat]]$name, collapse=', '), DISPLAY_TEXT_MAX_WIDTH)), collapse='\n'))
  }
  sprintf('\\preformatted{%s}', out)
}

# Format the table for text output
.formatIndexLatex <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term <- gsub('valued', 'val', gsub('binary', 'bin', df$Term))
  df$Term <- gsub('\n', ' \\\\newline ', df$Term) %>% gsub('`([^`]*)`', '\\1', .) %>%
    strsplit(' ') %>% sapply(., function(x) paste(sprintf('\\code{%s}', x), collapse=' ')) %>%
    gsub('\\\\code\\{\\(([^(]*)\\)\\}', '(\\1)', .) %>% gsub('\\\\code\\{\\\\newline\\}', '\\\\newline', .) %>%
    paste('\\\\raggedright \\\\allowbreak', .)
  df$Link <- NULL
  sprintf('\\out{%s}',
    knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1f\\textwidth}', DISPLAY_LATEX_INDEX_PCT_WIDTHS), vline="") %>%
      gsub(' *\n *', ' ', .) %>%
      gsub('\\\\ ', '\\\\\\\\ ', .))
}

.formatMatrixLatex <- function(df) {
  if(is.null(df)) return(NULL)

  df$Link <- NULL

  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  sprintf('\\out{%s}',
    knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1fcm}', DISPLAY_LATEX_TOC_PCT_WIDTHS(dim(df)[2] - 1)), vline="") %>%
      gsub(' *\n *', ' ', .) %>%
      gsub('\\\\ ', '\\\\\\\\ ', .))
}

.formatTocLatex <- function(toc) {
  if(is.null(toc)) return(NULL)

  out <- '\\out{\\noindent\\textbf{Terms by concepts}\n\n\\begin{description}'
  for (cat in names(toc)) {
    out <- sprintf('%s\n\\item[%s] %s', out, cat, paste(toc[[cat]]$name, collapse=', '))
  }
  paste(out, "\n\\end{description}}")
}

# Format the table for text output
.formatIndexHtml <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term <- gsub('valued', 'val', gsub('binary', 'bin', df$Term))

  # Hack! because HTML code has to be wrapped \out{} which prevents \link{} from being parsed, the link has
  # to be constructed with the actual address. To find the address, generate the correct link with 
  # \link[=absdiff-ergmTerm]{test} and check that it works with \out{<a href="../help/absdiff-ergmTerm">test</a>}.
  # This address may change from an upstream R-studio change

  df$Term <- sprintf(gsub('`([^`(]*)([^`]*)`', '<span class="code"><a href="../help/%1$s" id="%1$s">\\1\\2</a></span>', gsub('\n', '<br />', df$Term)), df$Link)
  df$Link <- NULL

  css <- '<style>.striped th,.striped td {padding:3px 10px} .striped tbody tr:nth-child(odd) {background: #eee} .striped .code {font-family: monospace; font-size:75\\%} .matrix td {align: center} .matrix th,.matrix td {padding-right:5px; width: 75px}</style>'
  sprintf('\\out{%s%s}', css, knitr::kable(df, 'html', escape=FALSE, table.attr='class="striped"'))
}

.formatMatrixHtml <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term <- sprintf('<a href="#%s">%s</a>', df$Link, df$Term)
  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], '&#10004;', '')
  }

  sprintf('\\out{%s}', knitr::kable(df, 'html', escape=FALSE, table.attr='class="matrix"'))
}

.formatTocHtml <- function(toc) {
  if(is.null(toc)) return(NULL)

  out <- paste('Jump to category:', paste(sprintf('<a href="#cat_%s">%s</a>', names(toc), names(toc)), collapse=' '))
  for (cat in names(toc)) {
    out <- sprintf('%s<h3><a id="%s">%s</a></h3>%s', out, cat, cat,
      paste(sprintf('<a href="#%s">%s</a>', toc[[cat]]$link, toc[[cat]]$name), collapse=' '))
  }
  sprintf('\\out{%s}', out)
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

#' Terms used in Exponential Family Random Graph Models
#'
#' @name ergmTerm
#' @aliases ergm-terms ergm.terms terms-ergm terms.ergm InitErgmTerm InitErgmWtTerm
#' @docType package
#' @description The function [`ergm`] is used to fit exponential random graph
#' models, in which the probability of a given network, \eqn{y}, on a set of
#' nodes is \deqn{h(y) \exp\{\eta(\theta) \cdot }{h(y) exp{eta(theta).g(y)} /
#' c(theta),}\deqn{ g(y)\}/c(\theta)}{h(y) exp{eta(theta).g(y)} / c(theta),}
#' where \eqn{h(y)} is the reference measure (for valued network models),
#' \eqn{g(y)} is a vector of network statistics for \eqn{y}, \eqn{\eta(\theta)}
#' is a natural parameter vector of the same length (with
#' \eqn{\eta(\theta)\equiv\theta}{eta(theta)=theta} for most terms),
#' \eqn{\cdot}{.} is the dot product, and \eqn{c(\theta)} is the normalizing
#' constant for the distribution.
#' 
#' The network statistics \eqn{g(y)} are entered as terms in the function call
#' to [`ergm`].  This page describes the possible terms (and hence
#' network statistics) included in [`ergm`][ergm-package] package.
#' 
#' A cross-referenced HTML version of the term documentation is available via
#' `vignette('ergm-term-crossRef')` and terms can also be searched via
#' [`search.ergmTerms`].
#'
#' @section Specifying models:
#' Terms to [`ergm`] are specified by a formula to represent the
#' network and network statistics. This is done via a `formula`, that is,
#' an formula object, of the form `y ~ <term 1> + <term 2> ...`, where
#' `y` is a network object or a matrix that can be coerced to a network
#' object, and `<term 1>`, `<term 2>`, etc, are each terms chosen
#' from the list given below.  To create a network object in , use the
#' [`network`] function, then add nodal attributes to it
#' using the `%v%` operator if necessary.
#' 
#' ## Operator terms
#' Operator terms like `B` and `F` take
#' formulas with other `ergm` terms as their arguments and transform them
#' by modifying their inputs (e.g., the network they evaluate) and/or their
#' outputs.
#' 
#' By convention, their names are capitalized and CamelCased.
#' 
#' ## Interactions
#' For binary ERGMs, interactions between [`ergm`] terms can be
#' specified in a manner similar to [`lm`] and others, as using the
#' `:` and `*` operators. However, they must be interpreted
#' carefully, especially for dyad-dependent terms. (Interactions involving
#' curved terms are not supported at this time.)
#' 
#' Generally, if term `a` has \eqn{p_a}{p[a]} statistics and `b` has
#' \eqn{p_b}{p[b]}, `a:b` will add \eqn{p_a \times p_b}{p[a]*p[b]}
#' statistics to the model, corresponding to each element of
#' \eqn{g_a(y)}{g[a](y)} interacted with each element of \eqn{g_b(y)}{g[b](y)}.
#' 
#' The interaction is defined as follows. Dyad-independent terms can be
#' expressed in the general form \eqn{g(y;x)=\sum_{i,j} }{sum[i,j]
#' x[i,j]*y[i,j]}\eqn{ x_{i,j}y_{i,j}}{sum[i,j] x[i,j]*y[i,j]} for some edge
#' covariate matrix \eqn{x}, \deqn{g_{a:b}(y)=\sum_{i,j}
#' x_{a,i,j}x_{b,i,j}y_{i,j}.}{g[a:b](y) = \sum[i,j] x[a,i,j]*x[b,i,j]*y[i,j].}
#' In other words, rather than being a product of their sufficient statistics
#' (\eqn{g_{a}(y)g_{b}(y)}{g[a](y)*g[b](y)}), it is a dyadwise product of their
#' dyad-level effects.
#' 
#' This means that an interaction between two dyad-independent terms can be
#' interpreted the same way as it would be in the corresponding logistic
#' regression for each potential edge. However, for undirected networks in
#' particular, this may lead to somewhat counterintuitive results. For example,
#' given two nodal covariates `"a"` and `"b"` (whose values for node
#' \eqn{i} are denoted \eqn{a_i}{a[i]} and \eqn{b_i}{b[i]}, respectively),
#' `nodecov("a")` adds one statistic of the form \eqn{\sum_{i,j}
#' (a_{i}+a_{j}) y_{i,j}}{sum[i,j] (a[i]+a[j])*y[i,j]} and analogously for
#' `nodecov("b")`, so `nodecov("a"):nodecov("b")` produces
#' \deqn{\sum_{i,j} (a_{i}+a_{j}) (b_{i}+b_{j}) y_{i,j}.}{sum[i,j]
#' (a[i]+a[j])*(b[i]+b[j])*y[i,j].}
#' 
#' ## Binary and valued ERGM terms
#' [`ergm`][ergm-package] functions such as [`ergm`] and
#' [`simulate`][simulate.formula] (for ERGMs) may operate in two
#' modes: binary and weighted/valued, with the latter activated by passing a
#' non-NULL value as the `response` argument, giving the edge attribute
#' name to be modeled/simulated.
#' 
#' ### Generalizations of binary terms
#' Binary ERGM statistics cannot be
#' used directly in valued mode and vice versa. However, a substantial number
#' of binary ERGM statistics --- particularly the ones with dyadic indepenence
#' --- have simple generalizations to valued ERGMs, and have been adapted in
#' [`ergm`][ergm-package]. They have the same form as their binary
#' ERGM counterparts, with an additional argument: `form`, which, at this
#' time, has two possible values: `"sum"` (the default) and
#' `"nonzero"`. The former creates a statistic of the form \eqn{\sum_{i,j}
#' x_{i,j} y_{i,j}}{sum[i,j] x[i,j]*y[i,j]}, where \eqn{y_{i,j}}{y[i,j]} is the
#' value of dyad \eqn{(i,j)} and \eqn{x_{i,j}}{x[i,j]} is the term's covariate
#' associated with it. The latter computes the binary version, with the edge
#' considered to be present if its value is not 0.  Valued version of some
#' binary ERGM terms have an argument `threshold`, which sets the value
#' above which a dyad is conidered to have a tie. (Value less than or equal to
#' `threshold` is considered a nontie.)
#' 
#' The `B()` operator term documented below can be used to pass other
#' binary terms to valued models, and is more flexible, at the cost of being
#' somewhat slower.
#' 
#' ## Nodal attribute levels and indices
#' Terms taking a categorical nodal covariate also take the `levels`
#' argument.  (There are analogous `b1levels` and `b2levels`
#' arguments for some terms that apply to bipartite networks, and the
#' `levels2` argument for mixing terms.)  The `levels` argument can
#' be used to control the set and the ordering of attribute levels.
#' 
#' Terms that allow the selection of nodes do so with the `nodes`
#' argument, which is interpreted in the same way as the `levels`
#' argument, where the categories are the relevant nodal indices themselves.
#' 
#' Both `levels` and `nodes` use the new level selection UI. (See
#' \link[=nodal_attributes]{Specifying Vertex attributes and Levels} (\verb{?
#' nodal_attributes}) for details.)
#' 
#' ### Legacy arguments
#' 
#' The legacy `base` and `keep` arguments are deprecated as of
#' version 3.10, and replaced by the `levels` UI. The `levels`
#' argument provides consistent and flexible mechanisms for specifying which
#' attribute levels to exclude (previously handled by `base`) and include
#' (previously handled by `keep`).  If `levels` or `nodes`
#' argument is given, then `base` and `keep` arguments are ignored.
#' The legacy arguments will most likely be removed in a future version.
#' 
#' Note that this exact behavior is new in version 3.10, and it differs
#' slightly from older versions: previously if both `levels` and
#' `base`/`keep` were given, `levels` argument was applied first
#' and then applied the `base`/`keep` argument. Since version 3.10,
#' `base`/`keep` would be ignored, even if old term behavior is
#' invoked (as described in the next section).
#' 
#' ## Term versioning
#' When a term's behavior has changed from prior version, it is often possible
#' to invoke the old behavior by setting and/or passing a `version` term
#' option, giving the verison (constructed by [`as.package_version`])
#' desired.
#' 
#' ## Custom `ergm` terms
#' Users and other packages may build custom terms, and package
#' \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} provides
#' tools for implementing them.
#' 
#' The current recommendation for any package implementing additional terms is
#' to document the term with Roxygen comments and a name in the form
#' termName-ergmTerm. This ensures that \code{help("ergmTerm")} will list ERGM
#' terms available from all loaded packages.
#'
#' @section Terms included in the [`ergm`][ergm-package] package:
#' As noted above, a cross-referenced HTML version of the term documentation is
#' also available via `vignette('ergm-term-crossRef')` and terms
#' can also be searched via [`search.ergmTerms`].
#'
#' ## Term index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm"))}}
#'
#' ## Frequently-used terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#'
#' ## Operator terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' 
#' ## All terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm"))}}
#' 
#' ## Terms by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmTerm"))}}
#'
#' @seealso [`ergm`][ergm-package] package, [`search.ergmTerms`], [`ergm`], [`network`], [`%v%`], [`%n%`]
#'
#' @references 
#' - Krivitsky P. N., Hunter D. R., Morris M., Klumb
#' C. (2021). "ergm 4.0: New features and improvements."
#' arXiv:2106.04997. \url{https://arxiv.org/abs/2106.04997}
#' 
#' - Bomiriya, R. P, Bansal, S., and Hunter, D. R. (2014).  Modeling
#' Homophily in ERGMs for Bipartite Networks.  Submitted.
#' 
#' - Butts, CT.  (2008).  "A Relational Event Framework for Social
#' Action." *Sociological Methodology,* 38(1).
#' 
#' - Davis, J.A. and Leinhardt, S.  (1972).  The Structure of Positive
#' Interpersonal Relations in Small Groups.  In J. Berger (Ed.),
#' *Sociological Theories in Progress, Volume 2*, 218--251.  Boston:
#' Houghton Mifflin.
#' 
#' - Holland, P. W. and S. Leinhardt (1981). An exponential family of
#' probability distributions for directed graphs.  *Journal of the
#' American Statistical Association*, 76: 33--50.
#' 
#' - Hunter, D. R. and M. S. Handcock (2006). Inference in curved
#' exponential family models for networks. *Journal of Computational and
#' Graphical Statistics*, 15: 565--583.
#' 
#' - Hunter, D. R. (2007). Curved exponential family models for social
#' networks. *Social Networks*, 29: 216--230.
#' 
#' - Krackhardt, D. and Handcock, M. S. (2007).  Heider versus Simmel:
#' Emergent Features in Dynamic Structures. *Lecture Notes in Computer
#' Science*, 4503, 14--27.
#' 
#' - Krivitsky P. N. (2012). Exponential-Family Random Graph Models for
#' Valued Networks. *Electronic Journal of Statistics*, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#' 
#' - Robins, G; Pattison, P; and Wang, P.  (2009).  "Closure,
#' Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models
#' for Directed Social Networks." *Social Networks,* 31:105-117.
#' 
#' - Snijders T. A. B., G. G. van de Bunt, and C. E. G. Steglich.
#' Introduction to Stochastic Actor-Based Models for Network Dynamics.
#' *Social Networks*, 2010, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}
#' 
#' - Morris M, Handcock MS, and Hunter DR. Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 2008, 24(4), 1-24.
#' \url{https://www.jstatsoft.org/v24/i04}
#' 
#' - Snijders, T. A. B., P. E. Pattison, G. L. Robins, and M. S. Handcock
#' (2006). New specifications for exponential random graph models,
#' *Sociological Methodology*, 36(1): 99-153.
#' 
#' @keywords models
#' 
#' @examples
#' \dontrun{
#' ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#' 
#' ergm(molecule ~ edges + kstar(2:3) + triangle
#'                       + nodematch("atomic type",diff=TRUE)
#'                       + triangle + absdiff("atomic type"))
#' }
NULL
#TODO: Write a valued example.

#' Sample Space Constraints for Exponential-Family Random Graph Models
#'
#' @name ergmConstraint
#' @docType package
#' @description [`ergm`] is used to fit exponential-family random graph models
#' (ERGMs), in which the probability of a given network, \eqn{y}, on a set of
#' nodes is \eqn{h(y) \exp\{\eta(\theta) \cdot g(y)\}/c(\theta)}, where
#' \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}), \eqn{g(y)} is a
#' vector of network statistics for \eqn{y}, \eqn{\eta(\theta)} is a natural
#' parameter vector of the same length (with \eqn{\eta(\theta)=\theta} for most
#' terms), and \eqn{c(\theta)} is the normalizing constant for the
#' distribution.
#' 
#' This page describes the constraints (the networks \eqn{y} for which
#' \eqn{h(y)>0}) that are included with the [`ergm`][ergm-package]
#' package. Other packages may add new constraints.
#'
#' @section Constraints formula:
#' A constraints formula is a one- or two-sided formula whose left-hand side is
#' an optional direct selection of the `InitErgmProposal` function and
#' whose right-hand side is a series of one or more terms separated by
#' `"+"` and `"-"` operators, specifying the constraint.
#' 
#' The sample space (over and above the reference distribution) is determined
#' by iterating over the constraints terms from left to right, each term
#' updating it as follows: 
#' - If the constraint introduces complex
#' dependence structure (e.g., constrains degree or number of edges in the
#' network), then this constraint always restricts the sample space. It may
#' only have a `"+"` sign.
#' 
#' - If the constraint only restricts the set of dyads that may vary in the
#' sample space (e.g., block-diagonal structure or fixing specific dyads at
#' specific values) and has a `"+"` sign, the set of dyads that may
#' vary is restricted to those that may vary according to this constraint
#' *and* all the constraints to date.
#' 
#' - If the constraint only restricts the set of dyads that may vary in the
#' sample space but has a `"-"` sign, the set of dyads that may
#' vary is expanded to those that may vary according to this constraint
#' *or* all the constraints up to date.
#' 
#' For example, a constraints formula `~a-b+c-d` with all constraints
#' dyadic will allow dyads permitted by either `a` or `b` but only if they are
#' also permitted by `c`; as well as all dyads permitted by `d`. If `A`, `B`,
#' `C`, and `D` were logical matrices, the matrix of variable dyads would be
#' equal to `((A|B)&C)|D`.
#' 
#' Terms with a positive sign can be viewed as "adding" a constraint
#' while those with a negative sign can be viewed as "relaxing" a constraint.
#' 
#' @section Constraints implemented in the [`ergm`][ergm-package] package:
#' \describe{
#'    \item{\code{.} or \code{NULL} (dyad-independent)}{
#'      A placeholder for no constraints: all networks of a
#'      particular size and type have non-zero probability.
#'      Cannot be combined with other constraints.
#'    }
#' }
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#'
#' ## All constraints
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmConstraint"))}}
#' 
#' ## Constraints by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmConstraint"))}}
#'
#' @references
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \pkg{statnet} Tutorial. *Journal of Statistical Software*, 24(8).
#' \url{https://www.jstatsoft.org/v24/i08/}.
#' 
#' - Hunter, D. R. and Handcock, M. S. (2006) *Inference in curved
#' exponential family models for networks*, Journal of Computational and
#' Graphical Statistics.
#' 
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  *Journal of Statistical Software*, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#' 
#' - Karwa V, Krivitsky PN, and Slavkovi\'c AB (2016). Sharing Social Network
#' Data: Differentially Private Estimation of Exponential-Family Random Graph
#' Models. *Journal of the Royal Statistical Society, Series C*, 66(3):
#' 481-500. \doi{10.1111/rssc.12185}
#' 
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#' 
#' - Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 24(4). \url{https://www.jstatsoft.org/v24/i04/}.
#' @keywords models
NULL

#' Reference Measures for Exponential-Family Random Graph Models
#'
#' @name ergmReference
#' @aliases ergm-references references-ergm ergm.references references.ergm
#' @docType package
#' @description This page describes the possible reference measures (baseline distributions)
#' for found in the [`ergm`][ergm-package] package, particularly the
#' default (Bernoulli) reference measure for binary ERGMs.
#' 
#' The reference measure is specified on the RHS of a one-sided formula passed
#' as the `reference` argument to [`ergm`].  See the
#' [`ergm`] documentation for a complete description of how
#' reference measures are specified.
#'
#' @section Possible reference measures to represent baseline distributions:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmReference"))}}
#'
#' ## All references
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmReference"))}}
#' 
#' ## References by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmReference"))}}
#'
#' @seealso [`ergm`][ergm-package], [`network`], `sna`, [`summary.ergm`], [`print.ergm`], `\%v\%`, `\%n\%`
#' 
#' @references
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b). \pkg{ergm}:
#' A Package to Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks. *Journal of Statistical Software*, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#' 
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' @keywords models
NULL
