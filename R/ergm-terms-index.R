#  File R/ergm-terms-index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons
#######################################################################

library(magrittr)
library(knitr)

# Return the index entry for a single term in the new format
.parseTerm <- function(name, pkg, pkg_name) {
    doc <- pkg[[name]]
    tags <- sapply(doc, attr, 'Rd_tag')

	raw_usage <- doc[tags == '\\usage'] %>% 
        unlist %>% 
        paste(collapse="") %>%
        gsub("\n *# *"," ", .) %>%
        strsplit("\n") %>%
        .[[1]] %>%
        trimws()
	usages <- list()
	for (usage_line in regmatches(raw_usage, regexec("^(binary|valued): *(.+)$", raw_usage))) {
		usages[[usage_line[2]]] <- usage_line[3]
	}

    ret <- list(
		name=substr(name, 1, nchar(name) - 3),
		alias=doc[tags == '\\alias'] %>% unlist,
		package=pkg_name,
		usages=usages,
		title=doc[tags == '\\title'] %>% unlist,
		description=doc[tags == '\\description'] %>% unlist %>% paste(collapse='') %>% trimws(),
		concepts=doc[tags == '\\concept'] %>% unlist,
		keywords=c())

    return(ret)
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
	cache <- list()
	pkglist <- character(0) # Current list of packages.

	# Reset the cache and update the list of watched packages.
	unload <- function(pkg_name) {
		pkglist <<-.packages(TRUE)

		for (term in names(cache)) {
			if (cache[[term]]$package == pkg_name) {
				cache[[term]] <<- NULL
			}
		}
	}

	# Crawl all loaded packages for terms, parse them once and store in the singleton store
	load <- function(pkg_name) {
		pkg <- tools::Rd_db(pkg_name)
		all_doco <- attributes(pkg)$names
		converted <- all_doco[which(endsWith(all_doco, '-ergmTerm.Rd'))]

		for (term in lapply(converted, .parseTerm, pkg, pkg_name)) {
			cache[[term$name]] <<- term
		}
	}

	# Check if new namespaces have been added.
	checknew <- function() {
		loaded_packages <- .packages(TRUE)
		for (pkg_name in loaded_packages) {
			if (!pkg_name %in% pkglist) {
				load(pkg_name)

				setHook(packageEvent(pkg_name, "detach"), unload)
				setHook(packageEvent(pkg_name, "onUnload"), unload)
			}
		}

		pkglist <<- loaded_packages
	}

	function (name=NULL) {
		checknew()

		if (is.null(name)) (cache) else (cache[[name]])
	}
})

.buildTermsDataframe <- function(terms) {
	df <- c()
	for (term in terms) {
		usage <- paste(sprintf('%s (%s)', term$usages, names(term$usages)), collapse=', ')
		df <- rbind(df, c(usage, term$package, term$title, paste(term$concepts, collapse=', ')))
	}

	df <- data.frame(df)
	colnames(df) <- c('Term', 'Package', 'Description', 'Concepts')

	return (df)
}

# Generate the dynamic index text
.generateDynamicIndex <- function(formatter) {
	df <- .buildTermsDataframe(ergmTermCache())
	formatter(df)
}

# Format the table for text output
.formatText <- function(df) {
	knitr::kable(df, 'simple')
}

# Format the table for text output
.formatLatex <- function(df) {
	knitr::kable(df, 'latex')
}

# Format the table for text output
.formatHtml <- function(df) {
	gsub(' *\n *', '', sprintf('\\out{%s}', knitr::kable(df, 'html')))
}

#' An index of Ergm terms
#'
#' @name ergmTerm
#' @docType package
#' @section Term index:
#'
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatHtml)}}
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatLatex)}}
#' \if{text}{\Sexpr[results=verbatim,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatText)}}
NULL
