#  File R/ergm-terms-index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons
#######################################################################

library(knitr)
library(magrittr)
library(stringr)

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
		usages[[length(usages) + 1]] <- list('type'=usage_line[2], 'usage'=usage_line[3])
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
		db <- hsearch_db()
		for (pkg_name in loaded_packages) {
			if (!pkg_name %in% pkglist) {
				if (!any(endsWith(db$Base$Topic[db$Base$Package == pkg_name], "-ergmTerm"))) next
				load(pkg_name)

				setHook(packageEvent(pkg_name, "detach"), unload)
				setHook(packageEvent(pkg_name, "onUnload"), unload)
			}
		}

		cache <<- cache[sort(names(cache))]

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
		usage <- paste(sprintf('`%s` (%s)',
			sapply(term$usages, "[[", 'usage') %>% gsub('\\$', '\\\\$', .) %>% gsub('`', '', .) %>% gsub(' *=[^,)]*(,|\\)) *', '\\1 ', .) %>% trimws,
			sapply(term$usages, "[[", 'type')), collapse='\n')
		df <- rbind(df, c(usage, term$package, term$title, paste(term$concepts, collapse='\n')))
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
	df$Term <- gsub('`', '', df$Term)
	for (c in colnames(df)) {
		df[[c]] <- as.character(df[[c]])
	}

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

	max_widths <- list('Term'=25, 'Pkg'=5, 'Description'=35, 'Concepts'=10)
	colnames(df)[[2]] <- 'Pkg'
	out <- sprintf('|%s|\n', paste(stringr::str_pad(colnames(df), max_widths, side='right', pad='-'), collapse='|'))
	empty_row <- sprintf('|%s|\n', paste(stringr::str_pad(rep('', length(max_widths)), max_widths), collapse='|'))

	r <- list()
	for (i in 1:dim(df)[1]) {
		for (c in colnames(df)) {
			r[[c]] <- line_wrap(df[i, c], max_widths[[c]])
		}

		max_lines <- max(sapply(r, length))
		for (c in colnames(df)) {
			r[[c]] <- pad_lines(r[[c]], max_lines)
		}

		for (j in 1:max_lines) {
			out <- sprintf('%s|%s|\n', out, paste(stringr::str_pad(sapply(r, "[[", j), max_widths, side='right'), collapse='|'))
		}
		out <- paste(out, empty_row, sep='')
	}

	sprintf('\\preformatted{%s}', out)
}

# Format the table for text output
.formatLatex <- function(df) {
	df$Term <- gsub('\n', '\\\\newline ', df$Term) %>% gsub('`', '', .) %>%
		strsplit(' ') %>% sapply(., function(x) paste(sprintf('\\code{%s}', x), collapse=' '))
	latex <- knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1f\\textwidth}', c(0.35, 0.05, 0.5, 0.1)), vline="") %>%
		gsub(' *\n *', ' ', .) %>%
		gsub('\\\\', '\\\\\\\\', .)
	sprintf('\\out{%s}', latex)
}

# Format the table for text output
.formatHtml <- function(df) {
	# Hack! because HTML code has to be wrapped \out{} which prevents \link{} from being parsed, the link has
	# to be constructed with the actual address. To find the address, generate the correct link with 
	# \link[=absdiff-ergmTerm]{test} and check that it works with \out{<a href="../help/absdiff-ergmTerm">test</a>}.
	# This address may change from an upstream R-studio change

	df$Term <- gsub('`([^`(]*)([^`]*)`', '<a href="../help/\\1-ergmTerm">\\1\\2</a>', gsub('\n', '<br />', df$Term))

	css <- '<style>.striped th,.striped td {padding:3px 10px} .striped tbody tr:nth-child(odd) {background: #eee} .striped tr td:nth-child(1) {font-family: monospace; font-size:75\\%}</style>'
	sprintf('\\out{%s%s}', css, knitr::kable(df, 'html', escape=FALSE, table.attr='class="striped"'))
}

#' An index of Ergm terms
#'
#' @name ergmTerm
#' @docType package
#' @section Term index:
#'
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatLatex)}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatText)}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex(ergm:::.formatHtml)}}
NULL
