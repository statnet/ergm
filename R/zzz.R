#==============================================================================
# This file contains the following 2 conventional functions used by R to load
# and unload packages:
#             <.First.lib>
#             <.Last.lib>
#==============================================================================



###############################################################################
# The <.First.lib> function loads the compiled ergm code and prints the
# copyright information; <.First.lib> is called when ergm is loaded by library()
#
# --PARAMETERS--
#   lib: the name of the library directory where 'pkg' is stored
#   pkg: the name of the package
#
# --RETURNED--
#   a libraryIQR object
#
###############################################################################

.First.lib <- function(lib, pkg){
  library.dynam("ergm", pkg, lib)
  DESCpath <- file.path(system.file(package="ergm"), "DESCRIPTION")
  info <- read.dcf(DESCpath)
  cat('\nergm:', info[,"Title"], 
      '\nVersion', info[,"Version"], 'created on', info[,"Date"], '\n') 
  cat(paste("copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
"                    David R. Hunter, Penn State University\n",
"                    Carter T. Butts, University of California-Irvine\n",
"                    Steven M. Goodreau, University of Washington\n",
"                    Pavel N. Krivitsky, Carnegie Mellon University\n",
"                    Martina Morris, University of Washington\n",sep=""))
  cat('Type help(package="ergm") to get started.\n\n')
  cat('Based on "statnet" project software (http://statnetproject.org).\n',
      'For license and citation information see http://statnetproject.org/attribution\n',
      'or type citation("ergm").\n')
#   cat('Please cite it when you use it!\n')
#   cat('To cite, see citation("ergm")\n')
#   require(network, quietly=TRUE)
}



#############################################################
# The <.Last.lib> function unloads the compiled ergm code
#
# --PARAMETERS--
#   libpath: the complete path to the package, as a string
#############################################################

.Last.lib <- function(libpath){
  library.dynam.unload("ergm",libpath)
}
