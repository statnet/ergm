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

# "Declare" myLibLoc

myLibLoc <- NULL

.onAttach <- function(lib, pkg){
  ops <- options(warn = -1)
  on.exit(options(ops))
  library.dynam("ergm", pkg, lib)
  DESCpath <- file.path(system.file(package="ergm"), "DESCRIPTION")
  info <- read.dcf(DESCpath)
  packageStartupMessage(
    paste('\nergm: ', info[,"Title"], 
          '\nVersion ', info[,"Version"], ' created on ', info[,"Date"], '\n',
          "copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
          "                    David R. Hunter, Penn State University\n",
          "                    Carter T. Butts, University of California-Irvine\n",
          "                    Steven M. Goodreau, University of Washington\n",
          "                    Pavel N. Krivitsky, Penn State University\n",
          "                    Martina Morris, University of Washington\n",
          'Type help(package="ergm") to get started.\n\n',
          'Based on "statnet" project software (http://statnet.org).\n',
          'For license and citation information see http://statnet.org/attribution\n',
          'or type citation("ergm").\n', sep="")
 ) 
  # Remember where this package is located, to later make sure we load
  # the same version on a cluster node.
  unlockBinding("myLibLoc", environment(ergm.getCluster))
  assign("myLibLoc",lib,envir=environment(.onAttach))
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
