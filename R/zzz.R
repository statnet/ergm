#.onLoad <- function(lib, pkg) {
# FIXME: Setting binding for myLibLoc using old code below does not work.
# See also ergm.getCluster function in R/parallel.utils.R 
#  # "Declare" myLibLoc
#  myLibLoc <- NULL
#  # Remember where this package is located, to later make sure we load
#  # the same version on a cluster node.
#  unlockBinding("myLibLoc", environment(ergm.getCluster))
#  assign("myLibLoc",lib,pos=environment(ergm.getCluster))
#}

.onAttach <- function(lib, pkg){
  info <- packageDescription("ergm")
  packageStartupMessage(
    paste('\nergm: version ', info$Version, ', created on ', info$Date, '\n',
          "Copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
          "                    David R. Hunter, Penn State University\n",
          "                    Carter T. Butts, University of California-Irvine\n",
          "                    Steven M. Goodreau, University of Washington\n",
          "                    Pavel N. Krivitsky, Penn State University\n",
          "                    Martina Morris, University of Washington\n",
          'Based on "statnet" project software (statnet.org).\n',
          'For license and citation information see statnet.org/attribution\n',
          'or type citation("ergm").\n', sep="")
 )
}



#############################################################
# The <.Last.lib> function unloads the compiled ergm code
#
# --PARAMETERS--
#   libpath: the complete path to the package, as a string
#############################################################
#
#.Last.lib <- function(libpath){
#  library.dynam.unload("ergm",libpath)
#}

