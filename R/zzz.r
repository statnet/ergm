######################################################################
#
# zzz.r
#
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# written December 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/statnet package
#
# .First.lib is run when the package is loaded with library(statnet)
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("statnet", pkg, lib)
    ehelp <- library(help="statnet",lib.loc=NULL,character.only=TRUE)$info[[1]]
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
               substring(ehelp[3],first=16),".\n", sep=""))
    cat(paste("copyright (c) 2003, Mark S. Handcock, University of Washington\n",
"                    David R. Hunter, Penn State University\n",
"                    Carter T. Butts, University of California-Irvine\n",
"                    Steven M. Goodreau, University of Washington\n",
"                    Martina Morris, University of Washington\n",sep=""))
    cat('Type help(package="statnet") to get started.\n')
    cat('To cite, see citation("statnet")\n')
    require(network, quietly=TRUE)
}
