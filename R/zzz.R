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


