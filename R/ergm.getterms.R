#  File R/ergm.getterms.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
##########################################################################
# The <ergm.getterms> function returns the terms of a given formula and
# ensures that the formula is  indeed a formula with the necessary
# ~ operator
#
# --PARAMETERS--
#   formula: a formula
#
#
# --RETURNED--
#   trms: the terms object associated with the formula and returned by the 
#         native R function <terms>. see '?terms.object' for details about
#         the components of 'trms'.
#          
###########################################################################

ergm.getterms<-function(formula) {
    if ((dc<-data.class(formula)) != "formula")
        stop (paste("Invalid formula of class ",dc))
    trms<-terms(formula)
    if (trms[[1]]!="~")
        stop ("Formula must be of form 'network ~ model'.")
    trms
}
