#  File ergm/R/ergm.getterms.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.getterms<-function(formula) {
    if ((dc<-data.class(formula)) != "formula")
        stop (paste("Invalid formula of class ",dc))
    trms<-terms(formula)
    if (trms[[1]]!="~")
        stop ("Formula must be of form 'network ~ model'.")
    trms
}
