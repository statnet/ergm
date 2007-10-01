ergm.getterms<-function(formula) {
    if ((dc<-data.class(formula)) != "formula")
        stop (paste("Invalid formula of class ",dc))
    trms<-terms(formula)
    if (trms[[1]]!="~")
        stop ("Formula must be of form 'network ~ model'.")
    trms
}
