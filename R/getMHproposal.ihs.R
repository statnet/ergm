MHproposals.ihs<-
  #         Class Constraint      Weights        MHP
  rbind(I(c("c", "",              "default",      "TNT")),
          c("c", "",              "TNT",          "TNT"),
          c("c", "",              "random",       "randomtoggle"),
          c("c", "",              "TNT10",        "TNT10"),
          c("c", "bd",            "default",       "TNT"),
          c("c", "bd",            "TNT",           "TNT"),
          c("c", "bd",            "random",       "randomtoggle"),
          c("c", "bd+edges",      "default",      "ConstantEdges"),
          c("c", "bd+edges",      "random",       "ConstantEdges"),          
          c("c", "",              "nonobserved",  "randomtoggleNonObserved"),
          c("c", "degrees",       "default",      "CondDegree"),
          c("c", "degrees",       "random",       "CondDegree"),

          c("c", "degreedist",    "default",      "CondDegreeDist"),
          c("c", "degreedist",    "random",       "CondDegreeDist"), 
          c("c", "indegreedist",  "default",      "CondInDegreeDist"),
          c("c", "indegreedist",  "random",       "CondInDegreeDist"), 

#          c("c", "indegrees",     "default",      "CondInDegree"),
#          c("c", "indegrees",     "random",       "CondInDegree"),
#          c("c", "outdegrees",    "default",      "CondOutDegree"),
#          c("c", "outdegrees",    "random",       "CondOutDegree"),
          c("c", "edges",         "default",      "ConstantEdges"),
          c("c", "edges",         "random",       "ConstantEdges"),
          c("c", "hamming",       "default",      "HammingTNT"),
          c("c", "hamming",       "random",       "HammingTNT"),
          c("c", "edges+hamming", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "random",       "HammingConstantEdges"),
          c("f", "",              "default",      "formationTNT"),
          c("f", "",              "TNT",          "formationTNT"),
          c("f", "",              "random",       "formation"),
          c("f", "bd",            "default",      "formationTNT"),
          c("f", "bd",            "TNT",          "formationTNT"),
          c("f", "bd",            "random",       "formation"),
          c("d", "",              "default",      "dissolution"),
          c("d", "",              "random",       "dissolution"),
          c("d", "bd",            "default",      "dissolution"),
          c("d", "bd",            "random",       "dissolution")
        )
MHproposals.ihs <- data.frame(I(MHproposals.ihs[,1]), I(MHproposals.ihs[,2]), 
                          I(MHproposals.ihs[,3]), I(MHproposals.ihs[,4]))  
colnames(MHproposals.ihs)<-c("Class","Constraints","Weights","MHP")
MHproposals <- MHproposals.ihs

MHproposal.ergm.ihs<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, model=NULL,weights=NULL,class="c"){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$prop.args
  if(is.null(nw)) nw<-object$network
  if(is.null(weights)) weights<-"default"
  if(is.null(model)){
    model<-if(class %in% c("c","f"))
      ergm.getmodel(object$formula,nw,...)
    else
      ergm.getmodel.dissolve(object$formula,nw,...)
  }  
  MHproposal(constraints,arguments=arguments,nw=nw,model=model,weights=weights,class=class)
}
