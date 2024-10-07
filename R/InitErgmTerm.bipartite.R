#  File R/InitErgmTerm.bipartite.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#  See InitErgmTerm.R for a general explanation 
#  of InitErgmTerm functions

# NOTE: a number of undocumented terms from this file have been removed 
# but the terms are retained on the experimental_terms svn branch

###################################### InitErgmTerm TERMS:  A
##########################################################



#########################################################

#' @templateVar name b1nodematch
#' @title Nodal attribute-based homophily effect for the first mode in a bipartite network
#' @description This term is introduced
#'   in Bomiriya et al (2014). With the default `alpha` and `beta` values, this term will
#'   simply be a homophily based two-star statistic. This term adds one statistic to the model
#'   unless `diff` is set to `TRUE` , in which case the term adds multiple network
#'   statistics to the model, one for each of (a subset of) the unique values of the `attr`
#'   attribute.
#'   
#' @details If an `alpha`
#'   discount parameter is used, each of these statistics gives the sum of
#'   the number of common second-mode nodes raised to the power `alpha` for each pair of
#'   first-mode nodes with that attribute. If a `beta` discount parameter is used, each
#'   of these statistics gives half the sum of the number of two-paths with two first-mode nodes
#'   with that attribute as the two ends of the two path raised to the power `beta` for each
#'   edge in the network.
#'
#' @usage
#' # binary: b1nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1, byb2attr=NULL,
#' #                     levels=NULL)
#'
#' @template ergmTerm-attr
#' @param diff by default, one statistic will be added to the model. If `diff` is set to `TRUE`, one statistic will be added for each unique value of the `attr` attribute
#' @param keep deprecated
#' @param alpha,beta optional discount parameters both of which take values from `[0, 1]`, only one should be
#'   set at one time
#' @param byb2attr specifies a
#'   second mode categorical attribute. Setting this argument
#'   will separate the orginal statistics based on the values of the set second mode attribute---
#'   i.e. for example, if `diff` is `FALSE` , then the sum of all the statistics for
#'   each level of this second-mode attribute will be equal to the original `b1nodematch`
#'   statistic where `byb2attr` set to `NULL` .
#' @templateVar explain select a subset of `attr` values to include.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @template ergmTerm-keep-dep
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b1nodematch	<-	function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = TRUE,
                varnames 		= c("attrname", "diff", "keep", "beta", "alpha", "byb2attr"), 				
                vartypes 		= c("character", "logical", "numeric", "numeric", "numeric", "character"), 	 
                defaultvalues = list(NULL, FALSE, NULL, 1, 1, NULL), 										 
                required 		= c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
				dep.inform = list(FALSE, FALSE, "levels", FALSE, FALSE, FALSE))
	attrarg <- a$attrname
  }else{
    a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = TRUE,
                varnames 		= c("attr", "diff", "keep", "beta", "alpha", "byb2attr", "levels"), 				
                vartypes 		= c(ERGM_VATTR_SPEC, "logical", "numeric", "numeric", "numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC), 	 
                defaultvalues = list(NULL, FALSE, NULL, 1, 1, NULL, NULL), 										 
                required 		= c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
				dep.inform = list(FALSE, FALSE, "levels", FALSE, FALSE, FALSE, FALSE)) 								
    attrarg <- a$attr
  }
  ### Process the arguments
  if (!is.numeric(a$beta) || a$beta>1 || a$beta<0)
    ergm_Init_stop("beta argument to b1nodematch must be between 0 and 1 inclusive.")
  if (!is.numeric(a$alpha) || a$alpha>1 || a$alpha<0)
    ergm_Init_stop("alpha argument to b1nodematch must be between 0 and 1 inclusive.")
  
  nodecov <- ergm_get_vattr(attrarg, nw, bip="b1")
  attrname <- attr(nodecov, "name")
  
  b1.len	<-	get.network.attribute(nw, "bipartite")# gives # of b1 nodes
  u 		<-  ergm_attr_levels(a$levels, nodecov, nw, sort(unique(nodecov)))		  # gives unique attrnames of b1
  if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && !is.null(a$keep)) u <- u[a$keep]
  
  #   Recode to numeric
  nodecov   <- match(nodecov, u, nomatch = length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov == (length(u)+1)
  nodecov[dontmatch] 	<- length(u) + (1:sum(dontmatch))
  ui 					<- seq(along = u)
  
  
   if (!is.null(a$byb2attr)) {										  			
  	b2nodecov 	<-	ergm_get_vattr(a$byb2attr, nw, bip="b2")
  	v 	 		<-  sort(unique(b2nodecov))	# unique byb2attrnames of b2 
  	b2attrsize	<-	length(v)	# to get the levels of the byb2attr				
  																		    	
  	#   Recode to numeric														
  	b2nodecov   				<- match(b2nodecov, v, nomatch = length(v)+1)	
  	# All of the "nomatch" should be given unique IDs so they never match:		
  	b2dontmatch 				<- b2nodecov == (length(v)+1)					
  	b2nodecov[b2dontmatch] 		<- length(v) + (1:sum(b2dontmatch))				
  	vi							<- seq(along = v)         						
  	if (length(vi) < 2) {ergm_Init_stop("byb2attr should have at least two levels")}		
  	
  }	else {b2attrsize <- 0}	# to indicate that an b2attr does not exist	
  
  			
  ### Construct the list to return
  if (a$diff) {
  	if(!is.null(a$byb2attr)) {																	
    	coef.names 		<- unlist(lapply(paste("b1nodematch", paste(attrname, collapse="."), u, sep = "."), function(x){paste(x, v, sep = ".")}))																										   	
    	inputs 			<- c(b2attrsize, a$beta, a$alpha, ui, nodecov, vi, b2nodecov)    	
    }else{																						
    	coef.names 		<- paste("b1nodematch", paste(attrname,collapse="."), u, sep=".")		
    	inputs 			<- c(b2attrsize, a$beta, a$alpha, ui, nodecov)				
    }																							
   	 	
  } else {
  	
  	if(!is.null(a$byb2attr)) {																	
    	coef.names 		<- paste("b1nodematch", paste(attrname,collapse="."), v, sep=".")		
    	inputs 			<- c(b2attrsize, a$beta, a$alpha, nodecov, vi, b2nodecov)
    }else{																						
    	coef.names 		<- paste("b1nodematch", paste(attrname,collapse="."), sep=".")		
   		inputs 			<- c(b2attrsize, a$beta, a$alpha, nodecov)												
    }																							
 																					
  }							
  list(name			=  "b1nodematch",                       #name: required
       coef.names 	=  coef.names,                          #coef.names: required
       inputs 		=  inputs
       )
}



##########################################################

#' @templateVar name b2nodematch
#' @title Nodal attribute-based homophily effect for the second mode in a bipartite network
#' @description This term is introduced in Bomiriya et al (2014).
#'   With the default `alpha` and `beta` values, this term will
#'   simply be a homophily based two-star statistic. This term adds one statistic to the model
#'   unless `diff` is set to `TRUE` , in which case the term adds multiple network
#'   statistics to the model, one for each of (a subset of) the unique values of the `attr`
#'   attribute.
#'   
#' @details If an `alpha`
#'   discount parameter is used, each of these statistics gives the sum of
#'   the number of common first-mode nodes raised to the power `alpha` for each pair of
#'   second-mode nodes with that attribute. If a `beta` discount parameter is used, each
#'   of these statistics gives half the sum of the number of two-paths with two second-mode nodes
#'   with that attribute as the two ends of the two path raised to the power `beta` for each
#'   edge in the network.
#'
#' @usage
#' # binary: b2nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1, byb1attr=NULL,
#' #                     levels=NULL)
#'
#' @template ergmTerm-attr
#' @param diff by default, one statistic will be added to the model. If `diff` is set to `TRUE`, one statistic will be added for each unique value of the `attr` attribute
#' @param keep deprecated
#' @param alpha,beta optional discount parameters both of which take values from `[0, 1]`, only one should be
#'   set at one time
#' @param byb2attr specifies a
#'   second mode categorical attribute. Setting this argument
#'   will separate the orginal statistics based on the values of the set second mode attribute---
#'   i.e. for example, if `diff` is `FALSE` , then the sum of all the statistics for
#'   each level of this second-mode attribute will be equal to the original `b1nodematch`
#'   statistic where `byb2attr` set to `NULL` .
#' @templateVar explain select a subset of `attr` values to include.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-keep-dep
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b2nodematch	<-	function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                  varnames = c("attrname", "diff", "keep", "beta", "alpha", "byb1attr"),# RPB - 10/03/2012 - added the new arg "byb1attr"
                  vartypes = c("character", "logical", "numeric", "numeric", "numeric", "character"),
                  defaultvalues 	= list(NULL, FALSE, NULL, 1, 1, NULL),
                  required		= c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
				  dep.inform = list(FALSE, FALSE, "levels", FALSE, FALSE, FALSE))
	attrarg <- a$attrname  
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                  varnames = c("attr", "diff", "keep", "beta", "alpha", "byb1attr", "levels"),# RPB - 10/03/2012 - added the new arg "byb1attr"
                  vartypes = c(ERGM_VATTR_SPEC, "logical", "numeric", "numeric", "numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                  defaultvalues 	= list(NULL, FALSE, NULL, 1, 1, NULL, NULL),
                  required		= c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
				  dep.inform = list(FALSE, FALSE, "levels", FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
  }
   ### Process the arguments
  if (!is.numeric(a$beta) || a$beta>1 || a$beta<0)
    ergm_Init_stop("beta argument to b2nodematch must be between 0 and 1 inclusive.")
  if (!is.numeric(a$alpha) || a$alpha>1 || a$alpha<0)
    ergm_Init_stop("alpha argument to b2nodematch must be between 0 and 1 inclusive.")
  
  nodecov <- ergm_get_vattr(attrarg, nw, bip="b2")
  attrname <- attr(nodecov, "name")
  
  b1.len 	<-	get.network.attribute(nw, "bipartite")# gives # of b1 nodes
  u 	 	<-  ergm_attr_levels(a$levels, nodecov, nw, sort(unique(nodecov)))	  # gives unique attrnames of b2's
   																  # gives unique byb1attrnames of b1's
  if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && !is.null(a$keep)) u <- u[a$keep]
																  
  #   Recode to numeric
  nodecov 				<- match(nodecov, u, nomatch = length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch 			<- nodecov ==(length(u)+1)
  nodecov[dontmatch] 	<- length(u) + (1:sum(dontmatch))
  ui 					<- seq(along = u)
  
  if (!is.null(a$byb1attr)) {										# gives unique byb1attrnames of b1's
  	
  	b1nodecov 	<-	ergm_get_vattr(a$byb1attr, nw, bip="b1")												# get byb1attr vals

  	v 	 		<-  sort(unique(b1nodecov))# gives unique byb1attrnames of b1's
  	b1attrsize	<-	length(v)	# to get the levels of the byb1attr			
  	
  	#   Recode to numeric														
  	b1nodecov   				<- match(b1nodecov, v, nomatch = length(v)+1)	
  	# All of the "nomatch" should be given unique IDs so they never match:		
  	b1dontmatch 				<- b1nodecov == (length(v)+1)					
  	b1nodecov[b1dontmatch] 		<- length(v) + (1:sum(b1dontmatch))				
  	vi							<- seq(along = v)         						
  	
  	if (length(vi) < 2){ergm_Init_stop("byb1attr should have at least two levels")}		
  
  } else {b1attrsize <- 0}	# to indicate that an b1attr does not exist 	
   
  # if(drop) { # Check for zero statistics, print -Inf messages if applicable	# removed the whole if statement
    # obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    # ew <- extremewarnings(obsstats)
    # u <- u[!ew]
    # ui <- ui[!ew]
  # }
  ### Construct the list to return
  if (a$diff) {
  	if(!is.null(a$byb1attr)) {																	
    	coef.names 		<- unlist(lapply(paste("b2nodematch", paste(attrname, collapse="."), u, sep = "."), function(x){paste(x, v, sep = ".")}))																		
    	inputs 			<- c(b1attrsize, a$beta, a$alpha, ui, nodecov, vi, b1nodecov)
    }else{																						
    	coef.names 		<- paste("b2nodematch", paste(attrname,collapse="."), u, sep=".")		
    	inputs 			<- c(b1attrsize, a$beta, a$alpha, ui, nodecov)				
    }																							
   	
 	
  } else {
  	
  	if(!is.null(a$byb1attr)) {																	
    	coef.names 		<- paste("b2nodematch", paste(attrname,collapse="."), v, sep=".")		
    	inputs 			<- c(b1attrsize, a$beta, a$alpha, nodecov, vi, b1nodecov)	
    }else{																						
    	coef.names 		<- paste("b2nodematch", paste(attrname,collapse="."), sep=".")		
    	inputs 			<- c(b1attrsize, a$beta, a$alpha, nodecov)					
    }																							
    																						
  }	
  list(name="b2nodematch",                               #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs
       )
}



###################################### InitErgmTerm TERMS:  B
#########################################################



###################################### InitErgmTerm TERMS:  E
#########################################################


###################################### InitErgmTerm TERMS:  G
#########################################################


###################################### InitErgmTerm TERMS:  M
#########################################################



