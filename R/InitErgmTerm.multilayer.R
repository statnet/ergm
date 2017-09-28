## TODO: LL-Constrained proposals.
## TODO: Check that noncommutative LL operators work as intended.
.unenv <- function(f){
  environment(f) <- NULL
  f
}

.despace <- function(s) gsub("[[:space:]]", "", s)

.lspec_coef.names <- function(Llist){
  reprs <- sapply(seq_along(Llist), function(l){
    name <- names(Llist)[l]
    L <- Llist[[l]]
    s <- switch(class(L),
                formula = .despace(deparse(if(length(L)==2) L[[2]] else L)),
                character = L,
                as.character(L))
    if(NVL(name,"")!="") s <- paste0(name,"=",s)
    s
  })
  paste0("L(",paste0(reprs, collapse=","),")")
}

#' A multilayer network representation.
#'
#' A function for specifying the LHS of a multilayer (a.k.a. multiplex
#' a.k.a. multirelational a.k.a. multivariate) ERGM.
#'
#' @param ... layer specification, in one of three formats:
#' 
#'   1. An (optionally named) list of identically-dimensioned and
#'      directed networks.
#'
#'   1. Several networks as (optionally named) arguments.
#'
#'   1. A single network and a chararacter vector. Then, the layers are
#'      values of the named edge attributes.
#'
#' @return A network object with layer metadata.
#'
#' @seealso [Help on model specification][ergm-terms] for specific terms.
#' 
#' @examples
#'
#' data(florentine)
#'
#' # Method 1: list of networks
#' flo <- Layer(list(m = flomarriage, b = flobusiness))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 2: networks as arguments
#' flo <- Layer(m = flomarriage, b = flobusiness)
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 3: edge attributes:
#' flo <- flomarriage
#' flo <- flomarriage | flobusiness
#' flo[,, names.eval="m"] <- as.matrix(flomarriage)
#' flo[,, names.eval="b"] <- as.matrix(flobusiness)
#' flo <- Layer(flo, c("m","b"))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' @export
Layer <- function(...){
  args <- list(...)
  if(length(args)==2 && is(args[[1]], "network") && is(args[[2]], "vector")){
    nwl <-
      lapply(args[[2]], function(eattr){
        nw <- args[[1]]
        el <- as.edgelist(nw, eattr)
        nael <- as.edgelist(is.na(nw))
        nw[,]<-0
        nw[el[,-3]] <- el[,3]
        nw[nael] <- NA
        nw
      })
    names(nwl) <- args[[2]]
  }else if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multilayer specification. See help for information.")

  if(any(!is.na(suppressWarnings(sapply(names(nwl), as.numeric))))) warning("Using numeric layer names is ambiguous.")

  constraintsl <- lapply(nwl, get.network.attribute, "constraints")
  if(!all_identical(lapply(constraintsl, .unenv))) stop("Layers have differing constraint structures. This is not supported at this time.")
  constraints.obsl <- lapply(nwl, get.network.attribute, "constraints.obs")
  if(!all_identical(lapply(constraints.obsl, .unenv))) stop("Layers have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, blockID.vattr=".LayerID", blockName.vattr=".LayerName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "constraints.obs"))
  nw %n% "constraints" <-
      if(NVL(nwl[[1]]%n%"constraints",~.)==~.)
        ~blockdiag(".LayerID")
      else
        append.rhs.formula(nwl[[1]]%n%"constraints", alist(blockdiag(".LayerBlocks")), TRUE)
  nw %n% "constraints.obs" <- nwl[[1]]%n%"constraints.obs"

  nw
}


## InitErgmTerm..layer.nets <- function(nw, arglist, response=NULL, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c(),
##                       vartypes = c(),
##                       defaultvalues = list(),
##                       required = c())
##   list(name="_layer_nets", coef.names=c(), inputs=unlist(.block_vertexmap(nw, ".LayerID", TRUE)), dependence=FALSE)
## }

.depends_on_layers <- function(commands){
  coms <- commands[-1]
  if(any(coms==0)) coms <- coms[-(which(coms==0)+1)] # Drop all numeric literals (i.e., numbers preceded by 0).
  coms <- coms[coms>=1] # Drop all commands.
  unique(coms)
}

InitErgmTerm..layer.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("L"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  
  nwl <- .split_constr_network(nw,".LayerID",".LayerName")
  nwnames <- names(nwl)

  # Process layer specification
  namemap <- seq_along(nwl)
  names(namemap) <- nwnames
  
  ll <- pack.LayerLogic_formula_as_double(a$L, namemap)
  # Terms on this logical layer will induce dyadic independence if its
  # value depends on more than one other layer value.
  dependence <- length(.depends_on_layers(ll))>1
  
  if(test_eval.LayerLogic(ll, FALSE)) stop("Layer specifications that produce edges on the output layer for empty input layers are not supported at this time.", call.=FALSE)
  
  list(name="_layer_net", coef.names=c(), inputs=c(unlist(.block_vertexmap(nw, ".LayerID", TRUE)),ll), dependence=dependence)
}

pack.LayerLogic_formula_as_double <- function(formula, namemap){
  OPMAP <- list(
    # Unary operators
    c(`(` = NA,
      `!` = -1,
      `+` = NA,
      `-` = -16,
      `abs` = -17,
      `round` = -20,
      `sign` = -22),
    # Binary operators
    c(`&` = -2,
      `&&` = -2,
      `|` = -3,
      `||` = -3,
      `xor` = -4,
      `==` = -5,
      `!=` = -6,
      `<` = -7,
      `>` = -8,
      `<=` = -9,
      `>=` = -10,
      `+` = -11,
      `-` = -12,
      `*` = -13,
      `/` = -14,
      `%%` = -15,
      `^` = -18,
      `%/%` = -19,
      `round` = -21))

    
  lidMap <- function(l){
    switch(class(l),
           numeric = c(0,l),
           character =,
           name = if(regexpr('^[0-9]+$',l)!=-1) as.integer(as.character(l))
                  else namemap[as.character(l)])
  }
  
  postfix <- function(call, coml=c()){
    if(is.call(call)){
      op <- call[[1]]
      for(i in seq_along(call[-1])+1){
        coml <- c(coml, postfix(call[[i]]))
      }
      coml <- c(coml, OPMAP[[length(call)-1]][as.character(op)])
    }else{
      coml <- c(coml, lidMap(call))
    }
    coml[!is.na(coml)]
  }
  
  com <- postfix(formula[[length(formula)]])
  c(sum(com!=0), com)
}

test_eval.LayerLogic <- function(commands, lv){
  coms <- commands[-1]
  lv <- rep(lv, length.out=max(coms))
  stack <- c()
  if(sum(coms!=0)!=commands[1]) stop("Layer specification command vector specifies incorrect number of commands.", call.=FALSE)
  for(i in 1:commands[1]){
    com <- coms[1]
    if(com==0){
      coms <- coms[-1]
      com <- coms[1]
      stack <- c(com, stack)
    }else if(com==-1){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(!x0, stack)
    }else if(com==-2){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 && y0, stack)
    }else if(com==-3){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 || y0, stack)
    }else if(com==-4){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(xor(x0, y0), stack)
    }else if(com==-5){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 == y0, stack)
    }else if(com==-6){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 != y0, stack)
    }else if(com==-7){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 < y0, stack)
    }else if(com==-8){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 > y0, stack)
    }else if(com==-9){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 <= y0, stack)
    }else if(com==-10){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 >= y0, stack)
    }else if(com==-11){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 + y0, stack)
    }else if(com==-12){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 - y0, stack)
    }else if(com==-13){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 * y0, stack)
    }else if(com==-14){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 / y0, stack)
    }else if(com==-15){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %% y0, stack)
    }else if(com==-16){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(-x0, stack)
    }else if(com==-17){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(abs(x0), stack)
    }else if(com==-18){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 ^ y0, stack)
    }else if(com==-19){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %/% y0, stack)
    }else if(com==-20){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(round(x0), stack)
    }else if(com==-21){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(round(x0, y0), stack)
    }else if(com==-22){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(sign(x0), stack)
    }else{
      stack <- c(lv[com], stack)
    }
    coms <- coms[-1]
  }
  if(length(stack)!=1) stop("Invalid layer specification command sequence.", call.=FALSE)
  stack
}

.all_layers_terms <- function(n, LHS=NULL){
  if(is.null(LHS))
    lapply(seq_len(n), function(i) as.formula(substitute(~i,list(i=as.name(i)))))
  else
    lapply(seq_len(n), function(i) as.formula(substitute(lhs~i,list(lhs=LHS, i=as.name(i)))))
}

.mk_.layer.net_auxform <- function(ll, nl){
  trmcalls <- .layers_expand_dot(ll, nl)
  # Get the formula as a list of term calls.
  trmcalls <- lapply(trmcalls, function(ltrm) call(".layer.net", ltrm))
  append.rhs.formula(~.,trmcalls)[-2]
}

.layers_expand_dot <- function(ll, nl){
  if(is(ll, "formula")) ll <- list(ll)
  # Replace . with all layers.
  do.call(c, lapply(ll, function(f) if(f[[length(f)]]=='.') .all_layers_terms(nl, LHS = if(length(f)==3) f[[2]]) else list(as.formula(f))))
}

InitErgmTerm.L <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "Ls"),
                      vartypes = c("formula", "formula,list"),
                      defaultvalues = list(NULL, ~.),
                      required = c(TRUE, FALSE))
  f <- a$formula

  nwl <- .split_constr_network(nw,".LayerID",".LayerID")
  nwnames <- names(nwl)

  Ls <- a$Ls
  if(is(Ls, "formula")) Ls <- list(Ls)
  Ls.dotexp <- .layers_expand_dot(Ls, length(nwl))
  

  auxiliaries <- .mk_.layer.net_auxform(Ls, length(nwl))  
  nltrms <- length(term.list.formula(auxiliaries[[2]]))

  w <- rep(1,nltrms)
  have.LHS <- sapply(Ls.dotexp, length)==3
  w[have.LHS] <- as.numeric(sapply(lapply(Ls.dotexp[have.LHS], "[[", 2), eval,environment(Ls[[1]])))
  
  nw1 <- nwl[[1]]
  
  if(length(f)==2) f <- nonsimp.update.formula(f, nw1~.)
  else nw1 <- ergm.getnetwork(f)
  
  m <- ergm.getmodel(f, nw1, response=response,...)

  dependence <- !is.dyad.independent(m) || !is.dyad.independent(nonsimp.update.formula(auxiliaries, nw~., from.new="nw"))

  
  Clist <- ergm.Cprepare(nw1, m, response=response)
  inputs <- pack.Clist_as_num(Clist)
  
  inputs <- c(nltrms, w, inputs)

  gs <- ergm.emptynwstats.model(m) * nltrms
  
  c(list(name="OnLayer", coef.names = paste0(.lspec_coef.names(list(a$Ls)),":",m$coef.names), inputs=inputs, dependence=dependence, emptynwstats = gs, auxiliaries = auxiliaries),
    passthrough.curved.ergm.model(m, function(x) paste0(.lspec_coef.names(list(a$Ls)),":",x)))
}

InitErgmTerm.lCMB <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls"),
                      vartypes = c("formula,list"),
                      defaultvalues = list(~.),
                      required = c(FALSE))

  nwl <- .split_constr_network(nw,".LayerID",".LayerID")
  Ls <- a$Ls
  auxiliaries <- .mk_.layer.net_auxform(Ls, length(nwl))
  nltrms <- length(term.list.formula(auxiliaries[[2]]))

  inputs <- c(nltrms)

  list(name="layerCMB", coef.names = paste0('lCMB(',.despace(deparse(Ls)),')'), inputs=inputs, dependence=TRUE, auxiliaries = auxiliaries)
}

################################################################################
InitErgmTerm.ldegree<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "Ls", "dir"),
                      vartypes = c("numeric", "character", "formula,list", "character"),
                      defaultvalues = list(NULL, NULL, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree() has only one value", call.=FALSE)
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste("degree",d,sep="")
    name <- "degree"
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    # See comment in d_degree_w_homophily function
    coef.names <- paste("deg", d, ".homophily.",byarg, sep="")
    name <- "degree_w_homophily"
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]], sep="")
    name <- "degree_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }

  if(is.null(a$Ls) || is.null(a$dir)) stop("A layer and direction specification is required for the ldegree term.")
  
  c(.process_layers_degree(nw, a, name=name,coef.names=coef.names, inputs=inputs, emptynwstats=emptynwstats), minval = 0)
}

################################################################################
InitErgmTerm.gwldegree<-function(nw, arglist,  ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff", "Ls", "dir"),
                      vartypes = c("numeric", "logical", "character","numeric", "formula,list", "character"),
                      defaultvalues = list(NULL, FALSE, NULL, 30, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff

  if(is.null(a$Ls) || is.null(a$dir)) stop("A layer and direction specification is required for the gwldegree term.")

  # d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed ) {
    stop("The gwldegree term cannot yet handle a nonfixed decay ",
            "term with an attribute. Use fixed=TRUE.", call.=FALSE)
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwldegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    c(.process_layers_degree(nw, a, name="odegree", coef.names=paste("gwldegree#",d,sep=""), inputs=c(d)),
      params=list(gwldegree=NULL,gwldegree.decay=decay),
      map=map, gradient=gradient, conflicts.constraints="degreedist")
  } else {
    if(is.null(a$decay)) stop("Term 'gwldegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwldegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwldegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwldegree_by_attr"
      coef.names <- paste("gwodeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwldegree"
      coef.names <- paste("gwodeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="degreedist")
  }
}
