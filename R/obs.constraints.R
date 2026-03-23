#  File R/obs.constraints.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

dot_sub_constraints <- function(nw, ..., default.dot = c("first", "last", "none")) {
  default.dot <- match.arg(default.dot)

  # Ensure that each argument is a term_list with LHS embedded.
  tll <- list(...) |> map(ergm_flatten_conterm_list) |> compact()
  # If no missing edges, remove the "observed" constraint.
  if (network.naedgecount(nw) == 0L) tll <- map(tll, .delete_term, "observed")

  if(all(lengths(tll) == 0L)) return(NULL)

  # Go through the constraint lists, substituting the earlier ones into the dots in the later ones.
  otl <- tll[[1]]
  for (i in seq_along(tll)[-1]) {
    ntl <- tll[[i]]

    # Find the substitution positions.
    pos <- which(ntl |> map_chr(\(x) as.character(x)[1]) == ".")

    # Don't substitute at negative dots.
    pos_sign <- attr(ntl, "sign")[pos]
    pos <- pos[pos_sign>0]

    # If no dots, use default behaviour unless there is a negative dot.
    if (!length(pos) && all(pos_sign > 0))
      ntl <- switch(default.dot,
                    first = c(otl, ntl),
                    last = c(ntl, otl),
                    none = ntl)

    # Substitute (working backwards, to prevent pos from changing).
    for(p in rev(pos))
      ntl <- c(ntl[seq_len(pos-1)], otl, ntl[-seq_len(pos)])

    otl <- ntl
  }

  # Delete remaining dots (including the negative ones).
  .delete_term(otl, ".")
}


.handle.auto.constraints <- function(nw,
                                     constraints=trim_env(~.),
                                     obs.constraints=trim_env(~.-observed),
                                     target.stats=NULL,
                                     control = control.ergm(),
                                     control.prop = "MCMC",
                                     default.dot=c("first", "last", "none")){
  default.dot <- match.arg(default.dot)

  tl <- dot_sub_constraints(nw,
                            nw%ergmlhs%"constraints",
                            c(enlist(constraints),
                              enlist(control[[paste0(control.prop, ".prop")]])),
                            default.dot = default.dot)

  obs.tl <- dot_sub_constraints(nw,
                                nw%ergmlhs%"obs.constraints",
                                c(enlist(obs.constraints),
                                  enlist(control[[paste0("obs.", control.prop, ".prop")]])),
                                default.dot = default.dot)

  # Do any of the observational constraints formulas have terms?
  obs.tl <- if (length(obs.tl)) {
              # Observation process handling only needs to happen if
              # the sufficient statistics are not specified. If the
              # sufficient statistics are specified, the nw's dyad
              # states are irrelevant.
              if(!is.null(target.stats)){
                message("Target statistics specified in a network with missing dyads and/or a nontrivial observation process. Since (by sufficiency) target statistics provide all the information needed to fit the model, missingness and observation process will not affect estimation.")
                if(network.naedgecount(nw)) nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
                NULL
              } else # Prepend the sample space constraint, but don't propagate .select().
                c(obs.tl, .delete_term(tl, ".select"))
            }

  list(nw = nw, conterms = tl, conterms.obs = obs.tl)
}

has.obs.constraints <- function(...) length(.handle.auto.constraints(...)$conterms.obs) > 0

#' @importFrom rlang rep_named
.embed.target.stats <- function(model, target.stats){
  replace(rep_named(param_names(model, canonical=TRUE), NA_real_),
          !model$etamap$offsetmap,
          target.stats)
}
