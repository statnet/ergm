#' @note While there is only one shared partner configuration in the undirected
#'   case, nine distinct configurations are possible for directed graphs, selected
#'   using the `type` argument. Currently, terms may be defined with respect to
#'   five of these configurations; they are defined here as follows (using
#'   terminology from Butts (2008) and the `relevent` package):
#' -  Outgoing Two-path (`"OTP"`): vertex \eqn{k} is an OTP shared partner of ordered
#'    pair \eqn{(i,j)} iff \eqn{i \to k \to j}{i->k->j}. Also known as "transitive
#'    shared partner".
#' - Incoming Two-path (`"ITP"`): vertex \eqn{k} is an ITP shared partner of ordered
#'    pair \eqn{(i,j)} iff \eqn{j \to k \to i}{j->k->i}. Also known as "cyclical shared
#'    partner"
#' - Reciprocated Two-path (`"RTP"`): vertex \eqn{k} is an RTP shared partner of ordered
#'    pair \eqn{(i,j)} iff \eqn{i \leftrightarrow k \leftrightarrow j}{i<->k<->j}.
#' - Outgoing Shared Partner (`"OSP"`): vertex \eqn{k} is an OSP shared partner of
#'    ordered pair \eqn{(i,j)} iff \eqn{i \to k, j \to k}{i->k, j->k}.
#' - Incoming Shared Partner (`"ISP"`): vertex \eqn{k} is an ISP shared partner of ordered
#'    pair \eqn{(i,j)} iff \eqn{k \to i, k \to j}{k->i, k->j}.
#' By default, outgoing two-paths (`"OTP"`) are calculated. Note that Robins et al. (2009)
#'   define closely related statistics to several of the above, using slightly different terminology.
