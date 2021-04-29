#' @note For directed networks, only outgoing two-path ("OTP") shared partners
#'   are counted.  In other words, for a #2 in a directed graph, the number of
#'   shared partners counted by `#1` is the number of nodes `k` that have edges
#'   `i -> k -> j`.  (These may also be called homogeneous shared partners.)  To
#'   count other types of shared partners instead, see `#3`.
