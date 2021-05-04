#' @note For directed networks, only outgoing two-path ("OTP") shared partners
#'   are counted.  In other words, for a <%= kind %> in a directed graph, the number of
#'   shared partners counted by `<%= fn %>` is the number of nodes `k` that have edges
#'   `i -> k -> j`.  (These may also be called homogeneous shared partners.)  To
#'   count other types of shared partners instead, see `<%= see %>`.
