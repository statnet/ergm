ergm.revisetheta0 <- function(m, theta0) {
  for(i in 1:length(m$terms)) {
    if (!is.null(m$terms[[i]]$params)) {
      n1=names(m$terms[[i]]$params)
      n2=names(theta0)
      overlap = (n2 %in% n1)
      if (length(n1)>sum(overlap)) {
        before = (cumsum(overlap)==0)
        after = (!before & overlap==0)
#  It is assumed that the terms in m$terms[[i]]$params match the
# terms in theta0 only in a contiguous region of theta0.  This region
# of overlap is expanded to include any new terms in m$terms[[i]]$params,
# along with their values.  This expansion is very naive, simply accomplished
# by concatenating the new terms onto the end of the overlapping terms.
        newoverlap=c(theta0[overlap],rep(0,length(n1)-sum(overlap)))
        names(newoverlap)=n1
        for (j in (1+sum(overlap)):length(n1)) {
          newoverlap[j] = m$terms[[i]]$params[[j]]
        }
        theta0=c(theta0[before],newoverlap,theta0[after])
      }
    }
  }
  theta0
}
      















































































