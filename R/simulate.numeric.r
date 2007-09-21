simulate.numeric <- function(row=NULL, col=NULL, nsim=1, seed=NULL, ...,
                          algorithm.control=list(),
                          verbose=FALSE) {
  out.list <- vector("list", nsim)
# if(is.null(col)){
#  col <- row
# }else{
#  row <- object
# }

  ## Defaults :
  con <- list(boundDeg=NULL, drop=FALSE,
              summarizestats=FALSE
             )

  con[(namc <- names(algorithm.control))] <- algorithm.control
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))

  verb <- match(verbose,
      c("FALSE","TRUE", "very"), nomatch=1)-1

  ##### Begin current code #####
  row.sums<-row
  col.sums<-col
# col.sums<-col.sums[(length(row.sums)+1):length(col.sums)]
  order.col.sums<-order(-col.sums)
  col.sums<-col.sums[order.col.sums]

  nrow<-length(row.sums)
  ncol<-length(col.sums)

  mat <- matrix(0, nrow=nrow, ncol=ncol)
  nedges <- sum(row.sums)
  nsim.in <- 1
  col.sums <- col.sums[col.sums>0]
  row.sums <- row.sums[row.sums>0]
  nrow.nonnull<-length(row.sums)
  ncol.nonnull<-length(col.sums)

  for(i in 1:nsim){
   if(ncol.nonnull > 0 & nrow.nonnull > 0){
    ############################################
    ############################################
    ## Call the C function to create a new    ##
    ## matrix fitting the current constraints ##
    ############################################
    ############################################

    s <- .C("sissim", 
          rowsums=as.integer(row.sums), colsums=as.integer(col.sums),
          newmat=integer(nrow.nonnull*ncol.nonnull), 
          nrow=as.integer(nrow.nonnull), 
          ncol=as.integer(ncol.nonnull), 
          as.double(nsim.in),
          as.integer(verb),
          prob=as.double(1), probvec=double(nsim.in),
          PACKAGE="ergm")

    mat[1:nrow.nonnull, 1:ncol.nonnull] <- matrix(s$newmat,
             ncol=ncol.nonnull)[,order(order.col.sums)]
   }
   out.list[[i]] <- as.network.matrix(mat, matrix.type="bipartite")
  }

  if(nsim > 1){
    out.list <- list(formula = as.formula("~ sis"), 
                     networks = out.list, 
                     stats = NULL, coef=NULL)
    class(out.list) <- "network.series"
  }else{
    out.list <- out.list[[1]]
  }
  return(out.list)
}
