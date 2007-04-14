"ergm.rpvm.run" <- function(nsamp, rpvmbasename,
   SLAVEDIR=NULL, SLAVENAME=NULL,
   PARTAG=1, RESTAG=2, EXITTAG=3,
   maxdeg=32){
  if(missing(SLAVEDIR)){
   SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")
  }
  if(missing(SLAVENAME)){
   SLAVENAME <- paste(rpvmbasename,".slave",sep="")
  }
# administration variables
  slaves <- matrix(0,nrow=0,ncol=2)  # slave table
  slfail <- matrix(0,nrow=0,ncol=2)  # table of failed slaves
  finished <- 0                      # number of finished slaves
# file for storing/changing maxdeg
  fpargrad <- paste(SLAVEDIR,"/.pargrad.",rpvmbasename,sep="")
# file for storing slave table
  fproc <- paste(SLAVEDIR,"/.proc.",rpvmbasename,sep="")
# file for storing failed slaves
  fprocfail <- paste(SLAVEDIR,"/.proc_fail.",rpvmbasename,sep="")
#
# maxdeg  maximal degree of parallelism
# store maxdeg
#
  write(maxdeg, file=fpargrad)

  out <- vector(length=nsamp, mode="numeric")

  # replication loop
  for (it in 1:nsamp) {
    sendflag <- TRUE
    newslave<-TRUE
    while (nrow(slaves) >= maxdeg) {   # maxdeg processes running
      
      # some administration in the waiting time
      failed <- ergm.FindFailedSlaves(slaves[,1]) # look for failed slaves 
      if (length(failed) > 0) {  # if failed slaves found, update tables
        slfail <- rbind(slfail,slaves[slaves[,1] %in% failed,])
        cat(paste("Slave",failed,"failed\n"))
        slaves <- matrix(slaves[!slaves[,1] %in% failed,],ncol=2)
      }
      # read maxdeg from a file and write tables to files
      maxdeg <- scan(file=fpargrad,what=integer(),nlines=1, quiet=TRUE)
      write.matrix(slaves, file=fproc)
      write.matrix(slfail, file=fprocfail)
      # end of administration
      
      if (nrow(slaves) < maxdeg) {newslave<-TRUE; break}  
                                          # wait for results from any slave
      bid <- .PVM.recv (tid = -1, msgtag = RESTAG) 
      info <- .PVM.bufinfo (bid)     # info about the finished slave
      tid <- info[["tid"]]
      rep <- slaves[slaves[,1]==tid,2] # get corresp. replicate
#     cat(paste("Received return from",rep,"\n"))
      result <- .PVM.upkintvec()
      out[rep] <- result[1]
  #
      finished <- finished + 1
#     cat(paste("Finished",finished,"\n"))
      newslave<-FALSE
      if (nrow(slaves) > maxdeg) {  # too many processes running         
        .PVM.send(tid, msgtag=EXITTAG)  # send signal to exit slave
      }                        # remove slave from the table
      slaves<-matrix(slaves[slaves[,1]!=tid,],ncol=2)
    }
    if (newslave) {              # free capacity for spawning slaves
      new <- .PVM.spawnR(slave=paste(SLAVENAME,".r",sep=""),
                         slavedir=SLAVEDIR, outdir=SLAVEDIR)
      if (length (new[new > 0]) == 0) { # spawn failed
        print ("Failed to spawn task.")
        sendflag <- FALSE
      }else{
        tid <- new[1]              # spawn successful
        cat(paste("Slave spawned successfully\n"))
      }
    }
    if (sendflag) {
      # add slave to the table
      slaves <- rbind(slaves,c(tid,it))
      # initialize the send buffer
      .PVM.initsend()
  #
      # pack parameters
      .PVM.pkintvec (c(nsamp,it))
      # send message to slave
      .PVM.send (tid, msgtag=PARTAG)
    }
  } # end of the replication loop
  
  write.matrix(slaves, file=fproc)
  write.matrix(slfail, file=fprocfail)
  # wait for remaining slaves
  nfail <- nrow(slfail)
  it <- finished + 1  
  while (it <= (nsamp-nfail)) {
    bid<-.PVM.recv (tid = -1, msgtag = RESTAG)
    info<-.PVM.bufinfo (bid)
    tid<-info[["tid"]]
    rep <- slaves[slaves[,1]==tid,2]
    result <- .PVM.upkintvec()
    out[rep] <- result[1]
  #
    .PVM.send (tid, msgtag=EXITTAG) # send signal to exit slave
    slaves<-matrix(slaves[slaves[,1]!=tid,],ncol=2) # remove slave from the table
    failed <- ergm.FindFailedSlaves(slaves[,1]) # look for failed slaves  
    if (length(failed) > 0) {  # if failed slaves found, update tables
      slfail<-rbind(slfail,slaves[slaves[,1] %in% failed,])
      slaves<-matrix(slaves[!slaves[,1] %in% failed,],ncol=2)
      nfail<-nrow(slfail)
    }
    write.matrix(slaves, file=fproc)
    write.matrix(slfail, file=fprocfail)
    it<-it+1
  }
  out
}

"ergm.rpvm.clean" <- function(rpvmbasename,
   SLAVEDIR=NULL, SLAVENAME=NULL,
   all=TRUE, Rout=FALSE, indata=FALSE, outdata=FALSE){
  if(missing(SLAVEDIR)){
   SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")
  }
  if(missing(SLAVENAME)){
   SLAVENAME <- paste(rpvmbasename,".slave",sep="")
  }
 if(all){
# cat(paste("/bin/rm -f ",SLAVEDIR,"/.pargrad.",rpvmbasename,"\n",sep=""))
# cat(paste("/bin/rm -f ",SLAVEDIR,"/.proc.",rpvmbasename,"\n",sep=""))
# cat(paste("/bin/rm -f ",SLAVEDIR,"/.proc_fail.",rpvmbasename,"\n",sep=""))
# cat(paste("/bin/rm -f ",SLAVENAME,".r","\n",sep=""))
# cat(paste("/bin/rm -f ",SLAVENAME,".r.Rout","\n",sep=""))
# cat(paste("/bin/rm -fr ",SLAVEDIR,"\n"))
# system(paste("/bin/rm -f ",SLAVEDIR,"/.pargrad.",rpvmbasename,sep=""))
# system(paste("/bin/rm -f ",SLAVEDIR,"/.proc.",rpvmbasename,sep=""))
# system(paste("/bin/rm -f ",SLAVEDIR,"/.proc_fail.",rpvmbasename,sep=""))
# system(paste("/bin/rm -f ",SLAVENAME,".r",sep=""))
# system(paste("/bin/rm -f ",SLAVENAME,".r.Rout",sep=""))
# system(paste("/bin/rmdir --ignore-fail-on-non-empty ",SLAVEDIR))
#
# cat(paste("unlink(",SLAVEDIR,", recursive=TRUE)\n"))
  unlink(SLAVEDIR, recursive=TRUE)
  unlink(paste(getwd(),"/",SLAVENAME,".r",sep=""))
# Next deletes from home
# unlink(paste(Sys.getenv("HOME"),"/",SLAVENAME,"*.r.Rout",sep=""))
  Routs <- paste(system(paste("ls ",Sys.getenv("HOME"),"/",
   SLAVENAME,"*.Rout",sep=""),intern=TRUE,ignore.stderr=TRUE))
  if(length(Routs)>0){
   system(paste("rm -f ",Routs,collapse=" "),ignore.stderr=TRUE)
  }
# system(paste("/bin/rmdir --ignore-fail-on-non-empty ",SLAVEDIR))
 }else{
  if(Rout){
#  system(paste("/bin/rm -f ",SLAVEDIR,"/*.Rout",sep=""))
   unlink(paste(SLAVEDIR,"/*.Rout",sep=""))
  }
  if(indata){
#  system(paste("/bin/rm -f ",SLAVEDIR,"/",SLAVENAME,"*.in.RData",sep=""))
   unlink(paste(SLAVEDIR,"/",SLAVENAME,"*.in.RData",sep=""))
  }
  if(outdata){
#  system(paste("/bin/rm -f ",SLAVEDIR,"/",SLAVENAME,"*.out.RData",sep=""))
   unlink(paste(SLAVEDIR,"/",SLAVENAME,"*.out.RData",sep=""))
  }
 }
 invisible()
}

ergm.FindFailedSlaves <- function (tids) {
  status <- .PVM.pstats(tids)
  tids <- tids[status != "OK"]
  return(tids)
}

"ergm.rpvm.setup" <- function(rpvmbasename,
   SLAVEDIR=NULL, SLAVENAME=NULL,
   verbose=FALSE, aux=FALSE, packagename="ergm", ...){
  if(missing(SLAVEDIR)){
   SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")
  }
  if(missing(SLAVENAME)){
   SLAVENAME <- paste(rpvmbasename,".slave",sep="")
  }
  system(paste("mkdir -p ",SLAVEDIR,sep=""),ignore.stderr=TRUE)
# system(paste("/bin/rm -f ",SLAVEDIR,"/*.Rout",sep=""))
  Routs <- paste(system(paste("ls ",SLAVEDIR,"/*.Rout",sep=""),intern=TRUE,ignore.stderr=TRUE))
  if(length(Routs)>0){
   unlink(paste(SLAVEDIR,"/*.Rout",sep=""))
  }
# system(paste("ln -fs ",getwd(),"/",SLAVENAME,".r ",SLAVEDIR,"/",SLAVENAME,".r",sep=""))
  file.symlink(paste( getwd(),"/",SLAVENAME,".r",sep=""),
               paste(SLAVEDIR,"/",SLAVENAME,".r",sep=""))
#
  if(verbose){
   cat("The slave directory is set to", SLAVEDIR, "\n")
  }
#
# Start PVM if necessary
#
  PVM.running <- try(.PVM.config(), silent=TRUE)
  if(inherits(PVM.running,"try-error")){
   hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
   .PVM.start.pvmd(hostfile)
   cat("no problem... PVM started by R...\n")
  }
#
# enrolling the process into PVM
#
  mytid <- .PVM.mytid ()
#
# Write the slave file
#
  write(matrix(c(
   'options(echo = TRUE)',
   paste('setwd("',getwd(),'")',sep=''),
   paste('rpvmbasename <- "',rpvmbasename,'"',sep=""),
   'library(rpvm)',
   'library(rsprng)',
   'library(MASS)',
   'PARTAG <- 1; RESTAG <- 2; EXITTAG <- 3',
   'SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")',
   'load(file=paste(SLAVEDIR,"/",rpvmbasename,".common.RData",sep=""))',
   'mytid <- .PVM.mytid ()',
   'parentid <- .PVM.parent ()',
   'repeat {',
   '  # wait for a message from the parent with any tag',
   '  bid <- .PVM.recv (tid = parentid, msgtag = -1)',
   '  # get info about the message',
   '  info<-.PVM.bufinfo (bid)',
   '  if (info[["msgtag"]] == EXITTAG) break  # exit the process',
   '#',
   '# Input parameters',
   '#',
   '  iparam <- .PVM.upkintvec()',
   '  # establish RNG state',
   '  r <- iparam[1]; it<-iparam[2]',
   '  oldrng<-init.sprng(nstream=r,streamno=it-1,seed=1)',
   '  cat("Iteration",it,"of",r,"\\n")',
#  paste('library(',packagename,')',sep=""),
   paste('library(','statnet',')',sep=""),
   '        z <- .C("MCMC_wrapper",',
   '               as.double(Clist$heads), as.double(Clist$tails), ',
   '               as.double(Clist$nedges), as.double(Clist$n),',
   '               as.integer(Clist$dir), as.double(Clist$bipartite),',
   '               as.integer(Clist$nterms), ',
   '               as.character(Clist$fnamestring),',
   '               as.character(Clist$snamestring),',
   '               as.character(MHproposal$name), as.character(MHproposal$package),',
   '               as.double(Clist$inputs), as.double(eta0),',
   '               as.double(MCMCparams.parallel$samplesize),',
   '               s = as.double(t(MCMCparams.parallel$stats)),',
   '               as.double(MCMCparams.parallel$burnin), as.double(MCMCparams.parallel$interval), ',
   '               newnw = integer(maxedges), ',
   '               as.integer(verbose), as.integer(BD$attribs), ',
   '               as.integer(BD$maxout), as.integer(BD$maxin),',
   '               as.integer(BD$minout), as.integer(BD$minin),',
   '               as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), ',
   '               as.double(maxedges),',
   '               as.double(MCMCparams.parallel$Clist.miss$heads),',
   '               as.double(MCMCparams.parallel$Clist.miss$tails),',
   '               as.double(MCMCparams.parallel$Clist.miss$nedges),',
   '               PACKAGE="statnet")',
   '# save the results',
   '  z <- list(s=z$s, newnw=z$newnw)',
   '  save(z,',
   '       file=paste(SLAVEDIR,"/",rpvmbasename,".out.",it,".RData",sep=""))',
   '  # initialize the send buffer',
   '  .PVM.initsend()',
   '  # pack and send the iteration back as a guide',
   '  .PVM.pkintvec(it)',
   '  .PVM.send (parentid, msgtag=RESTAG)',
   '  free.sprng(oldrng)',
   '} # end of repeat',
   '.PVM.exit()'
  ),ncol=1), file=paste(SLAVENAME,".r",sep=""), ncolumns=1)
 list(mytid=mytid, SLAVEDIR=SLAVEDIR, SLAVENAME=SLAVENAME,
       pvm.PARTAG=1, pvm.RESTAG=2, pvm.EXITTAG=3)
}

"ergm.rpvm.setup.dyn" <- function(rpvmbasename,
   SLAVEDIR=NULL, SLAVENAME=NULL,
   verbose=FALSE, aux=FALSE, packagename="ergm", ...){
  if(missing(SLAVEDIR)){
   SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")
  }
  if(missing(SLAVENAME)){
   SLAVENAME <- paste(rpvmbasename,".slave",sep="")
  }
  system(paste("mkdir -p ",SLAVEDIR,sep=""),ignore.stderr=TRUE)
# system(paste("/bin/rm -f ",SLAVEDIR,"/*.Rout",sep=""))
  Routs <- paste(system(paste("ls ",SLAVEDIR,"/*.Rout",sep=""),intern=TRUE,ignore.stderr=TRUE))
  if(length(Routs)>0){
   unlink(paste(SLAVEDIR,"/*.Rout",sep=""))
  }
# system(paste("ln -fs ",getwd(),"/",SLAVENAME,".r ",SLAVEDIR,"/",SLAVENAME,".r",sep=""))
  file.symlink(paste( getwd(),"/",SLAVENAME,".r",sep=""),
               paste(SLAVEDIR,"/",SLAVENAME,".r",sep=""))
#
  if(verbose){
   cat("The slave directory is set to", SLAVEDIR, "\n")
  }
#
# Start PVM if necessary
#
  PVM.running <- try(.PVM.config(), silent=TRUE)
  if(inherits(PVM.running,"try-error")){
   hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
   .PVM.start.pvmd(hostfile)
   cat("no problem... PVM started by R...\n")
  }
#
# enrolling the process into PVM
#
  mytid <- .PVM.mytid ()
#
# Write the slave file
#
  write(matrix(c(
   'options(echo = TRUE)',
   paste('setwd("',getwd(),'")',sep=''),
   paste('rpvmbasename <- "',rpvmbasename,'"',sep=""),
   'library(rpvm)',
   'library(rsprng)',
   'library(MASS)',
   'PARTAG <- 1; RESTAG <- 2; EXITTAG <- 3',
   'SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")',
   'load(file=paste(SLAVEDIR,"/",rpvmbasename,".common.RData",sep=""))',
   'mytid <- .PVM.mytid ()',
   'parentid <- .PVM.parent ()',
   'repeat {',
   '  # wait for a message from the parent with any tag',
   '  bid <- .PVM.recv (tid = parentid, msgtag = -1)',
   '  # get info about the message',
   '  info<-.PVM.bufinfo (bid)',
   '  if (info[["msgtag"]] == EXITTAG) break  # exit the process',
   '#',
   '# Input parameters',
   '#',
   '  iparam <- .PVM.upkintvec()',
   '  # establish RNG state',
   '  r <- iparam[1]; it<-iparam[2]',
   '  oldrng<-init.sprng(nstream=r,streamno=it-1,seed=1)',
   '  cat("Iteration",it,"of",r,"\\n")',
#  paste('library(',packagename,')',sep=""),
   paste('library(','statnet',')',sep=""),
   '      z <- .C("MCMCDyn_wrapper",',
   '          as.integer(Clist.dissolve$order.code),',
   '          as.double(Clist$heads), as.double(Clist$tails), ',
   '          as.double(Clist$nedges), as.double(Clist$n),',
   '          as.integer(Clist$dir), as.double(Clist$bipartite),',
   '          as.integer(Clist$nterms), ',
   '          as.character(Clist$fnamestring),',
   '          as.character(Clist$snamestring),',
   '          as.character(MHproposal$name), as.character(MHproposal$package),',
   '          as.double(Clist$inputs), as.double(eta0),',
   '          as.integer(Clist.dissolve$nterms),',
   '          as.character(Clist.dissolve$fnamestring),',
   '          as.character(Clist.dissolve$snamestring),',
   '          as.double(Clist.dissolve$inputs),',
   '          as.double(MCMCparams.parallel$samplesize),',
   '          s = as.double(t(MCMCparams.parallel$stats)),',
   '          as.double(MCMCparams.parallel$burnin), as.double(MCMCparams.parallel$interval),',
   '          newnwhead = integer(maxchanges), newnwtail = integer(maxchanges), ',
   '          diffnwtime = integer(maxchanges), diffnwhead = integer(maxchanges), diffnwtail = integer(maxchanges), ',
   '          dissnwtime = integer(maxchanges), dissnwhead = integer(maxchanges), dissnwtail = integer(maxchanges),',
   '          as.integer(verbose), ',
   '          as.double(MCMCparams.parallel$gamma), as.integer(MCMCparams.parallel$dyninterval), ',
   '          as.integer(BD$attribs), ',
   '          as.integer(BD$maxout), as.integer(BD$maxin),',
   '          as.integer(BD$minout), as.integer(BD$minin),',
   '          as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), ',
   '          as.double(maxchanges),',
   '          as.double(0.0), as.double(0.0),',
   '          as.double(0.0),',
   '          PACKAGE="statnet")',
   '# save the results',
   '  z <- list(s=z$s,newnwhead=z$newnwhead,newnwtail=z$newnwtail,diffnwhead=z$diffnwhead,dissnwhead=z$dissnwhead,',
   '            diffnwtail=z$diffnwtail, dissnwtail=z$dissnwtail,',
   '            diffnwtime=z$diffnwtime, dissnwtime=z$dissnwtime)',
   '  save(z,',
   '       file=paste(SLAVEDIR,"/",rpvmbasename,".out.",it,".RData",sep=""))',
   '  # initialize the send buffer',
   '  .PVM.initsend()',
   '  # pack and send the iteration back as a guide',
   '  .PVM.pkintvec(it)',
   '  .PVM.send (parentid, msgtag=RESTAG)',
   '  free.sprng(oldrng)',
   '} # end of repeat',
   '.PVM.exit()'
  ),ncol=1), file=paste(SLAVENAME,".r",sep=""), ncolumns=1)
 list(mytid=mytid, SLAVEDIR=SLAVEDIR, SLAVENAME=SLAVENAME,
       pvm.PARTAG=1, pvm.RESTAG=2, pvm.EXITTAG=3)
}
"ergm.rpvm.miss.setup" <- function(rpvmbasename,
   SLAVEDIR=NULL, SLAVENAME=NULL,
   verbose=FALSE,  aux=FALSE, packagename="ergm", ...){
  if(missing(SLAVEDIR)){
   SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")
  }
  if(missing(SLAVENAME)){
   SLAVENAME <- paste(rpvmbasename,".slave",sep="")
  }
  system(paste("mkdir -p ",SLAVEDIR,sep=""),ignore.stderr=TRUE)
# system(paste("/bin/rm -f ",SLAVEDIR,"/*.Rout",sep=""))
  unlink(paste(SLAVEDIR,"/*.Rout",sep=""))
# system(paste("ln -fs ../",SLAVENAME,".r ",SLAVEDIR,"/",SLAVENAME,".r",sep=""))
  file.symlink(paste( getwd(),"/",SLAVENAME,".r",sep=""),
               paste(SLAVEDIR,"/",SLAVENAME,".r",sep=""))
#
  if(verbose){
   cat("The slave directory is set to", SLAVEDIR, "\n")
  }
#
# Start PVM if necessary
#
  PVM.running <- try(.PVM.config(), silent=TRUE)
  if(inherits(PVM.running,"try-error")){
   hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
   .PVM.start.pvmd(hostfile)
   cat("no problem... PVM started by R...\n")
  }
#
# enrolling the process into PVM
#
  mytid <- .PVM.mytid ()
#
# Write the slave file
#
  write(matrix(c(
   'options(echo = TRUE)',
   paste('setwd("',getwd(),'")',sep=''),
   paste('rpvmbasename <- "',rpvmbasename,'"',sep=""),
   'library(rpvm)',
   'library(rsprng)',
   'PARTAG <- 1; RESTAG <- 2; EXITTAG <- 3',
   'SLAVEDIR <- paste(getwd(),"/",rpvmbasename,".run",sep="")',
   'load(file=paste(SLAVEDIR,"/",rpvmbasename,".common.RData",sep=""))',
   'mytid <- .PVM.mytid ()',
   'parentid <- .PVM.parent ()',
   'repeat {',
   '  # wait for a message from the parent with any tag',
   '  bid <- .PVM.recv (tid = parentid, msgtag = -1)',
   '  # get info about the message',
   '  info<-.PVM.bufinfo (bid)',
   '  if (info[["msgtag"]] == EXITTAG) break  # exit the process',
   '#',
   '# Input parameters',
   '#',
   '  iparam <- .PVM.upkintvec()',
   '  # establish RNG state',
   '  r <- iparam[1]; it<-iparam[2]',
   '  oldrng<-init.sprng(nstream=r,streamno=it-1,seed=1)',
   '  cat("Iteration",it,"of",r,"\\n")',
   paste('library(',packagename,')',sep=""),
  switch(aux+1,
    '       z <- .C("MCMC_wrapper",',
    '       z <- .C("MCMCaux_wrap",'
        ),
   '        as.double(Clist.expanded$heads), as.double(Clist.expanded$tails),',
   '        as.double(Clist.expanded$nedges), as.double(Clist.expanded$n),',
   '        as.integer(Clist.expanded$notdir), as.integer(Clist.expanded$noptions),',
   '        as.character(Clist.expanded$fnamestring),',
   '        as.character(Clist.expanded$snamestring),',
   '        as.double(Clist.expanded$inputs),',
   '        as.double(eta0),',
  switch(aux+1, NULL,
   '        as.double(etaprop0),'
        ),
   '        as.double(MCMCsamplesize.parallel),',
   '        s = double(MCMCsamplesize.parallel * Clist.expanded$nparam),',
  switch(aux+1, NULL,
   '        stheta = double(MCMCsamplesize.parallel * Clist.expanded$nparam),'
        ),
   '        as.double(burnin), as.double(interval),',
   '        as.double(distanceMetric), quote(prepmahal(stats, dimension)),',
   '        new.env(), newg = double(maxedges),',
   '        as.integer(1),',
   '        as.integer(attribs),',
   '        as.integer(maxout), as.integer(maxin), as.integer(minout),',
   '        as.integer(minin), as.integer(condAllDegExact),',
   '        as.integer(length(attribs)),',
   '        as.integer(multiplicity.constrained), as.integer(proposalnumber.constrained),',
   '        as.double(Clist.miss$heads), as.double(Clist.miss$tails),',
   '        as.double(Clist.miss$nedges), as.integer(Clist.miss$dir),',
   '        as.double(fixterms),',
   '        as.double(maxedges),',
#  paste('        PACKAGE="',packagename,'"',sep=""),
   '        PACKAGE="ergm"',
   '       )',
   '# save the results',
  switch(aux+1,
    '  z <- list(s=z$s, newg=z$newg)',
    '  z <- list(s=z$s, stheta=z$stheta, newg=z$newg)',
        ),
   '  save(z,',
   '       file=paste(SLAVEDIR,"/",rpvmbasename,".out.",it,".RData",sep=""))',
   '  # initialize the send buffer',
   '  .PVM.initsend()',
   '  # pack and send the iteration back as a guide',
   '  .PVM.pkintvec(it)',
   '  .PVM.send (parentid, msgtag=RESTAG)',
   '  free.sprng(oldrng)',
   '} # end of repeat',
   '.PVM.exit()'
  ),ncol=1), file=paste(SLAVENAME,".r",sep=""), ncolumns=1)
 list(mytid=mytid, SLAVEDIR=SLAVEDIR, SLAVENAME=SLAVENAME,
       pvm.PARTAG=1, pvm.RESTAG=2, pvm.EXITTAG=3)
}
