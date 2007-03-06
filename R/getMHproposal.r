#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

getMHproposal <- function(name, arguments, nw, model) {
  proposal <- eval(call(paste("InitMHP.", name, sep=""),
                        arguments, nw, model))
#  nwlist <- proposal$nwlist
proposal
}



