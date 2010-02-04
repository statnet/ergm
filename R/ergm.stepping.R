#  File ergm/R/ergm.stepping.R  Last updated by Ruth, 1 February 2010
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
#  MCMCparams=c(control,
#   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,stepMCMCsize=stepMCMCsize,
#        gridsize=gridsize, maxit=maxit,Clist.miss=Clist.miss, mcmc.precision=control$mcmc.precision))

ergm.stepping = function(theta0, nw, model, Clist, initialfit, 
#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
## approx=lognormapprox, filename.prefix=NULL, plots=FALSE,  # currently useless, but plots can be reimplemented
				MCMCparams=MCMCparams, 
				MHproposal=MHproposal, MHproposal.miss=MHproposal.miss, 
				verbose=FALSE, estimate=TRUE, sequential=TRUE,...){ 

#   preliminary, to set up structure. 
## Relevant only if maxit>1:					iteration <- 1
					nw.orig <- nw
					asyse=theta0-theta0
					mc.se=1+0.05*asyse
					mle.lik=initialfit$mle.lik
					v=list(coef=theta0)
					theta.original=theta0
					thetaprior=theta0-theta0
					stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)					
					stats[1,] <- Clist$obs - Clist$meanstats
					MCMCparams$stats <- stats
					MCMCparams$meanstats <- Clist$meanstats
					
###while(iteration <= MCMCparams$maxit){
					thetaprior <- theta0
					theta0 <- v$coef
					eta0 <- ergm.eta(theta0, model$etamap)
					
## Relevant only if maxit>1:
#					if(verbose){
#						cat("Iteration ",iteration," of at most ", MCMCparams$maxit,
#							" with parameter: \n", sep="")
#						print(theta0)
#					}else{
#						cat("Iteration ",iteration," of at most ", MCMCparams$maxit,": ",sep="")
#					}
					
#					Clist <- ergm.Cprepare(nw, model) # probably does not need to be here
						
#   Takes MCMC sample from the initial value, theta0 (technically, from its etamap) 
#   for the first iteration.
#####					z <- ergm.getMCMCsample(Clist, MHproposal, eta0, MCMCparams, verbose)
#####					statsmatrix <- z$statsmatrix
#####					colnames(statsmatrix) <- model$coef.names
#####					v$sample <- statsmatrix
	
					statsmatrix=simulate.formula(model$formula, nsim=MCMCparams$stepMCMCsize, 
											theta0=theta0, burnin=MCMCparams$burnin, 
											interval=MCMCparams$interval, statsonly=T)
					colnames(statsmatrix) <- model$coef.names
					v$sample <- statsmatrix
#					sampmeans <- colMeans(xsim)
					
#  Does the same, if missing edges:					
					if(MCMCparams$Clist.miss$nedges > 0){
							Clist <- ergm.Cprepare(nw, model)
							statsmatrix.miss <- ergm.getMCMCsample(Clist, MHproposal.miss, eta0, MCMCparams, verbose)$statsmatrix
							colnames(statsmatrix.miss) <- model$coef.names
							if(verbose){cat("Back from constrained MCMC...\n")}
						}						
					else{
							statsmatrix.miss <- NULL
							if(verbose){cat("Back from unconstrained MCMC...\n")}
						}					

#  ?? What is "sequential"?  (Default is that it is TRUE)  Is this relevant?
#					if(sequential & MCMCparams$Clist.miss$nedges == 0){
#							nw <- network.update(nw, z$newedgelist, "edgelist")
#							nw.obs <- summary(model$formula, basis=nw)
#							namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
#							MCMCparams$stats[1,] <- nw.obs[namesmatch]-MCMCparams$meanstats
#						}			
					
# Next iteration begins here:
#					iteration <- iteration + 1
	
#  [skipped some stuff from mainfitloop]  -- should the stepping happen here?  Dont think so.
#  Print some updates if required by "verbose":  
					if(verbose){
#                           (only use if sequential part is used and is generating a new "nw")
#							cat(paste("The density of the returned network is",
#									  network.density(nw),"\n"))
							cat(paste("The density of the original network is",
									  network.density(nw.orig),"\n"))
							cat("Summary of simulation, relative to observed network:\n")
							print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
							degreedist(nw)
							cat("Meanstats of simulation, relative to observed network:\n")
							print(summary(model$formula, basis=nw)-Clist$meanstats)
						}
					if(verbose){cat("Calling optimization routines...\n")}

# If no estimate is called for, skip otimizing and return the preliminary info.						
					if(verbose) {cat("statseval:  epsilon value is ",MCMCparams$epsilon,"\n")}
					if(!estimate){
							if(verbose){cat("Skipping optimization routines...\n")}
							l <- list(coef=theta0, mc.se=rep(NA,length=length(theta0)),
									  sample=statsmatrix, sample.miss=statsmatrix.miss,
									  iterations=1, MCMCtheta=theta0,
									  loglikelihood=NA, #mcmcloglik=NULL, 
									  mle.lik=NULL,
									  gradient=rep(NA,length=length(theta0)), #acf=NULL,
									  samplesize=MCMCparams$samplesize, failure=TRUE,
									  newnetwork = nw)
							return(structure (l, class="ergm"))
						}
						
#  Print some updates if required by "verbose":
					if(verbose){cat("Calling MCMLE Optimization...\n")}
					if(verbose){cat("Using Newton-Raphson Step ...\n")}

#  The actual optimization and stepping mechanisms:
#  (If not the last iteration do not compute all the extraneous
#   statistics that are not needed until output)
	thetaSTEPPED = theta0   #change to eta0?
					stepiter=0
					finished = FALSE
#check if obs.stats are in the convex hull of the sample taken already:
#					if(is.inCH(obs, statsmatrix)){finished=TRUE}
					if(is.inCH(Clist$obs, statsmatrix)){
						finished=TRUE
						cat(paste("Original observed statistics are already in the convex hull - ","\n"))
						cat(paste("Stepping procedure is unnecessary and is skipped.","\n"))					
					}
#Begins stepping loop:
					while (!finished) { # Iterate until alpha==1
						inCH=FALSE
						print(1000)
# If gyobs is not in convex hull, then use stepping procedure to move the fake gyobs, xi, closer:
						stepiter=stepiter+1
						hi <- MCMCparams$gridsize  # Goal: Let alpha be largest possible multiple of .01
						lo <- alph <- 0
						cat("Stepping Iteration #",stepiter,"\n")
						flush.console()
# Finds the steplength that puts xi in the convex hull of the sample:
						while (hi-lo>1 || hi > alph) {
							cat("  Trying alpha=", alph<-ceiling((hi+lo)/2), "/", MCMCparams$gridsize, "\n")
							flush.console()
							alpha = alph/MCMCparams$gridsize
							print(alph)
							xi = alpha*Clist$meanstats  + (1-alpha)*colMeans(statsmatrix[,])
#browser()							
							cat(paste("xi",xi,"\n"))
							
							inCH=is.inCH(xi, statsmatrix)
							if (inCH)
							lo=alph
## If 1/gridsize is not small enough, function fails and gives error message
							else if (alph==1)
							stop("alpha=", 1/MCMCparams$gridsize, "still not small enough for MLE.N.  Try using a", 
								 " gridsize larger than ", MCMCparams$gridsize)
							else
							hi=alph
						}
						if (!inCH) {# Last attempt was a fail, so decrease alph by one
							alph = alph-1
							alpha = alph/MCMCparams$gridsize
							xi = alpha*Clist$meanstats  + (1-alpha)*colMeans(statsmatrix[,])	
						}
						print(4000)
## When the stepped xi is in the convex hull, find the MLE for gyobs=xi
# I think I made this for the stepped xi instead of gyobs by changing ergm.estimate (now I can say "xobs=xi"),
# but how do I specify the correct theta0?  And am I specifying the right sample?
						v<-ergm.estimate(theta0=thetaSTEPPED, model=model, xobs=xi,  
										 statsmatrix=statsmatrix, 
										 statsmatrix.miss=statsmatrix.miss, 
										 epsilon=MCMCparams$epsilon,
										 nr.maxit=MCMCparams$nr.maxit,
										 nr.reltol=MCMCparams$nr.reltol,
										 calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
										 trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
										 compress=MCMCparams$compress, verbose=verbose,
										 estimateonly=TRUE)
						print(5000)
## Once a new stepped MLE is found, generate a new sample from it:
#   Takes MCMC sample from the initial value, theta0 (technically, from its etamap) 
#   for the first iteration. 
#?						theta0 <- v$theta0   ## does this take the MLE from the last xi?
#?						cat(paste("Names of new v", names(v),"\n"))
						
			thetaSTEPPED<-v$coef			
						cat(paste("v$coef", v$coef, "\n"))
						cat(paste("MCMCtheta", v$MCMCtheta, "\n"))
										
						
#?						eta0 <- ergm.eta(theta0, model$etamap)
#?						cat(paste("Eta", eta, "\n"))
#####						z <- ergm.getMCMCsample(Clist, MHproposal, eta0, MCMCparams, verbose)
#####						statsmatrix <- z$statsmatrix
## If alpha is 1, then latest eta is a valid MLE, but take final sample with larger MCMCparams$MCMCsamplesize sample size.  
## If alpha is not 1, this sample will be the starting sample for the next loop.
						nsimhere=MCMCparams$stepMCMCsize
						if(alph==MCMCparams$gridsize){nsimhere=MCMCparams$MCMCsamplesize
						cat("Now ending with one large sample.\n")
						flush.console()
						}
						statsmatrix=simulate.formula(model$formula, nsim=nsimhere,   ##should I be using the eta?
													theta0=thetaSTEPPED, burnin=MCMCparams$burnin, 
													interval=MCMCparams$interval, statsonly=T)
						colnames(statsmatrix) <- model$coef.names
						v$sample <- statsmatrix
#					sampmeans <- colMeans(xsim)						
#?						cat(paste("alph is ", alph, "\n"))
						
##  Does the same, if missing edges:					
#								if(MCMCparams$Clist.miss$nedges > 0){
#									Clist <- ergm.Cprepare(nw, model)
#									statsmatrix.miss <- ergm.getMCMCsample(Clist, MHproposal.miss, eta0, MCMCparams, verbose)$statsmatrix
#									colnames(statsmatrix.miss) <- model$coef.names
#									if(verbose){cat("Back from constrained MCMC...\n")}
#								}						
#								else{
#									statsmatrix.miss <- NULL
#									if(verbose){cat("Back from unconstrained MCMC...\n")}
#								}					
#								
##  ?? What is "sequential"?  (Default is that it is TRUE)
#								if(sequential & MCMCparams$Clist.miss$nedges == 0){
#									nw <- network.update(nw, z$newedgelist, "edgelist")
#									nw.obs <- summary(model$formula, basis=nw)
#									namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
#									MCMCparams$stats[1,] <- nw.obs[namesmatch]-MCMCparams$meanstats
#								}
						
## If alpha is still not 1, go back to next iteration
## If alpha is 1, then latest eta is a valid MLE and has a sample of size MCMCparams$MCMCsamplesize with it.  
						
						finished = (alph==MCMCparams$gridsize)
					} # Ends "while (!finished)"
					cat(paste("FINISHED","\n"))
					
#					}	# Ends "while(iteration <= MCMCparams$maxit)"


					print(7000)
## Do one last check that the observed stats are truly in the new larger sample convex hull, and end break from function if
# that is the case.
	inCH=is.inCH(Clist$obs, statsmatrix)
	if (!inCH)
	stop("Observed statistics are not in final sample convex hull.")
					v<-ergm.estimate(theta0=theta0, model=model,                 #### xobs?
									 statsmatrix=statsmatrix, statsmatrix.miss=statsmatrix.miss, 
									 epsilon=MCMCparams$epsilon,
									 nr.maxit=MCMCparams$nr.maxit,
									 nr.reltol=MCMCparams$nr.reltol,
									 calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
									 trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
									 compress=MCMCparams$compress, verbose=verbose)

					mle.lik <- mle.lik + abs(v$loglikelihood)
					v$newnetwork <- nw
					v$burnin <- MCMCparams$burnin
					v$samplesize <- MCMCparams$samplesize
					v$interval <- MCMCparams$interval
					v$network <- nw.orig
					v$newnetwork <- nw
					v$interval <- MCMCparams$interval
					v$theta.original <- theta.original
					v$mplefit <- initialfit
					v$parallel <- MCMCparams$parallel

					if(!v$failure & !any(is.na(v$coef))){
					#     asyse <- sqrt(diag(robust.inverse(-v$hessian)))
					#     asyse <- try(sqrt(diag(robust.inverse(cov(statsmatrix)))))
						asyse <- mc.se
						options(warn=-1)
						if(is.null(v$covar)){
							asyse[names(v$coef)] <- sqrt(diag(robust.inverse(-v$hessian)))
						}else{
							asyse[names(v$coef)] <- sqrt(diag(v$covar))
						}
						options(warn=0)
					}

					endrun <- MCMCparams$burnin+MCMCparams$interval*(MCMCparams$samplesize-1)
					attr(v$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
					attr(v$sample, "class") <- "mcmc"
					# v$null.deviance <- v$mplefit$glm$null.deviance
					v$null.deviance <- 2*network.dyadcount(nw.orig)*log(2)
					# v$mle.lik <- -v$mplefit$glm$deviance/2 + v$loglikelihood
					# v$mle.lik <- -v$null.deviance/2 + v$loglikelihood
					# v$mle.lik <- -initialfit$mle.lik/2 + v$loglikelihood
					v$mle.lik <- mle.lik
					# This next is the right one
					# v$mcmcloglik <- v$mcmcloglik - network.dyadcount(nw.orig)*log(2)
					v$etamap <- model$etamap
					# if(!is.na(v$mle.lik) && v$mle.lik>0){
					#   v$mle.lik <- -v$mplefit$glm$deviance/2 
					# }  #I am really not sure about this mle.lik value.
					v
				}  # Ends the whole function									
							
							

