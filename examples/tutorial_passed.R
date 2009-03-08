install.packages('statnet')
library(statnet)
#setwd('full.path.for.the.folder')
install.packages('coda')			# install package from CRAN
library(coda)					# attach installed package

a <- 3							# assignment
a								# evaluation

sqrt(a)						# perform an operation
b <- sqrt(a)						# perform operation and save
b							

ls()								# list objects in global enviro
help(sqrt)						# help w/ functions
?sqrt							# same thing
help.start()						# lots of help

rm(a)							# remove an object
library(ergm)
library(sna)

library(statnet)
data()							# tells us the datasets in our packages
data(florentine)					# loads flomarriage & flobusiness data
flomarriage						# Let's look at the flomarriage data
plot(flomarriage)					# Let's view the flomarriage plot
flomodel.01 <- ergm(flomarriage~edges)			# fit model
flomodel.01							# look at the model
summary(flomodel.01)					# look in more depth
flomodel.02 <- ergm(flomarriage~edges+triangle)		
summary(flomodel.02)
class(flomodel.02)				# R objects possess a class
names(flomodel.02)				# let's look straight at the ERGM obj.
flomodel.02$coef				# the $ allows you to pull an element out from
flomodel.02$mle.lik				# a list
flomodel.02$formula
wealth <- flomarriage %v% 'wealth'
wealth
plot(flomarriage, vertex.cex=wealth/25)

flomodel.03 <- ergm(flomarriage~edges+nodemain('wealth'))
summary(flomodel.03)
data(samplk)						# Let's try a model or two on
ls()								# directed data: Sampson's Monks
samplk3
plot(samplk3)
sampmodel.01 <- ergm(samplk3~edges+mutual)
summary(sampmodel.01)


data(faux.mesa.high)				# Let's try a larger network
mesa <- faux.mesa.high
plot(mesa)						
summary(mesa)					
plot(mesa,vertex.col='Grade')

fauxmodel.01 <- ergm(mesa~edges + nodematch('Grade',diff=T) +
	nodematch('Race',diff=T)
	)

summary(fauxmodel.01)
help('ergm-terms')
help('ergm-terms', chmhelp=FALSE)

flomodel.03.sim <- simulate(flomodel.03,nsim=10)
class(flomodel.03.sim)

names(flomodel.03.sim)

length(flomodel.03.sim$networks)			 

flomodel.03.sim$networks[[1]]		# double brackets pulls an element 
plot(flomodel.03.sim$networks[[1]])
fit <- ergm(flobusiness~edges+kstar(2:3)+triangle, interval=1,burnin=0,seed=2)

mcmc.diagnostics(fit, center=F)


fit <- ergm(flobusiness~edges+kstar(2:3)+triangle)

mcmc.diagnostics(fit, center=F)

data('faux.magnolia.high')
fmh <- faux.magnolia.high
fmh

fit <- ergm(fmh~edges+triangle,seed=1)
mcmc.diagnostics(fit, center=F)

fit <- ergm(fmh~edges+triangle,seed=1,verbose=T)


fit <- ergm(fmh~edges+triangle,seed=1,verbose=T,MCMCsamplesize=20000)
try(mcmc.diagnostics(fit, center=F))

fit <- ergm(fmh~edges+gwesp(0.5,fixed=T),seed=1,verbose=T)
mcmc.diagnostics(fit)

fit <- ergm(fmh~edges+gwesp(0.5,fixed=T)+nodematch('Grade')+nodematch('Race')+
	nodematch('Sex'),seed=1,verbose=T)

pdf('diagnostics.pdf')
mcmc.diagnostics(fit)
dev.off()

fit <- ergm(fmh~edges+gwesp(0.25,fixed=T)+nodematch('Grade')+nodematch('Race')+
	nodematch('Sex'),seed=1,verbose=T)

pdf('diagnostics.pdf')
mcmc.diagnostics(fit)
dev.off()
flomodel.03.gof <- gof(flomodel.03~degree,verbose=T)

flomodel.03.gof
plot(flomodel.03.gof)
fauxmodel.02 <- ergm(fmh~edges)
fauxmodel.02.gof <- gof(fauxmodel.02~distance,nsim=10,verbose=T)
plot(fauxmodel.02.gof)
library(network)                                 #Make sure that network is loaded
data(package='network')                       #List available datasets in network
data(flo)                             #Load a built-in data set; see ?flo for more
flo                                              #Examine the flo adjacency matrix
#Read an adjacency matrix
floadj <- read.table('floadj.txt',header=TRUE)
floadj                                                        #Examine the matrix

#Read a Pajek file
flopaj <- read.paj('flo.paj')
names(flopaj)                 #This is a project file, with networks and other data
names(flopaj$networks)                          #See which networks are in the file
nflo2 <- flopaj$networks[[1]]                            #Extract the marriage data
nflo2                                                   #Examine the network object


nflo <- network(flo, directed=FALSE)          #Create a network object based on flo
nflo                                          #Get a quick description of the data
nempty <- network.initialize(5)              #Create an empty graph with 5 vertices
nempty                                                          #Compare with nflo  


summary(nflo)                                              #Get an overall summary
network.size(nflo)                                      #How large is the network?
network.edgecount(nflo)                               #How many edges are present?
as.sociomatrix(nflo)                                     #Show it as a sociomatrix
nflo[,]                                                      #Another way to do it

plot(nflo,displaylabels=T,boxed.labels=F)                          #Plot with names

#Adding edges
g <- network.initialize(5)                                  #Create an empty graph
g[1,2] <- 1                                               #Add an edge from 1 to 2
g[2,] <- 1                                      #Add edges from 2 to everyone else
g                                                              #Examine the result


#Delete edges
g[3,5] <- 0                                           #Remove the edge from 3 to 5
g                                                                      #It's gone!
g[,] <- 0                                                        #Remove all edges
g                                                             #Now, an empty graph

#Testing adjacency
nflo[9,3]                                                    #Medici to Barbadori?
nflo[9,]                                                        #Entire Medici row
nflo[1:4,5:8]                                                #Subsets are possible
nflo[-9,-9]                                      #Negative numbers _exclude_ nodes

#Setting edge values
m <- matrix(1:16^2, nrow=16, ncol=16)             #Create a matrix of edge 'values'
nflo %e% 'boo' <- m                                       #Value the marriage ties

#Retrieving edge values
list.edge.attributes(nflo)                                   #See what's available
nflo %e% 'boo'                                               #Use the %e% operator

#For more information....
?network.extraction


#Add vertex attributes
nflo %v% 'woo' <- letters[1:16]                               #Letter the families
nflo
list.vertex.attributes(nflo)                           #List all vertex attributes
nflo %v% 'woo'                                      #Retrieve the vertex attribute


library(sna)                                                  #Load the sna library
library(help='sna')                                       #See also this for a list

betweenness(flo)
isolates(nflo)
plot(nflo,displaylabels=T)
kcycle.census(nflo,mode='graph',maxlen=5)

