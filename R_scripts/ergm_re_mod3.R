#' This project is part of my PhD dissertation project.
#' This section follows from the first section dedicated to the Mathematical Modeling of our co-authorship network (available [HERE](http://#)).
#'
#' The purpose of this section is the use of Exponential Random Graph Modeling (ERGM) as a  Statistical modeling to better understand and explain our malaria co-authorship network.
#' Loading the necessary packages...

#+ setup, include=FALSE
setwd('~/R/R_scripts')
library(igraph)
library(ergm)
library(parallel)
library(doParallel)
library(lubridate)
#/* library(Rmpi) */
library(snow)
library(methods)

#' Loading the variables from the last section...

#+ var-load
load('./Rdata/mnet.rda')
load('./Rdata/mnet2.rda')
load('./Rdata/mnet2.gc.rda')
load('./Rdata/auth_data.rda')
load('./Rdata/edges.rda')
load('./Rdata/mnet2.gc.rda')

# /* Loading a simplified version of the network */
mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
mnet.w <- read_graph('./graphs/CAnet_weight.graphml', format = 'graphml')

#' Since this processing is going to take a lot of time, we decide to parallelize all the computations:
#+ detect-cores
# /* create a cluster with fixed number of processes */
cores=detectCores() # get number of cores
# /* cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out */
# /* registerDoParallel(cl) */

#+ preproc
A <- get.adjacency(mnet3)
v.attrs <- get.data.frame(mnet3, what="vertices")
v.attrs$name <- V(mnet.w)$name
v.attrs$timesCited <- V(mnet.w)$timesCited
v.attrs$numPub <- V(mnet.w)$numPub
v.attrs$community <- V(mnet2)$community
v.attrs$degree <- V(mnet2)$degree
v.attrs$affiliation <- V(mnet.w)$place
v.attrs$city <- V(mnet.w)$affil
v.attrs$country <- V(mnet.w)$country

mnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet2)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet2)$numPub)
network::set.vertex.attribute(mnet3.simple, "community", V(mnet2)$community)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)
#+ formula
mnet3.ergm <- formula(mnet3.simple ~ edges + nodecov("timesCited") + 
                        nodecov("numPub") + match("community"))
summary.statistics(mnet3.ergm)

#+ model-fitting
t0 <- Sys.time()
mnet3.ergm.fit <- ergm(mnet3.ergm, iterations=100,
                       control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=100,
                                            MCMC.interval=500, MCMC.burnin = 20000, seed=25, MCMLE.maxit = 100))
Sys.time() - t0
save(mnet3.ergm.fit, file = './Rdata/mnet3.ergm.fit.rda')

#+ model-summary
anova.ergm(mnet3.ergm.fit)
summary(mnet3.ergm.fit)

#+ gof-myergm
t0 <- Sys.time()
diag <- gof(mnet3.ergm.fit, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
save(diag, file = './Rdata/diag.mnet3.ergm.fit.rda')
Sys.time() - t0

#+ mcm_diag,fig.width=20,fig.weight=20
mcmc.diagnostics(mnet3.ergm.fit)

#+ gof-viz,fig.width=20,fig.weight=80
par(mfrow=c(2, 2))
plot(diag)


