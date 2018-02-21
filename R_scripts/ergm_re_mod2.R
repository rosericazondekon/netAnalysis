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

#' Since this processing is going to take a lot of time, we decide to parallelize all the computations:
#+ detect-cores
# /* create a cluster with fixed number of processes */
cores=detectCores() # get number of cores
# /* cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out */
# /* registerDoParallel(cl) */

#+ preproc
A <- get.adjacency(mnet3)
v.attrs <- get.data.frame(mnet3, what="vertices")
v.attrs$name <- V(mnet2)$name
v.attrs$timesCited <- V(mnet2)$timesCited
v.attrs$numPub <- V(mnet2)$numPub
v.attrs$community <- V(mnet2)$community
v.attrs$degree <- V(mnet2)$degree
v.attrs$affiliation <- V(mnet2)$place
v.attrs$city <- V(mnet2)$affil
v.attrs$country <- V(mnet2)$country

mnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet2)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet2)$numPub)
network::set.vertex.attribute(mnet3.simple, "community", V(mnet2)$community)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)

network::set.edge.attribute(mnet3.simple, "pairCited", E(mnet2)$timesCited)
network::set.edge.attribute(mnet3.simple, "year", E(mnet2)$year)
network::set.edge.attribute(mnet3.simple, "nCollab", E(mnet2)$weight)
#+ formula
my.ergm <- formula(mnet3.simple ~ edges + triangle)
summary.statistics(my.ergm)

#+ model-fitting
t0 <- Sys.time()
my.ergm.fit <- ergm(my.ergm, iterations=100,
                    control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=100,
                                         MCMC.interval=500, MCMLE.maxit = 100, MCMC.burnin = 20000, seed=25))
Sys.time() - t0
save(my.ergm.fit, file = './Rdata/my.ergm.fit.rda')

#+ model-summary
anova.ergm(my.ergm.fit)
summary(my.ergm.fit)

#+ gof-myergm
t0 <- Sys.time()
diag <- gof(my.ergm.fit, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
save(diag, file = './Rdata/diag.my.erggm.fit.rda')
Sys.time() - t0

#+ mcmc-diag,,fig.width=20,fig.height=80
mcmc.diagnostics(my.ergm.fit)

#+ viz-diag,fig.width=20,fig.height=60
par(mfrow=c(2, 2))
plot(diag)

