#' This project is part of my PhD dissertation project.
#' This section follows from the first section dedicated to the Mathematical Modeling of our co-authorship network (available [HERE](http://#)).
#'
#' The purpose of this section is the use of Exponential Random Graph Modeling (ERGM) as a  Statistical modeling to better understand and explain our malaria co-authorship network.
#' Loading the necessary packages...

#+ setup, include=FALSE,cache=FALSE
setwd('~/R/R_scripts')
library(igraph)
library(ergm)
library(parallel)
library(doParallel)
library(lubridate)
#/* library(Rmpi) */
library(snow)
library(methods)
library(dplyr)
library(stringr)

#' Loading the variables from the last section...

#+ var-load,cache=FALSE
load('./Rdata/mnet.rda')
load('./Rdata/mnet2.rda')
load('./Rdata/mnet2.gc.rda')
load('./Rdata/auth_data.rda')
load('./Rdata/edges.rda')
load('./Rdata/mnet2.gc.rda')
load('./Rdata/communitySBM.rda')

# /* Loading a simplified version of the network */
mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
mnet.w <- read_graph('./graphs/CAnet_weight.graphml', format = 'graphml')

#' Some data manipulations: Defining Local, Regional and International actors
#+ manip, cahce=FALSE
newCountry<-V(mnet2)$country %>% str_replace("\n", "")
V(mnet2)$country<-newCountry
cont<-read.csv("continent.csv", header = T)
continent<-list()
ctnt<-as.vector(cont$Continent)
ctry<-as.vector(cont$Country)
for(i in 1:length(ctnt)){
  continent[[toupper(ctry[i])]]<-toupper(ctnt[i])
}
# V(mnet2)$continent<-numeric(length(V(mnet)$country))
newContinent<-c()
for(i in 1:length(V(mnet2)$country)){
  entry<-continent[[V(mnet2)$country[i]]]
  if(is.null(entry)){entry<-''}
  newContinent<-c(newContinent, entry)
}
V(mnet2)$continent<-newContinent

# Assigning collaboration scope
collabType<-c()
for(i in 1:length(V(mnet2)$continent)){
  collabScope=''
  if(V(mnet2)$country[i]=='BENIN' && V(mnet2)$continent[i]=='AFRICA'){
    collabScope<-'NATIONAL'
  } else if(V(mnet2)$country[i]!='BENIN' && V(mnet2)$continent[i]=='AFRICA'){
    collabScope<-'REGIONAL'
  } else if(V(mnet2)$continent[i]!='AFRICA' && V(mnet2)$continent[i]!='') {
    collabScope<-'INTERNATIONAL'
  } else{
    collabScope<-NA
  }
  collabType<-c(collabType,collabScope)
}
V(mnet2)$collabType<-collabType


#' Since this processing is going to take a lot of time, we decide to parallelize all the computations:
#+ detect-cores
# /* create a cluster with fixed number of processes */
cores=detectCores() # get number of cores
# cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out
# registerDoParallel(cl)
#'
#+ preproc,cache=FALSE
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
v.attrs$collabType <- V(mnet2)$collabType


mnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet.w)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet.w)$numPub)
network::set.vertex.attribute(mnet3.simple, "community", V(mnet2)$community)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)
network::set.vertex.attribute(mnet3.simple, "collabType", V(mnet2)$collabType)
network::set.vertex.attribute(mnet3.simple, "communitySBM", community2)

# pairCited <- get.adjacency(mnet2,attr = 'timesCited')
# nCollab <- as_adjacency_matrix(mnet2,attr = 'weight')
# 
# # Adding edge attributes
# mnet3.simple %e% 'pairCited' <- A
# mnet3.simple %e% 'nCollab' <- nCollab
mnet3.simple
#'
#+ formula1
#model1.f <- formula(mnet3.simple ~ edges + triangle + twopath)
model1.f <- formula(mnet3.simple ~ edges + nodecov("timesCited") + nodecov('degree') +
                        nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType"))
model1.f

summary.statistics(model1.f)
#'
#+ model-fitting1,cache=FALSE
t0 <- Sys.time()
model1 <- ergm(model1.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
Sys.time() - t0

#'
#+ model-summary1,cache=FALSE
anova.ergm(model1)
summary(model1)
#'
#+ save-mod1,cache=FALSE
save(model1, file = 'modData2/smodel1.rda')

#######################################################################################
#######################################################################################
#######################################################################################
#'
#+ formula2
#model2.f <- formula(mnet3.simple ~ edges + gwesp(0, fixed = TRUE) + twopath)
model2.f <- formula(mnet3.simple ~ edges + nodecov("timesCited") + nodecov('degree') +
                        nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + nodefactor("collabType"))
model2.f

summary.statistics(model2.f)
#'
#+ model-fitting2,cache=FALSE
t0 <- Sys.time()
model2 <- ergm(model2.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
Sys.time() - t0

#'
#+ model-summary2,cache=FALSE
anova.ergm(model2)
summary(model2)
#'
#+ save-mod2,cache=FALSE
save(model2, file = 'modData2/smodel2.rda')

#######################################################################################
#######################################################################################
#######################################################################################

#'
#+ formula3
model3.f <- formula(mnet3.simple ~ edges + triangle + twopath)
model3.f

summary.statistics(model3.f)
#'
#+ model-fitting3,cache=FALSE
t0 <- Sys.time()
model3 <- ergm(model3.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
Sys.time() - t0

#'
#+ model-summary3,cache=FALSE
anova.ergm(model3)
summary(model3)

#'
#+ mod3-mcmc_diag,cache=FALSE
mcmc.model3 <- mcmc.diagnostics(model3)
mcmc.model3
save(mcmc.model3, file = "modData/mcmc-smodel3.rda")
#'
#+ save-mod3,cache=FALSE
save(model3, file = 'modData2/smodel3.rda')

#######################################################################################
#######################################################################################
#######################################################################################
#'
#+ formula4
model4.f <- formula(mnet3.simple ~ edges + triangle + twopath + nodecov("timesCited") + nodecov('degree') 
                        + nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") 
                        + nodefactor("collabType"))
model4.f

summary.statistics(model4.f)
#'
#+ model-fitting4,cache=FALSE
t0 <- Sys.time()
model4 <- ergm(model4.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
Sys.time() - t0

#'
#+ model-summary4,cache=FALSE
anova.ergm(model4)
summary(model4)

#'
#+ mod4-mcmc_diag,cache=FALSE
mcmc.model4 <- mcmc.diagnostics(model4)
mcmc.model4
save(mcmc.model4, file = "modData/mcmc-smodel4.rda")
#'
#+ save-mod4,cache=FALSE
save(model4, file = 'modData2/smodel4.rda')

#'
#+ gof-model1,cache=FALSE
t0 <- Sys.time()
diag.model2 <- gof(model2, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
Sys.time() - t0
save(diag.model2, file = 'modData2/diag.smodel2.rda')
jpeg(filename = './modData2/gof-smodel2.jpg')
plot(diag.model2)
dev.off()

#'
#+ gof-model2,cache=FALSE
t0 <- Sys.time()
diag.model4 <- gof(model4, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
Sys.time() - t0
save(diag.model4, file = 'modData2/diag.smodel4.rda')
jpeg(filename = './modData2/gof-smodel4.jpg')
plot(diag.model4)
dev.off()

#######################################################################################
#################################### END ##############################################
#######################################################################################
