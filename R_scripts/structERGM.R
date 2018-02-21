#' This project is part of my PhD dissertation project.
#' This section follows from the first section dedicated to the Mathematical Modeling of our co-authorship network (available [HERE](http://#)).
#'
#' The purpose of this section is the use of Exponential Random Graph Modeling (ERGM) as a  Statistical modeling to better understand and explain our malaria co-authorship network.
#' Loading the necessary packages...

#+ setup, include=FALSE,cache=FALSE
setwd('~/R/tb_scripts')
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
load('./Rdata/TBnet.rda')
load('./Rdata/TBnet2.rda')
load('./Rdata/TBnet2.gc.rda')
load('./Rdata/TBauth_data.rda')
# load('./Rdata/edges.rda')
load('./Rdata/TBnet2.gc.rda')
load("Rdata/communitySBM.rda")

# /* Loading a simplified version of the network */
tbnet3 <- read_graph('./graphs/TBnet.graphml', format = 'graphml')
tbnet.w <- read_graph('./graphs/TBnet_weighted.graphml', format = 'graphml')

#' Some data manipulations: Defining Local, Regional and International actors
#+ manip, cahce=FALSE
newCountry<-V(tbnet2)$country %>% str_replace("\n", "")
V(tbnet2)$country<-newCountry
cont<-read.csv("continent.csv", header = T)
continent<-list()
ctnt<-as.vector(cont$Continent)
ctry<-as.vector(cont$Country)
for(i in 1:length(ctnt)){
  continent[[toupper(ctry[i])]]<-toupper(ctnt[i])
}
# V(tbnet2)$continent<-numeric(length(V(tbnet)$country))
newContinent<-c()
for(i in 1:length(V(tbnet2)$country)){
  entry<-continent[[V(tbnet2)$country[i]]]
  if(is.null(entry)){entry<-''}
  newContinent<-c(newContinent, entry)
}
V(tbnet2)$continent<-newContinent

# Assigning collaboration scope
collabType<-c()
for(i in 1:length(V(tbnet2)$continent)){
  collabScope=''
  if(V(tbnet2)$country[i]=='BENIN' && V(tbnet2)$continent[i]=='AFRICA'){
    collabScope<-'NATIONAL'
  } else if(V(tbnet2)$country[i]!='BENIN' && V(tbnet2)$continent[i]=='AFRICA'){
    collabScope<-'REGIONAL'
  } else if(V(tbnet2)$continent[i]!='AFRICA' && V(tbnet2)$continent[i]!='') {
    collabScope<-'INTERNATIONAL'
  } else{
    collabScope<-NA
  }
  collabType<-c(collabType,collabScope)
}
V(tbnet2)$collabType<-collabType


#' Since this processing is going to take a lot of time, we decide to parallelize all the computations:
#+ detect-cores
# /* create a cluster with fixed number of processes */
cores=detectCores() # get number of cores
# cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out
# registerDoParallel(cl)
#'
#+ preproc,cache=FALSE
A <- get.adjacency(tbnet3)
v.attrs <- get.data.frame(tbnet3, what="vertices")
v.attrs$name <- V(tbnet.w)$name
v.attrs$timesCited <- V(tbnet.w)$timesCited
v.attrs$numPub <- V(tbnet.w)$numPub
v.attrs$community <- V(tbnet2)$community
v.attrs$degree <- V(tbnet2)$degree
v.attrs$affiliation <- V(tbnet.w)$place
v.attrs$city <- V(tbnet.w)$affil
v.attrs$country <- V(tbnet.w)$country
v.attrs$collabType <- V(tbnet2)$collabType


tbnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(tbnet3.simple, "timesCited", V(tbnet.w)$timesCited)
network::set.vertex.attribute(tbnet3.simple, "numPub", V(tbnet.w)$numPub)
network::set.vertex.attribute(tbnet3.simple, "community", V(tbnet2)$community)
network::set.vertex.attribute(tbnet3.simple, "degree", V(tbnet2)$degree)
network::set.vertex.attribute(tbnet3.simple, "collabType", V(tbnet2)$collabType)
network::set.vertex.attribute(tbnet3.simple, "communitySBM", communitySBM)

# pairCited <- get.adjacency(tbnet2,attr = 'timesCited')
# nCollab <- as_adjacency_matrix(tbnet2,attr = 'weight')
# 
# # Adding edge attributes
# tbnet3.simple %e% 'pairCited' <- A
# tbnet3.simple %e% 'nCollab' <- nCollab
tbnet3.simple
#'
#+ formula1
# model1.f <- formula(tbnet3.simple ~ edges + triangle + twopath) 
# this model never converges... Instead, fit 
# model1.f <- formula(tbnet3.simple ~ edges + gwesp(1, fixed=TRUE))
model1.f <- formula(tbnet3.simple ~ edges + triangle)
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
# model2.f <- formula(tbnet3.simple ~ edges + gwesp(0, fixed = TRUE) + twopath)
# this model never converges... Instead, fit 
model2.f <- formula(tbnet3.simple ~ edges + gwesp(1, fixed = TRUE) + nodecov("timesCited") + nodecov('degree')
+ nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType")
+ nodefactor("collabType"))
model2.f

summary.statistics(model2.f)
#'
#+ model-fitting1,cache=FALSE
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
# model3.f <- formula(tbnet3.simple ~ edges + triangle + twopath + nodecov("timesCited") + nodecov('degree') +
#                         nodecov("numPub") + nodematch("community") + nodematch("collabType") + nodefactor("community"))
# this model never converges... Instead, fit 
model3.f <- formula(tbnet3.simple ~ edges + gwesp(1, fixed = TRUE) + nodecov("timesCited") + nodecov('degree')
+ nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + nodefactor("communitySBM")
+ nodefactor("collabType"))
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
model4.f <- formula(tbnet3.simple ~ edges + gwesp(0, fixed = TRUE) + twopath + nodecov("timesCited") + nodecov('degree') 
                        + nodecov("numPub") + nodematch("community") + nodematch("collabType") + nodefactor("community") 
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
diag.model1 <- gof(model1, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
Sys.time() - t0
save(diag.model1, file = 'modData2/diag.smodel1.rda')
jpeg(filename = './modData2/gof-smodel1.jpg')
plot(diag.model1)
dev.off()

#'
#+ gof-model2,cache=FALSE
t0 <- Sys.time()
diag.model2 <- gof(model2, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
Sys.time() - t0
save(diag.model2, file = 'modData2/diag.smodel2.rda')
jpeg(filename = './modData2/gof-smodel2.jpg')
plot(diag.model2)
dev.off()

#'
#+ gof-model3,cache=FALSE
t0 <- Sys.time()
diag.model3 <- gof(model3, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
Sys.time() - t0
save(diag.model3, file = 'modData2/diag.smodel3.rda')
jpeg(filename = './modData2/gof-smodel3.jpg')
plot(diag.model3)
dev.off()

#'
#+ gof-model4,cache=FALSE
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
