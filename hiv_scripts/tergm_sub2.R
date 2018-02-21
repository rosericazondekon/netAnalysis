library(igraph)
library(parallel)
library(doParallel)
library(lubridate)
#/* library(Rmpi) */
library(snow)
library(methods)
library(dplyr)
library(stringr)
library(sna)
library(statnet)
library(texreg)
library(xergm)
library(network)
setwd('~/R/hiv_scripts')
set.seed(10)
cores=detectCores() - 2 # get number of cores
cl <- parallel::makeCluster(cores, type = "FORK")
collabType<-c(get(load('./Rdata/collabType.rda')),NA)
community2<-c(get(load('./Rdata/communitySBM.rda')),26)
# load('./Rdata/collab.rda')
collab<-collabType
collab[which(collab=="INTERNATIONAL")]<-1
collab[which(collab=="NATIONAL")]<-2
collab[which(collab=="REGIONAL")]<-3
collab[which(is.na(collab))]<-4
collab<-as.numeric(collab)

# Working with subdivision 1
n1<-read_graph('./graphs/sub2/CAnet_1996-2001.graphml', format = 'graphml')
E(n1)$weight <- 1
net1 <- simplify(n1, edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                           key="concat",subject="concat",year="concat",wosid="concat",
                                           journal="concat",title="concat",doi="concat"))
V(net1)$community2<-community2
V(net1)$collabType<-collabType
V(net1)$collab<-collab
V(net1)$degree<-centralization.degree(net1)$res

n2<-read_graph('./graphs/sub2/CAnet_2002-2008.graphml', format = 'graphml')
E(n2)$weight <- 1
net2 <- simplify(n2, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net2)$community2<-community2
V(net2)$collabType<-collabType
V(net2)$collab<-collab
V(net2)$degree<-centralization.degree(net2)$res

n3<-read_graph('./graphs/sub2/CAnet_2009-2010.graphml', format = 'graphml')
E(n3)$weight <- 1
net3 <- simplify(n3, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net3)$community2<-community2
V(net3)$collabType<-collabType
V(net3)$collab<-collab
V(net3)$degree<-centralization.degree(net3)$res

n4<-read_graph('./graphs/sub2/CAnet_2011-2012.graphml', format = 'graphml')
E(n4)$weight <- 1
net4 <- simplify(n4, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net4)$community2<-community2
V(net4)$collabType<-collabType
V(net4)$collab<-collab
V(net4)$degree<-centralization.degree(net4)$res

n5<-read_graph('./graphs/sub2/CAnet_2013-2014.graphml', format = 'graphml')
E(n5)$weight <- 1
net5 <- simplify(n5, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net5)$community2<-community2
V(net5)$collabType<-collabType
V(net5)$collab<-collab
V(net5)$degree<-centralization.degree(net5)$res

n6<-read_graph('./graphs/sub2/CAnet_2015-2016.graphml', format = 'graphml')
E(n6)$weight <- 1
net6 <- simplify(n6, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net6)$community2<-community2
V(net6)$collabType<-collabType
V(net6)$collab<-collab
V(net6)$degree<-centralization.degree(net6)$res

# n7<-read_graph('./graphs/sub2/CAnet_2016.graphml', format = 'graphml')
# E(n7)$weight <- 1
# net7 <- simplify(n7, 
#                  edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
#                                        key="concat",subject="concat",year="concat",wosid="concat",
#                                        journal="concat",title="concat",doi="concat"))
# V(net7)$community2<-community2
# V(net7)$collabType<-collabType
# V(net7)$collab<-collab
# V(net7)$degree<-centralization.degree(net7)$res



net2 <- list(net1, net2, net3, net4, net5, net6)

list.net2<-list()

i=1
# N=1792
# T=length(net2)
# weightNet<-array(dim=c(N,N,T))
while(i<=length(net2)){
  V(net2[[i]])$name<-as.character(1:vcount(net2[[i]]))
  nt<-net2[[i]]
  # nt<-igraph::delete.vertices(net2[[i]], centralization.degree(net2[[i]])$res==0)
  A <- get.adjacency(nt)
  v.attrs <- get.data.frame(nt, what="vertices")
  v.attrs$name <- V(nt)$name
  v.attrs$timesCited <- V(nt)$timesCited
  v.attrs$numPub <- V(nt)$numPub
  v.attrs$community <- V(nt)$community
  v.attrs$degree <- V(nt)$degree
  v.attrs$affiliation <- V(nt)$place
  v.attrs$city <- V(nt)$affil
  v.attrs$country <- V(nt)$country
  v.attrs$collabType <- V(nt)$collabType
  v.attrs$collab <- V(nt)$collab
  
  net.simple <- network::as.network(as.matrix(A), directed=FALSE)
  network::set.vertex.attribute(net.simple, "timesCited", V(nt)$timesCited)
  network::set.vertex.attribute(net.simple, "numPub", V(nt)$numPub)
  network::set.vertex.attribute(net.simple, "community", V(nt)$community2)
  network::set.vertex.attribute(net.simple, "degree", V(nt)$degree)
  network::set.vertex.attribute(net.simple, "collabType", V(nt)$collabType)
  network::set.vertex.attribute(net.simple, "collab", V(nt)$collab)
  
  # network::set.edge.value(net.simple, "weight", get.adjacency(nt, attr = "weight"))
  network::set.network.attribute(net.simple, "weight", as.matrix(get.adjacency(nt, attr = "weight")))
  network::set.network.attribute(net.simple, "timesCited", as.matrix(get.adjacency(nt, attr = "timesCited")))
  
  # weightNet[,,i]<-as.matrix(get.adjacency(nt, sparse=T))
  
  list.net2[[i]]<-net.simple
  i<-i+1
}

# save(list.net2,file = "./tergm_sub2/listnet2.rda")

# # Running Model 2a
# t0 <- Sys.time()
# model.2a <- btergm(list.net2 ~ edges + triangle + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +  
#                      nodematch("community") + nodematch("collab") + nodefactor("collabType")
#                    , R = 5000, parallel = "snow", ncpus = cores)
# Sys.time() - t0
# save(model.2a, file = "tergm_sub2/model2a.rda")
# 
# 
# # Running Model 2b
# t0 <- Sys.time()
# model.2b <- btergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
#                      nodematch("community") + nodematch("collab") + memory(type = "stability") + nodefactor("collabType")
#                      , R = 5000, parallel = "snow", ncpus = cores)
# Sys.time() - t0
# save(model.2b, file = "tergm_sub2/model2b.rda")
# 
# 
# # Running Model 2c
# t0 <- Sys.time()
# model.2c <- btergm(list.net2 ~ edges + triangle + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
#                      nodematch("community") + nodematch("collab") + memory(type = "stability") + nodefactor("collabType") + 
#                      timecov(transform = function(t) t), R = 5000, parallel = "snow", ncpus = cores, cl=cl)
# Sys.time() - t0
# save(model.2c, file = "tergm_sub2/model2c.rda")
# 
# # Running Model 2d
# t0 <- Sys.time()
# model.2d <- btergm(list.net2 ~ edges + triangle + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
#                      nodematch("community") + nodematch("collab") + delrecip + nodefactor("collabType") + 
#                      timecov(transform = function(t) t), R = 5000, parallel = "snow", ncpus = cores, cl=cl)
# # model.2d <- mtergm(list.net2 ~ edges + triangle + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
# #                      nodematch("community") + nodematch("collab") + memory(type = "stability") 
# #                      , R = 1000, ncpus = cores,
# #                    control = control.ergm(parallel.type="PSOCK",MCMC.samplesize = 5000, MCMC.interval = 2000, parallel = "snow"))
# Sys.time() - t0
# save(model.2d, file = "tergm_sub2/model2d.rda")

## Bootstrapped TERGM (BTERGM) estimates are too large and warrant no reliability
## I decide to use MLE or MCMC estimates using the mtergm() function

# Null TERGM model
t0 <- Sys.time()
model.2d <- mtergm(list.net2 ~ edges #+ nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                   #nodematch("community") + nodematch("collab") + memory(type = "stability")
                   , R = 1000, ncpus = cores, cl=cl)
#,control = control.ergm(parallel.type="PSOCK",MCMC.samplesize = 5000, MCMC.interval = 2000, parallel = cores))
Sys.time() - t0

# Model with covariates only ===> uses MLE estimates
t0 <- Sys.time()
model.2e <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collab") #+ memory(type = "stability")
                   , R = 1000, ncpus = cores, cl=cl)
#,control = control.ergm(parallel.type="PSOCK",MCMC.samplesize = 5000, MCMC.interval = 2000, parallel = cores))
Sys.time() - t0
# save(model.2e, file = "tergm_sub2/model2e.rda")

# Model with covariates and nodefactor term
t0 <- Sys.time()
model.2f <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collab") + nodefactor("collabType")
                   # + memory(type = "stability") + timecov(transform = function(t) t)
                   , R = 1000, ncpus = cores, cl=cl,control = control.ergm(MCMLE.maxit = 1000,
                                                                           MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2f, file = "tergm_sub2/model2f.rda")


# Model with covariates, linear trends and dyadic stability term ===> uses MCMC estimates, increased computation time
t0 <- Sys.time()
model.2g <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                     memory(type = "stability") + timecov(transform = function(t) t)
                   , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
                                                                                            MCMC.burnin = 20000)) # delrecip does not have any effct on mdel improvement
Sys.time() - t0
# save(model.2g, file = "tergm_sub2/model2g.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000


# Model with covariates and polynomial time function
t0 <- Sys.time()
model.2h <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
                     memory(type = "autoregression") + timecov(transform = function(t) t)
                   , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
                                                                                            MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2h, file = "tergm_sub2/model2h.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000


# Model with covariates and polynomial time function
t0 <- Sys.time()
model.2hb <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
                     memory(type = "loss") + timecov(transform = function(t) t)
                   , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
                                                                                            MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2h, file = "tergm_sub2/model2h.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000

t0 <- Sys.time()
model.2hc <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                      nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
                      memory(type = "innovation") + timecov(transform = function(t) t)
                    , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
                                                                                             MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2h, file = "tergm_sub2/model2h.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000

# save(model.2h, file = "tergm_sub2/model2h.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000

# t0 <- Sys.time()
# model.2hc <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                       nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
#                       memory(type = "innovation") + timecov(transform = function(t) t)
#                     , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
#                                                                                              MCMC.burnin = 20000))
# Sys.time() - t0
# save(model.2h, file = "tergm_sub2/model2h.rda")
#, MCMC.samplesize = 5000, MCMC.interval = 2000

# t0 <- Sys.time()
# model.2hd <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                      nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
#                      memory(type="stability") + memory(type = "loss") + timecov(transform = function(t) t)
#                    , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
#                                                                                             MCMC.burnin = 20000))
# Sys.time() - t0

# t0 <- Sys.time()
# model.2i <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                      nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
#                      memory(type = "loss") + timecov(transform = function(t) (1 + t + t^2))
#                    , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
#                                                                                             MCMC.burnin = 20000))
# Sys.time() - t0
# 
# 
# t0 <- Sys.time()
# model.2j <- mtergm(list.net2 ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                      nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
#                      memory(type = "innovation") + timecov(transform = function(t) (1 + t + t^2))
#                    , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
#                                                                                             MCMC.burnin = 20000))
# Sys.time() - t0


######################################
# New Network Data (2 last years)
######################################

n10<-read_graph('./graphs/sub1/CAnet_2014-2015.graphml', format = 'graphml')
E(n10)$weight <- 1
net10 <- simplify(n10, 
                  edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                        key="concat",subject="concat",year="concat",wosid="concat",
                                        journal="concat",title="concat",doi="concat"))
V(net10)$community2<-community2
V(net10)$collabType<-collabType
V(net10)$collab<-collab
V(net10)$degree<-centralization.degree(net10)$res

##
V(net10)$name<-as.character(1:vcount(net10))
A <- get.adjacency(net10)
v.attrs <- get.data.frame(net10, what="vertices")
v.attrs$name <- V(net10)$name
v.attrs$timesCited <- V(net10)$timesCited
v.attrs$numPub <- V(net10)$numPub
v.attrs$community <- V(net10)$community
v.attrs$degree <- V(net10)$degree
v.attrs$affiliation <- V(net10)$place
v.attrs$city <- V(net10)$affil
v.attrs$country <- V(net10)$country
v.attrs$collabType <- V(net10)$collabType
v.attrs$collabType <- V(net10)$collab

net.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(net.simple, "timesCited", V(net10)$timesCited)
network::set.vertex.attribute(net.simple, "numPub", V(net10)$numPub)
network::set.vertex.attribute(net.simple, "community", V(net10)$community2)
network::set.vertex.attribute(net.simple, "degree", V(net10)$degree)
network::set.vertex.attribute(net.simple, "collabType", V(net10)$collabType)

######################################

############# MODEL 2F: BEST MODEL in terms of AIC, BIC and Log Likelihood ######################

# # Goodness-of-fit
# t0 <- Sys.time()
# gof.mod2e <- gof(model.2e, nsim = 20, target = list.net2[[7]],
#                  formula = list.net2[1:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                    nodematch("community") + nodematch("collab") + memory(type = "stability") + 
#                    timecov(transform = function(t) t),
#                  coef = coef(model.2e), statistics = c(esp, dsp, geodesic, deg, triad.undirected, rocpr),
#                  parallel ="multicore", ncpus = cores, cl=cl)
# Sys.time() - t0
# 
# save(gof.mod2e, file = "tergm_sub2/gof.mod2e.rda")
# plot(gof.mod2e, roc.rgraph = TRUE, pr.rgraph = TRUE)
# 
# gof.mod2e#AUC=0.7998484
# 
# 
# #GOF model.2f
# t0 <- Sys.time()
# gof.mod2f <- gof(model.2f, nsim = 20, target = list.net2[[7]],
#                  formula = list.net2[1:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
#                    nodematch("community") + nodematch("collab") + nodefactor("collabType") +
#                    memory(type = "stability") + timecov(transform = function(t) t),
#                  coef = coef(model.2f), statistics = c(esp, dsp, geodesic, deg, triad.undirected, rocpr),
#                  parallel ="multicore", ncpus = cores, cl=cl)
# Sys.time() - t0
# save(gof.mod2f, file = "tergm_sub2/gof.mod2f.rda")
# 
# plot(gof.mod2f, roc.rgraph = TRUE, pr.rgraph = TRUE)
# gof.mod2f #AUC=0.8001797

## TIE PREDICTION 1 ########################
# We recompute Model 2g, omitting the last time step in the estimation
t0 <- Sys.time()
model.2gpred <- mtergm(list.net2[1:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                         nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                         memory(type = "stability") + timecov(transform = function(t) t)
                       , R = 1000, ncpus = cores, cl=cl,control = control.ergm(MCMLE.maxit = 1000,
                                                                               MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2fpred, file = "tergm_sub2/model2fpred.rda")

t0 <- Sys.time()
gof.mod2gpred <- gof(model.2gpred, nsim = 1000, target = list.net2[[7]],
                     formula = list.net2[5:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                       nodematch("community") + nodematch("collab") + nodefactor("collabType") +
                       memory(type = "stability") + timecov(transform = function(t) t),
                     coef = coef(model.2gpred), statistics = c(rocpr),
                     parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
# esp, dsp, geodesic, deg, triad.undirected, rocpr
save(gof.mod2gpred, file = "tergm_sub2/gof.mod2fpred.rda")

plot(gof.mod2gpred, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2gpred
# gof.mod2fpred[[6]]$auc.roc # AUC.ROC=0.9865334
# gof.mod2fpred[[6]]$auc.pr # AUC.ROC=0.9510345


# Predicting net.simple
t0 <- Sys.time()
gof.mod2gpred2 <- gof(model.2gpred, nsim = 350, target = net.simple,
                      formula = list.net2[1:5] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                        nodematch("community") + nodematch("collab") + nodefactor("collabType") +
                        memory(type = "stability") + timecov(transform = function(t) t),
                      coef = coef(model.2gpred), statistics = c(rocpr),
                      parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
# esp, dsp, geodesic, deg, triad.undirected, rocpr
save(gof.mod2gpred2, file = "tergm_sub2/gof.mod2fpred2.rda")

plot(gof.mod2gpred2, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2gpred2



## TIE PREDICTION 2 ########################
t0 <- Sys.time()
model.2gpred2 <- mtergm(list.net2[2:5] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                          nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                          memory(type = "stability") + timecov(transform = function(t) t)
                        , R = 1000, ncpus = cores, cl=cl,control = control.ergm(MCMLE.maxit = 1000,
                                                                                MCMC.burnin = 20000))
Sys.time() - t0
# save(model.2fpred, file = "tergm_sub2/model2fpred.rda")

t0 <- Sys.time()
gof.mod2fpred2 <- gof(model.2gpred2, nsim = 100, target = net.simple,
                      formula = list.net2[2:5] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                        nodematch("community") + nodematch("collab") + nodefactor("collabType") +
                        memory(type = "stability") + timecov(transform = function(t) t),
                      coef = coef(model.2gpred2), 
                      statistics = c(fastgreedy.roc,fastgreedy.pr,walktrap.roc,walktrap.pr,
                                     edgebetweenness.roc,edgebetweenness.pr,rocpr),
                      parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
# esp, dsp, geodesic, deg, triad.undirected,
save(gof.mod2fpred2, file = "tergm_sub2/gof.mod2fpred2.rda")
plot(gof.mod2fpred2, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2fpred2
# ROC model ROC random  PR model  PR random
# 1 0.9016556  0.4993918 0.8088972 0.05374174



##################################
#GOF model.2h
t0 <- Sys.time()
gof.mod2h <- gof(model.2h, nsim = 20, target = list.net2[[7]],
                 formula = list.net2[1:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                   nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                   memory(type = "stability") + timecov(transform = function(t) (t + t^2)),
                 coef = coef(model.2h), statistics = c(esp, dsp, geodesic, deg, triad.undirected, rocpr),
                 parallel ="multicore", ncpus = cores, cl=cl) # No substantial improvement with polynomial time function
Sys.time() - t0
save(gof.mod2h, file = "tergm_sub2/gof.mod2h.rda")

plot(gof.mod2h, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2h #AUC=0.7996671


# Alittle more....
t0 <- Sys.time()
gof.mod2i <- gof(model.2f, nsim = 25, target = list.net2[[6]],
                 formula = list.net2[2:5] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                   nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                   memory(type = "stability") + timecov(transform = function(t) t),
                 coef = coef(model.2f), statistics = c(esp, dsp, geodesic, deg, triad.undirected, rocpr),
                 parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
plot(gof.mod2f, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2i



######################### FINAL MODEL is model.2f #############################
# Now let's simulate more to improve ROC curve
# We start with 50 simulations/time = 300 simulations for the first 6 time steps
# last time stamp is used for tie prediction
t0 <- Sys.time()
gof.mod2j <- gof(model.2f, nsim = 50, target = list.net2[[7]],
                 formula = list.net2[1:6] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                   nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                   memory(type = "stability") + timecov(transform = function(t) t),
                 coef = coef(model.2f), statistics = c(rocpr),
                 parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
plot(gof.mod2f, roc.rgraph = TRUE, pr.rgraph = TRUE)
gof.mod2j


stopCluster(cl)

# STERGM
# t0 <- Sys.time()
# model.2f <- stergm(list.net2,
#                    formation=~edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
#                      nodematch("community") + nodematch("collab") + nodefactor("collabType"),
#                    dissolution=~edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
#                      nodematch("community") + nodematch("collab") + nodefactor("collabType"),
#                    estimate="CMLE",control=control.stergm(parallel=46, parallel.type="PSOCK",
#                    CMLE.MCMC.burnin=20000))
# Sys.time() - t0
#, times=1:7

# Model table to LaTeX table
texreg(list(model.2d, model.2e, model.2f, model.2g, model.2h, model.2hb, model.2hc), single.row = TRUE,include.nobs = FALSE, 
       file = "tergm_sub2/hiv_tergm-coauthNet5.tex",
       caption = "Temporal ERGM of HIV/AIDS Co-authorship Network.",label = "tab:hiv_tergm"
       , custom.model.names =c("Model~2a", "Model~2b", "Model~2c", "Model~2d~(stab)", "Model~2e~(auto)", "Model~2f~(loss)", "Model~2g~(innov)"), use.packages = FALSE,
       booktabs = TRUE, dcolumn = TRUE)


# Model table to LaTeX table
texreg(list(model.2d, model.2e, model.2f, model.2g), single.row = TRUE,include.nobs = FALSE, 
       file = "tergm_sub2/hiv_tergm-coauthNet.tex",
       caption = "Temporal ERGM of HIV/AIDS Co-authorship Network.",label = "tab:hiv_tergm"
       , custom.model.names =c("Model~2a", "Model~2b", "Model~2c", "Model~2d"), use.packages = FALSE,
       booktabs = TRUE, dcolumn = TRUE)

####################
texreg(list(model.2f, model.2g, model.2h, model.2i, model.2j), single.row = TRUE,include.nobs = FALSE, 
       file = "tergm_sub2/malaria-coauthNet4.tex",
       caption = "Temporal ERGM of Malaria Co-authorship Network.",label = "tab:tergm"
       , custom.model.names =c("Model~2f", "Model~2g(stab)", "Model~2h(autoreg)", "Model~2i(loss)", "Model~2j(innov)"), use.packages = FALSE,
       booktabs = TRUE, dcolumn = TRUE)

# Model diagnostic 2a
t0 <- Sys.time()
gof.2a <- gof(model.2a, nsim = 50, statistics = c(esp, geodesic, deg), parallel ="multicore", ncpus = cores)
Sys.time() - t0
save(gof.2a, file = "tergm_sub2/gof2a.rda")


# Model diagnostic 2b
t0 <- Sys.time()
gof.2b <- gof(model.2b, nsim = 50, statistics = c(esp, geodesic, deg), parallel ="multicore", ncpus = cores)
Sys.time() - t0
save(gof.2b, file = "tergm_sub2/gof2b.rda")

########################################################
list.hivnet<-list.net2
t0 <- Sys.time()
hivtergm <- mtergm(list.hivnet ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                     nodematch("community") + nodematch("collab") + nodefactor("collabType") + 
                     memory(type = "stability") + timecov(transform = function(t) t)
                   , R = 1000, ncpus = cores, cl=cl, iterations=1000,control = control.ergm(MCMLE.maxit = 1000,
                                                                                            MCMC.burnin = 20000)) # delrecip does not have any effct on mdel improvement
Sys.time() - t0

#Saving tergm model....
saveRDS(hivtergm,"tergm_sub2/hivtergm.rds")

saveRDS(list.hivnet,"tergm_sub2/listhivnet.rds")
