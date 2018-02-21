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
setwd('~/R/R_scripts')
set.seed(10)
cores=detectCores() - 2 # get number of cores
load('./Rdata/collabType.rda')
load('./Rdata/community2.rda')
load('./Rdata/collab.rda')

# Working with subdivision 1
n1<-read_graph('./graphs/sub1/CAnet_1996-1997.graphml', format = 'graphml')
E(n1)$weight <- 1
net1 <- simplify(n1, edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                           key="concat",subject="concat",year="concat",wosid="concat",
                                           journal="concat",title="concat",doi="concat"))
V(net1)$community2<-community2
V(net1)$collabType<-collabType
V(net1)$collab<-collab
V(net1)$degree<-centralization.degree(net1)$res

n2<-read_graph('./graphs/sub1/CAnet_1998-1999.graphml', format = 'graphml')
E(n2)$weight <- 1
net2 <- simplify(n2, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net2)$community2<-community2
V(net2)$collabType<-collabType
V(net2)$collab<-collab
V(net2)$degree<-centralization.degree(net2)$res

n3<-read_graph('./graphs/sub1/CAnet_2000-2001.graphml', format = 'graphml')
E(n3)$weight <- 1
net3 <- simplify(n3, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net3)$community2<-community2
V(net3)$collabType<-collabType
V(net3)$collab<-collab
V(net3)$degree<-centralization.degree(net3)$res

n4<-read_graph('./graphs/sub1/CAnet_2002-2003.graphml', format = 'graphml')
E(n4)$weight <- 1
net4 <- simplify(n4, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net4)$community2<-community2
V(net4)$collabType<-collabType
V(net4)$collab<-collab
V(net4)$degree<-centralization.degree(net4)$res

n5<-read_graph('./graphs/sub1/CAnet_2004-2005.graphml', format = 'graphml')
E(n5)$weight <- 1
net5 <- simplify(n5, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net5)$community2<-community2
V(net5)$collabType<-collabType
V(net5)$collab<-collab
V(net5)$degree<-centralization.degree(net5)$res

n6<-read_graph('./graphs/sub1/CAnet_2006-2007.graphml', format = 'graphml')
E(n6)$weight <- 1
net6 <- simplify(n6, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net6)$community2<-community2
V(net6)$collabType<-collabType
V(net6)$collab<-collab
V(net6)$degree<-centralization.degree(net6)$res

n7<-read_graph('./graphs/sub1/CAnet_2008-2009.graphml', format = 'graphml')
E(n7)$weight <- 1
net7 <- simplify(n7, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net7)$community2<-community2
V(net7)$collabType<-collabType
V(net7)$collab<-collab
V(net7)$degree<-centralization.degree(net7)$res

n8<-read_graph('./graphs/sub1/CAnet_2010-2011.graphml', format = 'graphml')
E(n8)$weight <- 1
net8 <- simplify(n8, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net8)$community2<-community2
V(net8)$collabType<-collabType
V(net8)$collab<-collab
V(net8)$degree<-centralization.degree(net8)$res

n9<-read_graph('./graphs/sub1/CAnet_2012-2013.graphml', format = 'graphml')
E(n9)$weight <- 1
net9 <- simplify(n9, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net9)$community2<-community2
V(net9)$collabType<-collabType
V(net9)$collab<-collab
V(net9)$degree<-centralization.degree(net9)$res

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

n11<-read_graph('./graphs/sub1/CAnet_2016.graphml', format = 'graphml')
E(n11)$weight <- 1
net11 <- simplify(n11, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))

V(net11)$community2<-community2
V(net11)$collabType<-collabType
V(net11)$collab<-collab
V(net11)$degree<-centralization.degree(net11)$res

net <- list(net1, net2, net3, net4, net5, net6, net7, net8, net9, net10, net11)

list.net<-list()

i=1

while(i<=length(net)){
  V(net[[i]])$name<-as.character(1:vcount(net[[i]]))
  A <- get.adjacency(net[[i]])
  v.attrs <- get.data.frame(net[[i]], what="vertices")
  v.attrs$name <- V(net[[i]])$name
  v.attrs$timesCited <- V(net[[i]])$timesCited
  v.attrs$numPub <- V(net[[i]])$numPub
  v.attrs$community <- V(net[[i]])$community
  v.attrs$degree <- V(net[[i]])$degree
  v.attrs$affiliation <- V(net[[i]])$place
  v.attrs$city <- V(net[[i]])$affil
  v.attrs$country <- V(net[[i]])$country
  v.attrs$collabType <- V(net[[i]])$collabType
  v.attrs$collabType <- V(net[[i]])$collab
  
  net.simple <- network::as.network(as.matrix(A), directed=FALSE)
  network::set.vertex.attribute(net.simple, "timesCited", V(net[[i]])$timesCited)
  network::set.vertex.attribute(net.simple, "numPub", V(net[[i]])$numPub)
  network::set.vertex.attribute(net.simple, "community", V(net[[i]])$community2)
  network::set.vertex.attribute(net.simple, "degree", V(net[[i]])$degree)
  network::set.vertex.attribute(net.simple, "collabType", V(net[[i]])$collabType)
  list.net[[i]]<-net.simple
  i<-i+1
}

# Running Model 1a
t0 <- Sys.time()
model.1a <- btergm(list.net ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +  
                     nodematch("community") + nodematch("collab"), R = 1000, parallel = "snow", ncpus = cores)
Sys.time() - t0
save(model.1a, file = "tergm_sub1/model1a.rda")


# Running Model 1b
t0 <- Sys.time()
model.1b <- btergm(list.net ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") + 
                     nodematch("community") + nodematch("collab") + memory(type = "stability") + 
                     timecov(transform = function(t) t), R = 1000, parallel = "snow", ncpus = cores)
Sys.time() - t0
save(model.1b, file = "tergm_sub1/model1b.rda")

# Model comparison to LaTeX file
texreg(list(model.1a, model.1b), single.row = TRUE,include.nobs = FALSE, file = "tergm_sub1/malaria-coauthNet.tex",
       caption = "Temporal ERGM of Malaria Co-authorship Network.",label = "tab:tergm"
       , custom.model.names =c("Model~1a", "Model~1b"), use.packages = FALSE,
       booktabs = TRUE, dcolumn = TRUE)

# Model diagnostic 1
t0 <- Sys.time()
gof.1a <- gof(model.1a, nsim = 50, statistics = c(esp, geodesic, deg), parallel ="multicore", ncpus = cores)
Sys.time() - t0
save(gof.1a, file = "tergm_sub1/gof1a.rda")

# Model diagnostic 2
t0 <- Sys.time()
gof.1b <- gof(model.1b, nsim = 50, statistics = c(esp, geodesic, deg), parallel ="multicore", ncpus = cores)
Sys.time() - t0
save(gof.1b, file = "tergm_sub1/gof1b.rda")

