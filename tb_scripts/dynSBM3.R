library(mixer)
library(igraph)
library(doParallel)
library(parallel)
library(dynsbm)

setwd('~/R/R_scripts')
set.seed(123)

# Loading a simplified version of the network
mnet <- read_graph('./graphs/CAnet.graphml', format = 'graphml')

# Aggregating network by time periods
mnet2.96.06 <- subgraph.edges(mnet,which(E(mnet)$year<=2006), delete.vertices = F)
V(mnet2.96.06)$degree<-degree(mnet2.96.06)
# E(mnet2.96.06)$weight <- 1
# mnet2.96.06 <- simplify(mnet2.96.06, edge.attr.comb = list(weight="sum"))

mnet2.07.09 <- subgraph.edges(mnet,which(E(mnet)$year<=2009 & E(mnet)$year>=2007), delete.vertices = F)
V(mnet2.07.09)$degree<-degree(mnet2.07.09)
# E(mnet2.07.09)$weight <- 1
# mnet2.07.09 <- simplify(mnet2.07.09, edge.attr.comb = list(weight="sum"))

mnet2.10.11 <- subgraph.edges(mnet,which(E(mnet)$year<=2011 & E(mnet)$year>=2010), delete.vertices = F)
V(mnet2.10.11)$degree<-degree(mnet2.10.11)
# E(mnet2.10.11)$weight <- 1
# mnet2.10.11 <- simplify(mnet2.10.11, edge.attr.comb = list(weight="sum"))

mnet2.12.13 <- subgraph.edges(mnet,which(E(mnet)$year<=2013 & E(mnet)$year>=2012), delete.vertices = F)
V(mnet2.12.13)$degree<-degree(mnet2.12.13)
# E(mnet2.12.13)$weight <- 1
# mnet2.12.13 <- simplify(mnet2.12.13, edge.attr.comb = list(weight="sum"))

mnet2.14 <- subgraph.edges(mnet,which(E(mnet)$year==2014), delete.vertices = F)
V(mnet2.14)$degree<-degree(mnet2.14)
# E(mnet2.14)$weight <- 1
# mnet2.14 <- simplify(mnet2.14, edge.attr.comb = list(weight="sum"))

mnet2.15 <- subgraph.edges(mnet,which(E(mnet)$year==2015), delete.vertices = F)
V(mnet2.15)$degree<-degree(mnet2.15)
# E(mnet2.15)$weight <- 1
# mnet2.15 <- simplify(mnet2.15, edge.attr.comb = list(weight="sum"))

mnet2.16 <- subgraph.edges(mnet,which(E(mnet)$year>=2016), delete.vertices = F)
V(mnet2.16)$degree<-degree(mnet2.16)
# E(mnet2.16)$weight <- 1
# mnet2.16 <- simplify(mnet2.16, edge.attr.comb = list(weight="sum"))

V(mnet2.96.06)$name <- as.character(1:vcount(mnet2.96.06))
V(mnet2.07.09)$name <- as.character(1:vcount(mnet2.07.09))
V(mnet2.10.11)$name <- as.character(1:vcount(mnet2.10.11))
V(mnet2.12.13)$name <- as.character(1:vcount(mnet2.12.13))
V(mnet2.14)$name <- as.character(1:vcount(mnet2.14))
V(mnet2.15)$name <- as.character(1:vcount(mnet2.15))
V(mnet2.16)$name <- as.character(1:vcount(mnet2.16))
list.mnet <- list(mnet2.96.06,mnet2.07.09,mnet2.10.11,mnet2.12.13,mnet2.14,mnet2.15,mnet2.16)

all.vertices.name <- unique(sort(sapply(list.mnet, function(g) V(g)$name)))
T <- length(list.mnet)
N <- length(all.vertices.name)
Y <- array(dim=c(T,N,N))
for (t in 1:T){
  g.tmp <- graph.empty(n=N, directed=F)
  V(g.tmp)$name <- all.vertices.name
  Y[t,,] <- as.matrix(get.adjacency(str(g.tmp %u% list.mnet[[t]]), sparse=F))
}

# dynSBM
t0 <- Sys.time()
# list.dynsbm <- select.dynsbm(Y, Qmin=2, Qmax=4, edge.type="binary", nstart = 1)
list.dynsbm <- select.dynsbm(Y, Qmin=1, Qmax=50, edge.type="continuous", nstart = 1)
Sys.time() - t0
dynsbm <- list.dynsbm[[2]]
# Connectivity Plot
connectivity.plot(dynsbm, Y)

p<-alluvial.plot(dynsbm)
