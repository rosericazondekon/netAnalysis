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
mnet2.96.97 <- subgraph.edges(mnet,which(E(mnet)$year<=1997), delete.vertices = F)
V(mnet2.96.97)$degree<-degree(mnet2.96.97)
# E(mnet2.96.97)$weight <- 1
# mnet2.96.97 <- simplify(mnet2.96.97, edge.attr.comb = list(weight="sum"))

mnet2.98.99 <- subgraph.edges(mnet,which(E(mnet)$year<=1999 & E(mnet)$year>=1998), delete.vertices = F)
V(mnet2.98.99)$degree<-degree(mnet2.98.99)
# E(mnet2.98.99)$weight <- 1
# mnet2.98.99 <- simplify(mnet2.98.99, edge.attr.comb = list(weight="sum"))

mnet2.00.01 <- subgraph.edges(mnet,which(E(mnet)$year<=2001 & E(mnet)$year>=2000), delete.vertices = F)
V(mnet2.00.01)$degree<-degree(mnet2.00.01)
# E(mnet2.00.01)$weight <- 1
# mnet2.00.01 <- simplify(mnet2.00.01, edge.attr.comb = list(weight="sum"))

mnet2.02.03 <- subgraph.edges(mnet,which(E(mnet)$year<=2003 & E(mnet)$year>=2002), delete.vertices = F)
V(mnet2.02.03)$degree<-degree(mnet2.02.03)
# E(mnet2.02.03)$weight <- 1
# mnet2.02.03 <- simplify(mnet2.02.03, edge.attr.comb = list(weight="sum"))

mnet2.04.05 <- subgraph.edges(mnet,which(E(mnet)$year<=2005 & E(mnet)$year>=2004), delete.vertices = F)
V(mnet2.04.05)$degree<-degree(mnet2.04.05)
# E(mnet2.04.05)$weight <- 1
# mnet2.04.05 <- simplify(mnet2.04.05, edge.attr.comb = list(weight="sum"))

mnet2.06.07 <- subgraph.edges(mnet,which(E(mnet)$year<=2007 & E(mnet)$year>=2006), delete.vertices = F)
V(mnet2.06.07)$degree<-degree(mnet2.06.07)
# E(mnet2.06.07)$weight <- 1
# mnet2.06.07 <- simplify(mnet2.06.07, edge.attr.comb = list(weight="sum"))

mnet2.08.09 <- subgraph.edges(mnet,which(E(mnet)$year<=2009 & E(mnet)$year>=2008), delete.vertices = F)
V(mnet2.08.09)$degree<-degree(mnet2.08.09)
# E(mnet2.08.09)$weight <- 1
# mnet2.08.09 <- simplify(mnet2.08.09, edge.attr.comb = list(weight="sum"))

mnet2.10.11 <- subgraph.edges(mnet,which(E(mnet)$year<=2011 & E(mnet)$year>=2010), delete.vertices = F)
V(mnet2.10.11)$degree<-degree(mnet2.10.11)
# E(mnet2.10.11)$weight <- 1
# mnet2.10.11 <- simplify(mnet2.10.11, edge.attr.comb = list(weight="sum"))

mnet2.12.13 <- subgraph.edges(mnet,which(E(mnet)$year<=2013 & E(mnet)$year>=2012), delete.vertices = F)
V(mnet2.12.13)$degree<-degree(mnet2.12.13)
# E(mnet2.12.13)$weight <- 1
# mnet2.12.13 <- simplify(mnet2.12.13, edge.attr.comb = list(weight="sum"))

mnet2.14.15 <- subgraph.edges(mnet,which(E(mnet)$year<=2015 & E(mnet)$year>=2014), delete.vertices = F)
V(mnet2.14.15)$degree<-degree(mnet2.14.15)
# E(mnet2.14.15)$weight <- 1
# mnet2.14.15 <- simplify(mnet2.14.15, edge.attr.comb = list(weight="sum"))

mnet2.16 <- subgraph.edges(mnet,which(E(mnet)$year>=2016), delete.vertices = F)
V(mnet2.16)$degree<-degree(mnet2.16)
# E(mnet2.16)$weight <- 1
# mnet2.16 <- simplify(mnet2.16, edge.attr.comb = list(weight="sum"))

V(mnet2.96.97)$name <- as.character(1:vcount(mnet2.96.97))
V(mnet2.98.99)$name <- as.character(1:vcount(mnet2.98.99))
V(mnet2.00.01)$name <- as.character(1:vcount(mnet2.00.01))
V(mnet2.02.03)$name <- as.character(1:vcount(mnet2.02.03))
V(mnet2.04.05)$name <- as.character(1:vcount(mnet2.04.05))
V(mnet2.06.07)$name <- as.character(1:vcount(mnet2.06.07))
V(mnet2.08.09)$name <- as.character(1:vcount(mnet2.08.09))
V(mnet2.10.11)$name <- as.character(1:vcount(mnet2.10.11))
V(mnet2.12.13)$name <- as.character(1:vcount(mnet2.12.13))
V(mnet2.14.15)$name <- as.character(1:vcount(mnet2.14.15))
V(mnet2.16)$name <- as.character(1:vcount(mnet2.16))
list.mnet <- list(mnet2.96.97,mnet2.98.99,mnet2.00.01,mnet2.02.03,
                  mnet2.04.05,mnet2.06.07,mnet2.08.09,mnet2.10.11,
                  mnet2.12.13,mnet2.14.15,mnet2.16)

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
list.dynsbm <- select.dynsbm(Y, Qmin=2, Qmax=4, edge.type="continuous", nstart = 1)
Sys.time() - t0
dynsbm <- list.dynsbm[[2]]
# Connectivity Plot
connectivity.plot(dynsbm, Y)

p<-alluvial.plot(dynsbm)
