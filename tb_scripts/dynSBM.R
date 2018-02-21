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
mnet2.96.01 <- subgraph.edges(mnet,which(E(mnet)$year<=2001), delete.vertices = F)
V(mnet2.96.01)$degree<-degree(mnet2.96.01)
# E(mnet2.96.01)$weight <- 1
# mnet2.96.01 <- simplify(mnet2.96.01, edge.attr.comb = list(weight="sum"))

mnet2.02.06 <- subgraph.edges(mnet,which(E(mnet)$year<=2006 & E(mnet)$year>=2002), delete.vertices = F)
V(mnet2.02.06)$degree<-degree(mnet2.02.06)
# E(mnet2.02.06)$weight <- 1
# mnet2.02.06 <- simplify(mnet2.02.06, edge.attr.comb = list(weight="sum"))

mnet2.07.11 <- subgraph.edges(mnet,which(E(mnet)$year<=2011 & E(mnet)$year>=2007), delete.vertices = F)
V(mnet2.07.11)$degree<-degree(mnet2.07.11)
# E(mnet2.07.11)$weight <- 1
# mnet2.07.11 <- simplify(mnet2.07.11, edge.attr.comb = list(weight="sum"))

mnet2.12.16 <- subgraph.edges(mnet,which(E(mnet)$year>=2012), delete.vertices = F)
V(mnet2.12.16)$degree<-degree(mnet2.12.16)
# E(mnet2.12.16)$weight <- 1
# mnet2.12.16 <- simplify(mnet2.12.16, edge.attr.comb = list(weight="sum"))

V(mnet2.96.01)$name <- as.character(1:vcount(mnet2.96.01))
V(mnet2.02.06)$name <- as.character(1:vcount(mnet2.96.01))
V(mnet2.07.11)$name <- as.character(1:vcount(mnet2.96.01))
V(mnet2.12.16)$name <- as.character(1:vcount(mnet2.96.01))
list.mnet <- list(mnet2.96.01, mnet2.02.06, mnet2.07.11, mnet2.12.16)

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


###############################
#UNEVEN GROUPING
##############################

mnet1 <- subgraph.edges(mnet,which(E(mnet)$year<=2014 &
                                     E(mnet)$year>=2012), delete.vertices = F)
mnet2 <- subgraph.edges(mnet,which(E(mnet)$year<=2015 &
                                     E(mnet)$year>=2015), delete.vertices = F)
mnet3 <- subgraph.edges(mnet,which(E(mnet)$year<=2016 &
                                     E(mnet)$year>=2016), delete.vertices = F)
V(mnet1)$name <- as.character(1:vcount(mnet1))
V(mnet2)$name <- as.character(1:vcount(mnet1))
V(mnet3)$name <- as.character(1:vcount(mnet1))
list.mnet2 <- list(mnet1,mnet2,mnet3)


T <- length(list.mnet2)
N <- length(all.vertices.name)
Y2 <- array(dim=c(T,N,N))
for (t in 1:T){
  g.tmp <- graph.empty(n=N, directed=F)
  V(g.tmp)$name <- all.vertices.name
  Y2[t,,] <- as.matrix(get.adjacency(str(g.tmp %u% list.mnet2[[t]]), sparse
                                    = F))
}

absent <- which(sapply(1:N, function(i) sum(Y[,i,])) == 0)
Y2 <- Y2[,-absent,-absent]
N <- dim(Y2)[2]

list.dynsbm2 <- select.dynsbm(Y2, Qmin=1, Qmax=5, edge.type="binary", nstart=1)
dynsbm2 <- list.dynsbm2[[2]]
# Connectivity Plot
connectivity.plot(dynsbm2, Y2)

alluvial.plot(dynsbm2)
