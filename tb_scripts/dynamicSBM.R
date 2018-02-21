library(dynsbm)
library(igraph)

setwd('~/R/R_scripts')

# Loading a simplified version of the network
mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
mnet2.96.01 <- read_graph('./graphs/CAnet_96_01s.graphml', format = 'graphml')
mnet2.02.06 <- read_graph('./graphs/CAnet_02_06s.graphml', format = 'graphml')
mnet2.07.11 <- read_graph('./graphs/CAnet_07_11s.graphml', format = 'graphml')
mnet2.12.16 <- read_graph('./graphs/CAnet_12_16s.graphml', format = 'graphml')

list.mnet <- list(mnet2.96.01, mnet2.02.06, mnet2.07.11, mnet2.12.16)

# converting to R array
all.vertices.name <- unique(sort(sapply(list.g, function(g) V(mnet3)$id)))
T <- length(list.g)
N <- length(all.vertices.name)
Y <- array(dim=c(T,N,N))
for (t in 1:T){
  g.tmp <- graph.empty(n=N, directed=F)
  V(g.tmp)$name <- all.vertices.name
  Y[t,,] <- as.matrix(get.adjacency(union(g.tmp,list.mnet[[t]])))
}

# dynSBM
list.dynsbm <- select.dynsbm(Y, Qmin=2, Qmax=3, edge.type="binary")
dynsbm <- list.dynsbm[[2]]
connectivity.plot(dynsbm, Y)
alluvial.plot(dynsbm)