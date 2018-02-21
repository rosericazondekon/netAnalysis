# detachAllPackages <- function() {
#   basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
#   package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#   package.list <- setdiff(package.list,basic.packages)
#   if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
# }
# detachAllPackages()

library(igraph)
# library(parallel)
# library(doParallel)
# library(lubridate)
# #/* library(Rmpi) */
# library(snow)
# library(methods)
# library(dplyr)
# library(stringr)
# library(sna)
# library(statnet)
# library(texreg)
# library(xergm)
# library(network)
library(RSiena)
setwd('~/R/R_scripts')
set.seed(10)
cores=detectCores() - 2 # get number of cores
# cl <- parallel::makeCluster(cores, type = "FORK")
# parallel::stopCluster(cl=NULL)
load('./Rdata/collabType.rda')
load('./Rdata/community2.rda')
load('./Rdata/collab.rda')

# Working with subdivision 1
n1<-read_graph('./graphs/sub2/CAnet_1996-2006.graphml', format = 'graphml')
E(n1)$weight <- 1
net1 <- simplify(n1, edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                           key="concat",subject="concat",year="concat",wosid="concat",
                                           journal="concat",title="concat",doi="concat"))
V(net1)$community2<-community2
V(net1)$collabType<-collabType
V(net1)$collab<-collab
V(net1)$degree<-centralization.degree(net1)$res

n2<-read_graph('./graphs/sub2/CAnet_2007-2009.graphml', format = 'graphml')
E(n2)$weight <- 1
net2 <- simplify(n2, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net2)$community2<-community2
V(net2)$collabType<-collabType
V(net2)$collab<-collab
V(net2)$degree<-centralization.degree(net2)$res

n3<-read_graph('./graphs/sub2/CAnet_2010-2011.graphml', format = 'graphml')
E(n3)$weight <- 1
net3 <- simplify(n3, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net3)$community2<-community2
V(net3)$collabType<-collabType
V(net3)$collab<-collab
V(net3)$degree<-centralization.degree(net3)$res

n4<-read_graph('./graphs/sub2/CAnet_2012-2013.graphml', format = 'graphml')
E(n4)$weight <- 1
net4 <- simplify(n4, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net4)$community2<-community2
V(net4)$collabType<-collabType
V(net4)$collab<-collab
V(net4)$degree<-centralization.degree(net4)$res

n5<-read_graph('./graphs/sub2/CAnet_2014.graphml', format = 'graphml')
E(n5)$weight <- 1
net5 <- simplify(n5, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net5)$community2<-community2
V(net5)$collabType<-collabType
V(net5)$collab<-collab
V(net5)$degree<-centralization.degree(net5)$res

n6<-read_graph('./graphs/sub2/CAnet_2015.graphml', format = 'graphml')
E(n6)$weight <- 1
net6 <- simplify(n6, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net6)$community2<-community2
V(net6)$collabType<-collabType
V(net6)$collab<-collab
V(net6)$degree<-centralization.degree(net6)$res

n7<-read_graph('./graphs/sub2/CAnet_2016.graphml', format = 'graphml')
E(n7)$weight <- 1
net7 <- simplify(n7, 
                 edge.attr.comb = list(weight="sum",timesCited="sum",numPub="sum",
                                       key="concat",subject="concat",year="concat",wosid="concat",
                                       journal="concat",title="concat",doi="concat"))
V(net7)$community2<-community2
V(net7)$collabType<-collabType
V(net7)$collab<-collab
V(net7)$degree<-centralization.degree(net7)$res



lnet2 <- list(net1, net2, net3, net4, net5, net6, net7)


# Make a function that symmetrizes the network
sym.min<-function(x){
  tx<-t(x)
  return(pmin(x[],tx[]))
}


# pubs<-data.frame(Date=as.Date(character()),File=character(), User=character(), stringsAsFactors=FALSE)
# tcited<-data.frame(Date=as.Date(character()),File=character(), User=character(), stringsAsFactors=FALSE) 
pubs<-matrix(, nrow = 1792, ncol = 7)
tCited<-matrix(, nrow = 1792, ncol = 7)
for(i in 1:length(lnet2)){
  pubs[,i]<-V(lnet2[[i]])$numPub
  tCited[,i]<-V(lnet2[[i]])$timesCited
  assign(paste0("mnet",i),sym.min(as.matrix(get.adjacency(lnet2[[i]]))))
}

# Prepare data for RSiena.
mnet <- sienaNet( array( c( mnet1,mnet2,mnet3,mnet4,mnet5,mnet6,mnet7 ),
                               dim = c( 1792, 1792, 7 ) ) )# create dependent variable

community <- coCovar( community2 )# create constant covariate
collaboration <- coCovar( collab )# create constant covariate

timesCited <- varCovar( tCited )# create time varying covariate
numPub <- varCovar( pubs )# create time varying covariate

mydata <- sienaDataCreate( mnet, community, collaboration, timesCited, numPub )# define data
myeff <- getEffects( mydata )# create effects structure
print01Report( mydata, myeff, modelname = 'malaria_net' )# siena01 for reports

# Which effects are available for the symmetric network?
effectsDocumentation(myeff)

# Run the model in the various different specifications,

#myalgo  <- sienaAlgorithmCreate(projname = "malaria_net")
myalgo1 <- sienaAlgorithmCreate(projname = "malaria_net", modelType = 2)
myalgo2 <- sienaAlgorithmCreate(projname = "malaria_net", modelType = 3)
myalgo3 <- sienaAlgorithmCreate(projname = "malaria_net", modelType = 4)
myalgo4 <- sienaAlgorithmCreate(projname = "malaria_net", modelType = 5, firstg = 0.02)
myalgo5 <- sienaAlgorithmCreate(projname = "malaria_net", modelType = 6)


(ans <- siena07( myalgo1, data = mydata, effects = myeff, useCluster=TRUE, 
                 nbrNodes=48,clusterString=nbrNodes,parallelTesting=F,clusterType="FORK"))
(ans2 <- siena07( myalgo2, data = mydata, effects = myeff, useCluster=TRUE, 
                  nbrNodes=48,clusterString=48,parallelTesting=F,clusterType="FORK"))
(ans3 <- siena07( myalgo3, data = mydata, effects = myeff, useCluster=TRUE, 
                  nbrNodes=48,clusterString=48,parallelTesting=F,clusterType="FORK"))
# In my case, convergence here was not satisfactory; therefore a second try:
(ans3 <- siena07( myalgo3, data = mydata, effects = myeff, prevAns=ans3, useCluster=TRUE, 
                  nbrNodes=48,clusterString=48,parallelTesting=TRUE,clusterType="FORK"))
(ans4 <- siena07( myalgo4, data = mydata, effects = myeff, useCluster=TRUE, 
                  nbrNodes=48,clusterString=48,parallelTesting=F,clusterType="FORK"))
(ans5 <- siena07( myalgo5, data = mydata, effects = myeff, verbose=TRUE, useCluster=TRUE, 
                  nbrNodes=48,clusterString=48,parallelTesting=TRUE,clusterType="FORK"))


