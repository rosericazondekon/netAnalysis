library(latentnet)
library(igraph)
library(network)
library(stringr)

set.seed(25)

setwd('~/R/R_scripts')
mnet.w <- read_graph('./graphs/CAnet_weight.graphml', format = 'graphml')
mnet2 <- get(load("./Rdata/mnet2.gc.rda"))

#====================================
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

#=====================================
collabType <- V(mnet2)$collabType
collabType <- as.factor(collabType)
levels(collabType)<-1:length(levels(collabType))
collab <- as.numeric(collabType)

V(mnet2)$collab <- collab

V(mnet2)$communitySBM <- get(load("./Rdata/communitySBM2.rda"))
# lazega<-mnet2
A <- get.adjacency(mnet2, sparse=FALSE)

mnet3.simple <- network::as.network(A, directed=FALSE,
                                    matrix.type="a",ignore.eval=FALSE, names.eval="collaborations")

network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet.w)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet.w)$numPub)
network::set.vertex.attribute(mnet3.simple, "communitySBM", V(mnet2)$communitySBM)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)
network::set.vertex.attribute(mnet3.simple, "collab", V(mnet2)$collab)

mnet3.simple

t0 <- Sys.time()
modl1<-ergmm(mnet3.simple~euclidean(d=2), verbose=TRUE,tofit="mcmc",
             control=control.ergmm(burnin=10000,sample.size=500,interval=10,threads=20))
Sys.time() - t0

plot(modl1,what="mcmc", pie=TRUE,labels=TRUE)