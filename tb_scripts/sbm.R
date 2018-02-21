library(mixer)
library(igraph)

setwd('~/R/tb_scripts')

# Loading a simplified version of the network
tbnet <- read_graph('./graphs/TBnet_weighted.graphml', format = 'graphml')

newCountry<-V(tbnet)$country %>% str_replace("\n", "")
V(tbnet)$country<-newCountry
cont<-read.csv("continent.csv", header = T)
continent<-list()
ctnt<-as.vector(cont$Continent)
ctry<-as.vector(cont$Country)
for(i in 1:length(ctnt)){
  continent[[toupper(ctry[i])]]<-toupper(ctnt[i])
}
newContinent<-c()
for(i in 1:length(V(tbnet)$country)){
  entry<-continent[[V(tbnet)$country[i]]]
  if(is.null(entry)){entry<-''}
  newContinent<-c(newContinent, entry)
}
V(tbnet)$continent <- newContinent

# Assigning collaboration scope
collabType<-c()
for(i in 1:length(V(tbnet)$continent)){
  collabScope=''
  if(V(tbnet)$country[i]=='BENIN' && V(tbnet)$continent[i]=='AFRICA'){
    collabScope<-'NATIONAL'
  } else if(V(tbnet)$country[i]!='BENIN' && V(tbnet)$continent[i]=='AFRICA'){
    collabScope<-'REGIONAL'
  } else if(V(tbnet)$continent[i]!='AFRICA' && V(tbnet)$continent[i]!='') {
    collabScope<-'INTERNATIONAL'
  } else{
    collabScope<-NA
  }
  collabType<-c(collabType,collabScope)
}
V(tbnet)$collabType<-collabType


tbnet.sbm <- mixer(as.matrix(get.adjacency(tbnet)),
                   qmin=2, qmax=15)

tbnet.sbm.output <- getModel(tbnet.sbm)
names(tbnet.sbm.output)

tbnet.sbm.output$q

tbnet.sbm.output$alphas

tbnet.sbm.output$Taus[, 1:3]

my.ent <- function(x) { -sum(x*log(x, 2)) }
apply(tbnet.sbm.output$Taus[, 1:3], 2, my.ent)

log(tbnet.sbm.output$q, 2)

summary(apply(tbnet.sbm.output$Taus, 2, my.ent))

# plot(tbnet.sbm, classes=as.factor(V(tbnet)$community))
plot(tbnet.sbm, classes=as.factor(V(tbnet)$collabType))


communitySBM <- c()
A <- tbnet.sbm.output$Taus

for(i in 1:length(V(tbnet))){
  communitySBM <- c(communitySBM, which(A[,i] == max(A[,i]), arr.ind = TRUE))
}

save(communitySBM,file = './Rdata/communitySBM.rda')

dat<-table(V(tbnet)$collabType,as.factor(communitySBM))
barplot(dat,legend = rownames(dat), args.legend = list(x = 30),xlab = "Detected communites by SBM", ylab = "Number of authors")
abline(h=10,col='red',lty=2)
abline(h=5,col='blue',lty=2)