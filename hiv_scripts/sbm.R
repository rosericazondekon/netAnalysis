library(mixer)
library(igraph)

setwd('~/R/hiv_scripts')

# Loading a simplified version of the network
hivnet <- read_graph('./graphs/HIVnet_weighted.graphml', format = 'graphml')

newCountry<-V(hivnet)$country %>% str_replace("\n", "")
V(hivnet)$country<-newCountry
cont<-read.csv("continent.csv", header = T)
continent<-list()
ctnt<-as.vector(cont$Continent)
ctry<-as.vector(cont$Country)
for(i in 1:length(ctnt)){
  continent[[toupper(ctry[i])]]<-toupper(ctnt[i])
}
newContinent<-c()
for(i in 1:length(V(hivnet)$country)){
  entry<-continent[[V(hivnet)$country[i]]]
  if(is.null(entry)){entry<-''}
  newContinent<-c(newContinent, entry)
}
V(hivnet)$continent <- newContinent

# Assigning collaboration scope
collabType<-c()
for(i in 1:length(V(hivnet)$continent)){
  collabScope=''
  if(V(hivnet)$country[i]=='BENIN' && V(hivnet)$continent[i]=='AFRICA'){
    collabScope<-'NATIONAL'
  } else if(V(hivnet)$country[i]!='BENIN' && V(hivnet)$continent[i]=='AFRICA'){
    collabScope<-'REGIONAL'
  } else if(V(hivnet)$continent[i]!='AFRICA' && V(hivnet)$continent[i]!='') {
    collabScope<-'INTERNATIONAL'
  } else{
    collabScope<-NA
  }
  collabType<-c(collabType,collabScope)
}
V(hivnet)$collabType<-collabType


hivnet.sbm <- mixer(as.matrix(get.adjacency(hivnet)),
                   qmin=2, qmax=30)

hivnet.sbm.output <- getModel(hivnet.sbm)
names(hivnet.sbm.output)

hivnet.sbm.output$q

hivnet.sbm.output$alphas

hivnet.sbm.output$Taus[, 1:3]

my.ent <- function(x) { -sum(x*log(x, 2)) }
apply(hivnet.sbm.output$Taus[, 1:3], 2, my.ent)

log(hivnet.sbm.output$q, 2)

summary(apply(hivnet.sbm.output$Taus, 2, my.ent))

# plot(hivnet.sbm, classes=as.factor(V(hivnet)$community))
plot(hivnet.sbm, classes=as.factor(V(hivnet)$collabType), quantile.val=1)


communitySBM <- c()
A <- hivnet.sbm.output$Taus

for(i in 1:length(V(hivnet))){
  communitySBM <- c(communitySBM, which(A[,i] == max(A[,i]), arr.ind = TRUE))
}

save(communitySBM, file = "./Rdata/communitySBM.rda")

dat<-table(V(hivnet)$collabType,as.factor(communitySBM))
barplot(dat,legend = rownames(dat), args.legend = list(x = 30),xlab = "Detected communites by SBM", ylab = "Number of authors")
abline(h=20,col='red',lty=2)
abline(h=10,col='blue',lty=2)

################################
# hivnet.sbm2 <- mixer(as.matrix(get.adjacency(hivnet)),
#                     qmin=2, qmax=6)
# jpeg(filename = "figure/hiv_sbm_g6.jpeg",width = 1000,height = 420)
plot(hivnet.sbm, q = 6, frame = c(2,4), classes=as.factor(V(hivnet)$collabType), quantile.val=0.5)
# dev.off()
