library(mixer)
library(igraph)

setwd('~/R/R_scripts')

# Loading a simplified version of the network
# mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
mnet3 <- get(load("./Rdata/mnet2.gc.rda"))

# mnet3.sbm <- mixer(as.matrix(get.adjacency(mnet2)), qmin=2, qmax=40)
mnet3.sbm2 <- mixer(as.matrix(get.adjacency(mnet2)), qmin=10, qmax=40)

# mnet3.sbm.output <- getModel(mnet3.sbm)
mnet3.sbm2.output <- getModel(mnet3.sbm2)
# names(mnet3.sbm.output)
names(mnet3.sbm2.output)

# mnet3.sbm.output$q
mnet3.sbm2.output$q

# mnet3.sbm.output$alphas
mnet3.sbm2.output$alphas

# mnet3.sbm.output$Taus[, 1:3]
mnet3.sbm2.output$Taus[, 1:3]

my.ent <- function(x) { -sum(x*log(x, 2)) }
# apply(mnet3.sbm.output$Taus[, 1:3], 2, my.ent)
apply(mnet3.sbm2.output$Taus[, 1:3], 2, my.ent)

# log(mnet3.sbm.output$q, 2)
log(mnet3.sbm2.output$q, 2)

# summary(apply(mnet3.sbm.output$Taus, 2, my.ent))
summary(apply(mnet3.sbm2.output$Taus, 2, my.ent))

# plot(mnet3.sbm, classes=as.factor(V(mnet3)$community))
plot(mnet3.sbm2, classes=as.factor(V(mnet2)$collabType))

save(mnet3.sbm2, file = 'Rdata/modSBM2GC.rda')
#####################
# # Get membership
# community2 <- c()
# A <- mnet3.sbm2.output$Taus
# 
# for(i in 1:length(V(mnet3))){
#   community2 <- c(community2, which(A[,i] == max(A[,i]), arr.ind = TRUE))
# }
# 
# save(community2, file = "./Rdata/communitySBM2.rda")
# V(mnet3)$communitySBM <- community2
# dat<-table(V(mnet2)$collabType,as.factor(community2))
# barplot(dat,legend = rownames(dat), args.legend = list(x = 30),xlab = "Detected communites by SBM", ylab = "Number of authors")
# abline(h=50,col='red',lty=2)


community3 <- c()
B <- mnet3.sbm3.output$Taus

for(i in 1:length(V(mnet3))){
  community3 <- c(community3, which(B[,i] == max(B[,i]), arr.ind = TRUE))
}

dat<-table(V(mnet2)$collabType,as.factor(community3))
barplot(dat,legend = rownames(dat), args.legend = list(x = 30),xlab = "Detected communites by SBM", ylab = "Number of authors")
abline(h=50,col='red',lty=2)
abline(h=25,col='blue',lty=2)
