library(igraph)
library(eigenmodel)

## For parallelization
# library(igraph)
# library(foreach)
# library(doParallel)
# library(lubridate)

set.seed(25)

setwd('~/R/R_scripts')

mnet2 <- get(load("./Rdata/sbm_mnet2.rda"))
# mnet2 <- get(load("./Rdata/mnet2.gc.rda"))

#====================================
# newCountry<-V(mnet2)$country %>% str_replace("\n", "")
# V(mnet2)$country<-newCountry
# cont<-read.csv("continent.csv", header = T)
# continent<-list()
# ctnt<-as.vector(cont$Continent)
# ctry<-as.vector(cont$Country)
# for(i in 1:length(ctnt)){
#   continent[[toupper(ctry[i])]]<-toupper(ctnt[i])
# }
# # V(mnet2)$continent<-numeric(length(V(mnet)$country))
# newContinent<-c()
# for(i in 1:length(V(mnet2)$country)){
#   entry<-continent[[V(mnet2)$country[i]]]
#   if(is.null(entry)){entry<-''}
#   newContinent<-c(newContinent, entry)
# }
# V(mnet2)$continent<-newContinent
# 
# # Assigning collaboration scope
# collabType<-c()
# for(i in 1:length(V(mnet2)$continent)){
#   collabScope=''
#   if(V(mnet2)$country[i]=='BENIN' && V(mnet2)$continent[i]=='AFRICA'){
#     collabScope<-'NATIONAL'
#   } else if(V(mnet2)$country[i]!='BENIN' && V(mnet2)$continent[i]=='AFRICA'){
#     collabScope<-'REGIONAL'
#   } else if(V(mnet2)$continent[i]!='AFRICA' && V(mnet2)$continent[i]!='') {
#     collabScope<-'INTERNATIONAL'
#   } else{
#     collabScope<-NA
#   }
#   collabType<-c(collabType,collabScope)
# }
# V(mnet2)$collabType<-collabType

#=====================================

V(mnet2)$communitySBM <- get(load("./Rdata/communitySBM.rda"))
# V(mnet2)$communitySBM <- get(load("./Rdata/communitySBM2.rda"))
# lazega<-mnet2
A <- get.adjacency(mnet2, sparse=FALSE)

# eigenmodel_setup <- function (R = 0, seed = 1, em_env = .GlobalEnv) 
# {
#   n <- dim(Y)[1]
#   uRanks <- 1:length(unique(c(Y[!is.na(Y)])))
#   Ranks <- matrix(match(Y, sort(unique(c(Y)))), n, n)
#   if (!exists("X")) {
#     X <- array(dim = c(n, n, 0))
#   }
#   p <- dim(X)[3]
#   if (p > 0) {
#     x <- NULL
#     for (k in seq(1, dim(X)[3], length = p)) {
#       x <- cbind(x, c((X[, , k])[upper.tri(X[, , k])]))
#     }
#     xtx <- t(x) %*% x
#     tx <- t(x)
#     assign("xtx", xtx, em_env)
#     assign("tx", tx, em_env)
#   }
#   set.seed(seed)
#   RR <- rank(c(Y), ties.method = "random", na.last = "keep")
#   Z <- matrix(qnorm(RR/(sum(!is.na(RR)) + 1)), n, n)
#   Z[is.na(Z)] <- rnorm(sum(is.na(Z)))
#   Z <- Z * upper.tri(Z) + t(Z) * lower.tri(Z)
#   E <- Z
#   b <- rep(0, p)
#   if (p > 0) {
#     b <- lm(E[upper.tri(E)] ~ -1 + t(tx))$coef
#     E <- E - XB(X, b)
#   }
#   E[is.na(E)] <- Z[is.na(E)]
#   tmp <- eigen(E)
#   L <- diag(tmp$val[order(-abs(tmp$val))[seq(1, R, length = R)]]/n, 
#             nrow = R)
#   U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
#                drop = FALSE] * sqrt(n)
#   UL <- list(U = U, L = L)
#   assign("Z", Z, em_env)
#   assign("b", b, em_env)
#   assign("UL", UL, em_env)
#   assign("R", R, em_env)
#   assign("output", NULL, em_env)
#   assign("n", n, em_env)
#   assign("uRanks", uRanks, em_env)
#   assign("Ranks", Ranks, em_env)
#   assign("X", X, em_env)
#   assign("p", p, em_env)
#   assign("pp_b", diag(1/100, nrow = p), em_env)
#   assign("pm_b", rep(0, p), em_env)
#   assign("pp_zq", 1/10, em_env)
#   assign("pp_l", rep(1/100, R), em_env)
#   assign("pp_mu", rep(1/100, R), em_env)
#   assign("pm_mu", rep(0, R), em_env)
#   assign("var_u", rep(1, R), em_env)
#   assign("mean_u", rep(0, R), em_env)
# }


# pp_zq = 1
# Y<-A
# eigenmodel_setup(R=2)

t0 <- Sys.time()
mod1 <- eigenmodel_mcmc(A, R=2, S=1100,burn=1000)
Sys.time() - t0

save(mod1,file = "Rdata/lnm_mod1.rda")

collabType <- V(mnet2)$collabType
collabType <- as.factor(collabType)
levels(collabType)<-1:length(levels(collabType))
collab <- as.numeric(collabType)

V(mnet2)$collab <- collab

v.attr.mnet2 <- get.data.frame(mnet2, what="vertices")
v.attr.mnet2$numPub <- V(mnet2)$numPub
v.attr.mnet2$collabType <- V(mnet2)$collabType
v.attr.mnet2$collab <- V(mnet2)$collab

# CHUNK 24
# same.prac.op <- as.numeric(v.attr.lazega$collabType) %o% as.numeric(v.attr.lazega$collabType)
same.comm.op <- v.attr.mnet2$communitySBM %o% v.attr.mnet2$communitySBM
# same.comm <- matrix(as.numeric(same.comm.op %in% c(1, 4, 9)), 1792, 1792)
same.comm <- matrix(as.numeric(same.comm.op %in% (1:40)^2), 1792, 1792)
same.comm <- array(same.comm,dim=c(1792, 1792, 1))



# CHUNK 25
t0 <- Sys.time()
mod2 <- eigenmodel_mcmc(A, same.comm, R=2,S=1100,burn=1000)
Sys.time() - t0

save(mod2,file = "Rdata/lnm_mod2.rda")

# CHUNK 26
same.aff.op <- v.attr.mnet2$collab %o% v.attr.mnet2$collab
same.aff <- matrix(as.numeric(same.aff.op %in% c(1, 4, 9)), 1792, 1792)
same.aff <- array(same.aff,dim=c(1792, 1792, 1))

t0 <- Sys.time()
mod3 <- eigenmodel_mcmc(A, same.aff, R=2, S=1100, burn=1000)
Sys.time() - t0

save(mod3,file = "Rdata/lnm_mod3.rda")

# CHUNK 27
lat.sp.1 <- eigen(mod1$ULU_postmean)$vec[, 1:2]
lat.sp.2 <- eigen(mod2$ULU_postmean)$vec[, 1:2]
lat.sp.3 <- eigen(mod3$ULU_postmean)$vec[, 1:2]

# CHUNK 28
colbar <- c("red", "dodgerblue", "goldenrod")
colbar[is.na(colbar)] <- "darkgreen"
v.colors <- colbar[V(mnet2)$collab]
v.shapes <- c("circle", "square","csquare")[V(mnet2)$collab]
v.shapes[is.na(v.shapes)] <- "none"
# v.size <- 0.5*sqrt(V(mnet2)$betweenness/1000)
v.size <- 4

par(mfrow=c(3,1))
plot(mnet2, layout=lat.sp.1, vertex.color=v.colors, vertex.label=V(mnet.w)$id, edge.width = 0.07,
     vertex.shape=v.shapes, vertex.size=v.size, vertex.frame.color=NA, vertex.label.size=5)
plot(mnet2, layout=lat.sp.2, vertex.color=v.colors, edge.width = 0.07,
     vertex.shape=v.shapes, vertex.label=NA, vertex.size=v.size, vertex.frame.color=NA)
plot(mnet2, layout=lat.sp.3, vertex.color=v.colors, edge.width = 0.07,
     vertex.shape=v.shapes, vertex.label=NA, vertex.size=v.size, vertex.frame.color=NA)
par(mfrow=c(1,1))

# CHUNK 29
apply(mod1$L_postsamp, 2, mean)
# ---
## [1] 0.2603655 1.0384032
# ---
apply(mod2$L_postsamp, 2, mean)
# ---
## [1]  0.8898394 -0.1156671
# ---
apply(mod3$L_postsamp, 2, mean)
#5 ---
## [1] 0.5970403 0.3112896
# ---

# CHUNK 30
perm.index <- sample(1:1604736)
nfolds <- 5
nmiss <- round(1604736/nfolds,0)
Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))
Avec.pred2 <- numeric(length(Avec))
Avec.pred3 <- numeric(length(Avec))

# CHUNK 31
for(i in seq(1,nfolds)){
  # Index of missing values.
  miss.index <- seq(((i-1) * nmiss + 1), (i*nmiss), 1)
  A.miss.index <- perm.index[miss.index]
  
  # Fill a new Atemp appropriately with NA's.
  Avec.temp <- Avec
  Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
  Avec.temp <- as.numeric(Avec.temp)
  Atemp <- matrix(0, 1792, 1792)
  Atemp[lower.tri(Atemp)] <- Avec.temp
  Atemp <- Atemp + t(Atemp)
  
  # Now fit model and predict.
  Y <- Atemp
  
  model1.fit <- eigenmodel_mcmc(Y, R=2, S=11000, burn=10000)
  model2.fit <- eigenmodel_mcmc(Y, same.comm, R=2, S=11000, burn=10000)
  model3.fit <- eigenmodel_mcmc(Y, same.aff, R=2, S=11000, burn=10000)
  
  model1.pred <- model1.fit$Y_postmean
  model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
  Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
  
  model2.pred <- model2.fit$Y_postmean
  model2.pred.vec <- model2.pred[lower.tri(model2.pred)]
  Avec.pred2[A.miss.index] <- model2.pred.vec[A.miss.index]
  
  model3.pred <- model3.fit$Y_postmean
  model3.pred.vec <- model3.pred[lower.tri(model3.pred)]
  Avec.pred3[A.miss.index] <- model3.pred.vec[A.miss.index]
}

# ## For parallelization -- core detection
# cores=detectCores() # get number of cores=12
# cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out
# registerDoParallel(cl)
# 
# ## PARALLEL COMPUTING...
# start.time <- Sys.time() # Capturing Parallel Computing processing start time
# results <- foreach(1:ntrials, .combine='comb', .multicombine=TRUE,
#                    .init=list(list(),list(),list())) %dopar% {
#   library(igraph)
#   library(eigenmodel)
#   # Index of missing values.
#   miss.index <- seq(((i-1) * nmiss + 1), (i*nmiss), 1)
#   A.miss.index <- perm.index[miss.index]
#   
#   # Fill a new Atemp appropriately with NA's.
#   Avec.temp <- Avec
#   Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
#   Avec.temp <- as.numeric(Avec.temp)
#   Atemp <- matrix(0, 1792, 1792)
#   Atemp[lower.tri(Atemp)] <- Avec.temp
#   Atemp <- Atemp + t(Atemp)
#   
#   # Now fit model and predict.
#   Y <- Atemp
#   
#   model1.fit <- eigenmodel_mcmc(Y, R=2, S=11000, burn=10000)
#   model2.fit <- eigenmodel_mcmc(Y, same.comm, R=2, S=11000, burn=10000)
#   model3.fit <- eigenmodel_mcmc(Y, same.aff, R=2, S=11000, burn=10000)
#   
#   model1.pred <- model1.fit$Y_postmean
#   model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
#   Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
#   
#   model2.pred <- model2.fit$Y_postmean
#   model2.pred.vec <- model2.pred[lower.tri(model2.pred)]
#   Avec.pred2[A.miss.index] <- model2.pred.vec[A.miss.index]
#   
#   model3.pred <- model3.fit$Y_postmean
#   model3.pred.vec <- model3.pred[lower.tri(model3.pred)]
#   Avec.pred3[A.miss.index] <- model3.pred.vec[A.miss.index]
#   
#   list(
#     Avec.pred1,
#     Avec.pred2,
#     Avec.pred3
#   )
# }
# stopCluster(cl) # Stop cluster
# end.time <- Sys.time() # Capturing Parallel Computing processing end time
# time.taken <- end.time - start.time
# time.taken
# 
# Avec.pred1 <- results[[1]]
# Avec.pred2 <- results[[2]]
# Avec.pred3 <- results[[3]]

# CHUNK 32
library(ROCR)
pred1 <- prediction(Avec.pred1, Avec)
perf1 <- performance(pred1, "tpr", "fpr")
plot(perf1, col="blue", lwd=3)

pred2 <- prediction(Avec.pred2, Avec)
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2, col="red", lwd=3)

pred3 <- prediction(Avec.pred3, Avec)
perf3 <- performance(pred3, "tpr", "fpr")
plot(perf3, col="green", lwd=3)

# CHUNK 33
perf1.auc <- performance(pred1, "auc")
slot(perf1.auc, "y.values")

perf2.auc <- performance(pred2, "auc")
slot(perf2.auc, "y.values")

perf3.auc <- performance(pred3, "auc")
slot(perf3.auc, "y.values")
