setwd('~/R/R_scripts')
library(igraph)
# library(sand)
mnet <- load(file = './Rdata/sbm_mnet2.rda')
fblog <- get(mnet)
nv <- vcount(fblog)
ncn <- numeric()
A <- get.adjacency(fblog)

t0 <- Sys.time()
for(i in (1:(nv-1))){
  ni <- neighborhood(fblog, 1, i)
  nj <- neighborhood(fblog, 1, (i+1):nv)
  nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
  temp <- unlist(lapply(nbhd.ij, length)) - 
    2*A[i, (i+1):nv]
  ncn <- c(ncn, temp)
}
Sys.time() - t0

# CHUNK 2
par(mfrow=c(1,1))
library(vioplot)
Avec <- A[lower.tri(A)]
vioplot(ncn[Avec==0], ncn[Avec==1], 
        names=c("No Edge", "Edge"))
title(ylab="Number of Common Neighbors")

# CHUNK 3

library(ROCR)
pred <- prediction(ncn, Avec)
perf <- performance(pred, "auc")
slot(perf, "y.values")
# ---
## [[1]]
perf2 <- performance(pred, "tpr", "fpr")
plot(perf2, col='blue', lwd=3)
