cl<-makeCluster(6)
clusterExport(cl=cl, c("list.net2","model.2g"))
dyads<-parLapply(cl, 2:length(list.net2),
              function(t) {
                set.seed(10)
                library(btergm)#;print("...")
                library(network)#;print("...")
                collab <- get.vertex.attribute(list.net2[[t]], "collab")#;print("OK...")
                x <- sample(dim(as.matrix(list.net2[[t]]))[1], 10, replace=F)#;print("OK...")
                mat <- as.matrix(list.net2[[t]])[x,x]#;print("OK...")
                for (i in 1:nrow(mat)) {
                  for (j in 1:ncol(mat)) {
                    if (is.na(collab[i]) == F && is.na(collab[j]) == F && i != j && collab[i] == collab[j]) {
                      c(i, j, t, interpret(model.2g,
                                           type = "tie", i = i, j = j, t = t - 1), collab[i])#;print("probCheck...")
                    }
                  }
                }
              })

dyads <- do.call(rbind, dyads)
dyads <- as.data.frame(dyads)
colnames(dyads) <- c("i", "j", "t", "prob", "collabType")


samplesize <- 10000
results <- list()
for (t in 2:length(friendship)) {
  for (s in 2:3) {
    label <- ifelse(s == 2, paste0("f", t), paste0("m", t))
    d <- dyads[dyads$collab == s & dyads$t == t, ]
    n <- nrow(d)
    means <- sapply(1:samplesize, function(x) {
      samp <- sample(1:n, n, replace = TRUE)
      mean(d[samp, ]$prob)
    })
    results[[label]] <- means
  }
}

# cl<-makeCluster(21)
# registerDoParallel(cl, type="FORK")
# prob<- function(t) {
#   print("Ready...")
#   foreach(t = 2:7, .combine = list, .multicombine = TRUE, .export = c("friendship","model.2g"))  %dopar%  {
#     library(btergm);print("OK...")
#     library(network)
#     collab <- get.vertex.attribute(friendship[[t]], "collab");print("...")
#      x <- sample(1792, 50, replace=F);print("OK...")
#      mat <- as.matrix(friendship[[t]])[x,x];print("OK...")
#      print(t)
#      for (i in 1:nrow(mat)) {
#        for (j in 1:ncol(mat)) {
#          if (is.na(collab[i]) == F && is.na(collab[j]) == F && i != j && collab[i] == collab[j]) {
#           c(i, j, t, interpret(model.2g,
#                                type = "tie", i = i, j = j, t = t - 1), collab[i])
#          }
#        }
#        remove(mat)
#      }
#   }
# }
# dyads<-prob()
stopCluster(cl)
