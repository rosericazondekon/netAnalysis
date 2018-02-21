t0 <- Sys.time()
gof.mod2g <- gof(model.2g, nsim = 200, target = list.net2[6],
                      formula = list.net2[2:5] ~ edges + nodecov("timesCited") + nodecov('degree') + nodecov("numPub") +
                        nodematch("community") + nodematch("collabType") + nodefactor("collabType") +
                        memory(type = "stability") + timecov(transform = function(t) t),
                      coef = coef(model.2g), statistics = c(esp, dsp, geodesic, deg, triad.undirected,rocpr, 
                                                            walktrap.modularity, walktrap.roc, edgebetweenness.modularity),
                      parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0
# esp, dsp, geodesic, deg, triad.undirected, rocpr

t0 <- Sys.time()
gof.mod2g <- gof(model.2g, nsim = 100, statistics = c(esp, dsp, geodesic, triad.undirected, rocpr,
                                                      walktrap.modularity, walktrap.roc, edgebetweenness.modularity),
                 parallel ="multicore", ncpus = cores, cl=cl)
Sys.time() - t0


##################################
friendship <- list.net2
set.seed(1)
t0 <- Sys.time()
dyads <- list()
for (t in 2:length(friendship)) {
  print(t)
  collab <- get.vertex.attribute(friendship[[t]], "collabType")
  #x <- sample(dim(as.matrix(friendship[[t]]))[1], 100, replace=F)
  mat <- as.matrix(friendship[[t]])#[x,x]
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      if (is.na(collab[i]) == F && is.na(collab[j]) == F && i != j && collab[i] == collab[j]) {
        dyads[[length(dyads) + 1]] <- c(i, j, t, interpret(model.2g,
                                                           type = "tie", i = i, j = j, t = t - 1), collab[i])
      }
    }
  }
}
Sys.time() - t0
dyads <- do.call(rbind, dyads)
dyads <- as.data.frame(dyads)
colnames(dyads) <- c("i", "j", "t", "prob", "collabType")
save(dyads,file = "./tergm_sub2/dyads.rda")
# write.csv(dyads,file = "./tergm_sub2/dyads.csv")

library(readr)
# dyads2 <- read_csv("~/R/hiv_scripts/dyads.csv", col_types = cols(prob = col_number()))
dyads2 <- read_csv("~/R/R_scripts/tergm_sub2/dyads.csv", col_types = cols(collabType = col_number(),
                                                                          i = col_number(), j = col_number(),
                                                                          prob = col_number(), t = col_number()))

save(dyads,file = "./tergm_sub2/dyads.rda")

samplesize <- 10000
results <- list()
for (t in 2:length(friendship)) {
  print(t)
  for (s in c("INTERNATIONAL","NATIONAL","REGIONAL")) {
  # for (s in 1:3) {
    label <- ifelse(s == "INTERNATIONAL", paste0("INTERNATIONAL", t), 
                    ifelse(s=="NATIONAL", paste0("NATIONAL", t), paste0("REGIONAL",t)))
    d <- dyads2[dyads2$collabType == s & dyads2$t == t, ]
    # means<-mean(d$prob)
    n <- nrow(d)
    means <- sapply(1:samplesize, function(x) {
    samp <- sample(1:n, n, replace = T)
    mean(d[samp, ]$prob)
    })
  results[[label]] <- means
  }
}

quantiles <- sapply(results, function(x) {
    return(c(quantile(x, 0.025,na.rm = F), quantile(x, 0.5,na.rm = F), quantile(x, 0.975,na.rm = F)))
})

# colType<-rep(c("INTERNATIONAL","NATIONAL","REGIONAL"),5)
# yr<-c(rep("1996-2001",3),rep("2002-2008",3),rep("2009-2010",3),rep("2011-2012",3),rep("2013-2014",3))
# data<-as.data.frame(cbind(yr,colType,"prob"=quantiles[1,],"lb"=quantiles[2,],"ub"=quantiles[3,]))
# rownames(data)<-NULL
# # colnames(dyads) <- c("yr", "colType", "prob", "lb", "ub")
# write.csv(data,file = "barplotData.csv")
# 
# data <- read_csv("~/R/hiv_scripts/barplotData.csv",col_types = cols(lb = col_number(), prob = col_number(),ub = col_number()))


# ggplot(data, aes(x=yr, y=prob, fill=colType)) + 
#   geom_bar(position=position_dodge(), stat="identity") +
#   geom_errorbar(aes(ymin=ub-prob, ymax=lb+prob)
#                 ,width=.2,                    # Width of the error bars
#                 position=position_dodge(.9))


#####################
library("gplots")
# library(ggplot2)
par(mfrow=c(1,1))
# pdf("interpret.pdf")
barplot2(quantiles[2, ], col = c("lightpink", "cornflowerblue","darkgoldenrod1"),
         # col = c("lightpink", "cornflowerblue","tomato"),
         plot.ci = TRUE, ci.l = quantiles[1, ], ci.u = quantiles[3, ],
         # names = c("", "1996-2001","","", "2002-2008","","", "2009-2010","","", "2011-2012","","", "2013-2014",""),
         # names = c("", "1996-2006","","", "2007-2009","","", "2010-2011","","", "2012-2013","","", "2014","","", "2015",""),
         names = c("", "1996-2008","","", "2009-2011","","", "2012-2013","","", "2014-2015","","", "2016",""),#tb
         # names = c("1996-2001","2002-2008","2009-2010","2011-2012","2013-2014"), 
         # space = c(0.0,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0),#hiv
         # space = c(0.0,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0,0),#malaria
         space = c(0.0,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0,0.6,0.0,0.0),#tb
         xlab = "Years",ylim = c(0,0.15),
         # width=0.5,
         ylab = "Median edge probability (with 95 percent CI)", ci.lwd = 2,
         main = "Same collaboration type probabilities over time")
legend("top",fill=c("lightpink", "cornflowerblue","darkgoldenrod1"), legend=c("INTERNATIONAL", "NATIONAL","REGIONAL"))
# dev.off()