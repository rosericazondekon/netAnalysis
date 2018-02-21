This project is part of my PhD dissertation project.
This section follows from the first section dedicated to the Mathematical Modeling of our co-authorship network (available [HERE](http://#)).

The purpose of this section is the use of Exponential Random Graph Modeling (ERGM) as a  Statistical modeling to better understand and explain our malaria co-authorship network.
Loading the necessary packages...



Loading the variables from the last section...


```r
load('./Rdata/mnet.rda')
load('./Rdata/mnet2.rda')
load('./Rdata/mnet2.gc.rda')
load('./Rdata/auth_data.rda')
load('./Rdata/edges.rda')
load('./Rdata/mnet2.gc.rda')
load('./Rdata/communitySBM.rda')

mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
mnet.w <- read_graph('./graphs/CAnet_weight.graphml', format = 'graphml')
```

Some data manipulations: Defining Local, Regional and International actors


```r
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
```

Since this processing is going to take a lot of time, we decide to parallelize all the computations:


```r
cores=detectCores() # get number of cores
# cl <- makeCluster(cores[1]-2) #I decide to run leave 2 cores out
# registerDoParallel(cl)
```




```r
A <- get.adjacency(mnet3)
v.attrs <- get.data.frame(mnet3, what="vertices")
v.attrs$name <- V(mnet.w)$name
v.attrs$timesCited <- V(mnet.w)$timesCited
v.attrs$numPub <- V(mnet.w)$numPub
v.attrs$community <- V(mnet2)$community
v.attrs$degree <- V(mnet2)$degree
v.attrs$affiliation <- V(mnet.w)$place
v.attrs$city <- V(mnet.w)$affil
v.attrs$country <- V(mnet.w)$country
v.attrs$collabType <- V(mnet2)$collabType


mnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet.w)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet.w)$numPub)
network::set.vertex.attribute(mnet3.simple, "community", V(mnet2)$community)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)
network::set.vertex.attribute(mnet3.simple, "collabType", V(mnet2)$collabType)
network::set.vertex.attribute(mnet3.simple, "communitySBM", community2)

# pairCited <- get.adjacency(mnet2,attr = 'timesCited')
# nCollab <- as_adjacency_matrix(mnet2,attr = 'weight')
# 
# # Adding edge attributes
# mnet3.simple %e% 'pairCited' <- A
# mnet3.simple %e% 'nCollab' <- nCollab
mnet3.simple
```

```
##  Network attributes:
##   vertices = 1792 
##   directed = FALSE 
##   hyper = FALSE 
##   loops = FALSE 
##   multiple = FALSE 
##   bipartite = FALSE 
##   total edges= 95707 
##     missing edges= 0 
##     non-missing edges= 95707 
## 
##  Vertex attribute names: 
##     collabType community communitySBM degree numPub timesCited vertex.names 
## 
##  Edge attribute names not shown
```




```r
#model1.f <- formula(mnet3.simple ~ edges + triangle + twopath)
model1.f <- formula(mnet3.simple ~ edges + nodecov("timesCited") + nodecov('degree') +
                        nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType"))
model1.f
```

```
## mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType")
```

```r
summary.statistics(model1.f)
```

```
##                  edges     nodecov.timesCited         nodecov.degree 
##                  95707               33146116               56641204 
##         nodecov.numPub nodematch.communitySBM   nodematch.collabType 
##                 495605                  75410                  58334
```




```r
t0 <- Sys.time()
model1 <- ergm(model1.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## Evaluating log-likelihood at the estimate.
```

```r
Sys.time() - t0
```

```
## Time difference of 9.187222 secs
```




```r
anova.ergm(model1)
```

```
## Analysis of Variance Table
## 
## Model 1: mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType")
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                   1604736          0                 
## Model 1:  6  -220952   1604730     220952    < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(model1)
```

```
## 
## ==========================
## Summary of model fit
## ==========================
## 
## Formula:   mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType")
## 
## Iterations:  8 out of 1000 
## 
## Monte Carlo MLE Results:
##                          Estimate Std. Error MCMC % p-value    
## edges                  -7.988e+00  1.866e-02      0  <1e-04 ***
## nodecov.timesCited     -1.459e-02  1.069e-04      0  <1e-04 ***
## nodecov.degree          1.549e-02  6.681e-05      0  <1e-04 ***
## nodecov.numPub          1.293e-01  8.258e-04      0  <1e-04 ***
## nodematch.communitySBM  5.576e+00  1.775e-02      0  <1e-04 ***
## nodematch.collabType    4.617e-01  1.242e-02      0  <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2224636  on 1604736  degrees of freedom
##  Residual Deviance:  220952  on 1604730  degrees of freedom
##  
## AIC: 220964    BIC: 221038    (Smaller is better.)
```




```r
save(model1, file = 'modData2/smodel1.rda')

#######################################################################################
#######################################################################################
#######################################################################################
```




```r
#model2.f <- formula(mnet3.simple ~ edges + gwesp(0, fixed = TRUE) + twopath)
model2.f <- formula(mnet3.simple ~ edges + nodecov("timesCited") + nodecov('degree') +
                        nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + nodefactor("collabType"))
model2.f
```

```
## mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + 
##     nodefactor("collabType")
```

```r
summary.statistics(model2.f)
```

```
##                          edges             nodecov.timesCited 
##                          95707                       33146116 
##                 nodecov.degree                 nodecov.numPub 
##                       56641204                         495605 
##         nodematch.communitySBM           nodematch.collabType 
##                          75410                          58334 
## nodefactor.collabType.NATIONAL nodefactor.collabType.REGIONAL 
##                           6124                          30809
```




```r
t0 <- Sys.time()
model2 <- ergm(model2.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## Evaluating log-likelihood at the estimate.
```

```r
Sys.time() - t0
```

```
## Time difference of 12.10573 secs
```




```r
anova.ergm(model2)
```

```
## Analysis of Variance Table
## 
## Model 1: mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + 
##     nodefactor("collabType")
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                   1604736          0                 
## Model 1:  8  -217010   1604728     217010    < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(model2)
```

```
## 
## ==========================
## Summary of model fit
## ==========================
## 
## Formula:   mnet3.simple ~ edges + nodecov("timesCited") + nodecov("degree") + 
##     nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") + 
##     nodefactor("collabType")
## 
## Iterations:  9 out of 1000 
## 
## Monte Carlo MLE Results:
##                                  Estimate Std. Error MCMC % p-value    
## edges                          -8.2209929  0.0210040      0  <1e-04 ***
## nodecov.timesCited             -0.0135754  0.0001080      0  <1e-04 ***
## nodecov.degree                  0.0147995  0.0000676      0  <1e-04 ***
## nodecov.numPub                  0.1284750  0.0008353      0  <1e-04 ***
## nodematch.communitySBM          5.6831643  0.0180715      0  <1e-04 ***
## nodematch.collabType            0.6146095  0.0128654      0  <1e-04 ***
## nodefactor.collabType.NATIONAL -0.3279974  0.0167088      0  <1e-04 ***
## nodefactor.collabType.REGIONAL  0.5783559  0.0102558      0  <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2224636  on 1604736  degrees of freedom
##  Residual Deviance:  217010  on 1604728  degrees of freedom
##  
## AIC: 217026    BIC: 217125    (Smaller is better.)
```




```r
save(model2, file = 'modData2/smodel2.rda')

#######################################################################################
#######################################################################################
#######################################################################################
```




```r
model3.f <- formula(mnet3.simple ~ edges + triangle + twopath)
model3.f
```

```
## mnet3.simple ~ edges + triangle + twopath
```

```r
summary.statistics(model3.f)
```

```
##    edges triangle  twopath 
##    95707  9074443 28224895
```




```r
t0 <- Sys.time()
model3 <- ergm(model3.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
```

```
## Starting maximum likelihood estimation via MCMLE:
## Iteration 1 of at most 1000: 
## The log-likelihood improved by 2.139 
## Iteration 2 of at most 1000: 
## The log-likelihood improved by 3.069 
## Iteration 3 of at most 1000: 
## The log-likelihood improved by 2.423 
## Iteration 4 of at most 1000: 
## The log-likelihood improved by 1.642 
## Iteration 5 of at most 1000: 
## The log-likelihood improved by 1.848 
## Iteration 6 of at most 1000: 
## The log-likelihood improved by 2.767 
## Iteration 7 of at most 1000: 
## The log-likelihood improved by 2.1 
## Iteration 8 of at most 1000: 
## The log-likelihood improved by 2.329 
## Iteration 9 of at most 1000: 
## The log-likelihood improved by 2.469 
## Iteration 10 of at most 1000: 
## The log-likelihood improved by 2.954 
## Iteration 11 of at most 1000: 
## The log-likelihood improved by 2.249 
## Iteration 12 of at most 1000: 
## The log-likelihood improved by 1.838 
## Iteration 13 of at most 1000: 
## The log-likelihood improved by 1.807 
## Iteration 14 of at most 1000: 
## The log-likelihood improved by 2.713 
## Iteration 15 of at most 1000: 
## The log-likelihood improved by 3.215 
## Iteration 16 of at most 1000: 
## The log-likelihood improved by 2.264 
## Iteration 17 of at most 1000: 
## The log-likelihood improved by 1.785 
## Iteration 18 of at most 1000: 
## The log-likelihood improved by 2.012 
## Iteration 19 of at most 1000: 
## The log-likelihood improved by 1.844 
## Iteration 20 of at most 1000: 
## The log-likelihood improved by 1.444 
## Iteration 21 of at most 1000: 
## The log-likelihood improved by 1.835 
## Iteration 22 of at most 1000: 
## The log-likelihood improved by 0.3141 
## Step length converged once. Increasing MCMC sample size.
## Iteration 23 of at most 1000: 
## The log-likelihood improved by 1.96 
## Step length converged twice. Stopping.
## Evaluating log-likelihood at the estimate. Using 20 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 .
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
Sys.time() - t0
```

```
## Time difference of 1.998849 hours
```




```r
anova.ergm(model3)
```

```
## Analysis of Variance Table
## 
## Model 1: mnet3.simple ~ edges + triangle + twopath
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                   1604736          0                 
## Model 1:  3  -658762   1604733     658762    < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(model3)
```

```
## 
## ==========================
## Summary of model fit
## ==========================
## 
## Formula:   mnet3.simple ~ edges + triangle + twopath
## 
## Iterations:  23 out of 1000 
## 
## Monte Carlo MLE Results:
##            Estimate Std. Error MCMC % p-value    
## edges    -4.1386030  0.0126549      1  <1e-04 ***
## triangle  0.0937380  0.0003766      1  <1e-04 ***
## twopath  -0.0069399  0.0001070      1  <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2224636  on 1604736  degrees of freedom
##  Residual Deviance:  658762  on 1604733  degrees of freedom
##  
## AIC: 658768    BIC: 658805    (Smaller is better.)
```




```r
mcmc.model3 <- mcmc.diagnostics(model3)
```

```
## Sample statistics summary:
## 
## Iterations = 20000:61500
## Thinning interval = 500 
## Number of chains = 48 
## Sample size per chain = 84 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##             Mean      SD Naive SE Time-series SE
## edges      180.6   111.4    1.755          6.521
## triangle  4949.5  3210.4   50.559        135.293
## twopath  23261.7 14975.9  235.849        532.779
## 
## 2. Quantiles for each variable:
## 
##              2.5%   25%   50%   75% 97.5%
## edges      -35.22   106   178   250   406
## triangle -1089.58  2386  5365  7286 10921
## twopath  -3866.15 11361 24372 34533 49837
## 
## 
## Sample statistics cross-correlations:
##              edges  triangle   twopath
## edges    1.0000000 0.2781246 0.6871841
## triangle 0.2781246 1.0000000 0.5324242
## twopath  0.6871841 0.5324242 1.0000000
## 
## Sample statistics auto-correlation:
## Chain 1 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9745953 0.9524093 0.9623175
## Lag 1000 0.9467054 0.9062735 0.9283242
## Lag 1500 0.9193117 0.8667812 0.9026776
## Lag 2000 0.8916326 0.8250904 0.8733591
## Lag 2500 0.8626808 0.7872278 0.8492227
## Chain 2 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9713419 0.9600303 0.9219179
## Lag 1000 0.9411031 0.9212978 0.8446570
## Lag 1500 0.9139982 0.8812358 0.7738417
## Lag 2000 0.8867630 0.8405053 0.7076666
## Lag 2500 0.8592553 0.8008749 0.6413243
## Chain 3 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9391771 0.9632380 0.9485334
## Lag 1000 0.8779955 0.9252629 0.8993309
## Lag 1500 0.8220507 0.8883211 0.8599444
## Lag 2000 0.7805685 0.8540119 0.8318046
## Lag 2500 0.7422811 0.8221278 0.8006467
## Chain 4 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9148918 0.9575516 0.8878855
## Lag 1000 0.8491321 0.9161611 0.8129421
## Lag 1500 0.7935754 0.8757411 0.7455165
## Lag 2000 0.7275686 0.8338180 0.6600754
## Lag 2500 0.6730387 0.7887688 0.5816836
## Chain 5 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9465954 0.9618027 0.9413567
## Lag 1000 0.8876280 0.9211129 0.8799198
## Lag 1500 0.8235435 0.8833583 0.8193447
## Lag 2000 0.7638815 0.8469208 0.7698009
## Lag 2500 0.6980688 0.8109526 0.7095089
## Chain 6 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9417677 0.9635470 0.9143483
## Lag 1000 0.8774879 0.9251452 0.8253400
## Lag 1500 0.8206521 0.8854890 0.7517602
## Lag 2000 0.7648756 0.8439840 0.6868254
## Lag 2500 0.7119323 0.8011477 0.6164680
## Chain 7 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9581502 0.9634797 0.9409631
## Lag 1000 0.9213117 0.9271814 0.8828344
## Lag 1500 0.8815010 0.8918177 0.8207982
## Lag 2000 0.8458021 0.8586286 0.7671280
## Lag 2500 0.8104058 0.8235528 0.7312018
## Chain 8 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9663385 0.9617605 0.9468522
## Lag 1000 0.9292870 0.9212559 0.8892337
## Lag 1500 0.8945022 0.8818406 0.8369147
## Lag 2000 0.8613024 0.8409887 0.7881496
## Lag 2500 0.8242777 0.7987253 0.7388603
## Chain 9 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9651377 0.9630779 0.9677550
## Lag 1000 0.9294298 0.9247260 0.9331703
## Lag 1500 0.8925475 0.8852599 0.8971003
## Lag 2000 0.8555484 0.8461688 0.8582799
## Lag 2500 0.8151712 0.8133827 0.8200570
## Chain 10 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9737649 0.9691319 0.9784129
## Lag 1000 0.9453137 0.9370340 0.9514213
## Lag 1500 0.9169910 0.9077804 0.9253416
## Lag 2000 0.8867246 0.8777791 0.8981219
## Lag 2500 0.8509512 0.8463526 0.8709248
## Chain 11 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9243773 0.9623270 0.9271797
## Lag 1000 0.8718125 0.9255064 0.8793237
## Lag 1500 0.8062865 0.8899341 0.8187823
## Lag 2000 0.7354008 0.8525093 0.7540340
## Lag 2500 0.6573787 0.8132380 0.6847085
## Chain 12 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9624188 0.9629014 0.9638979
## Lag 1000 0.9212997 0.9243045 0.9273333
## Lag 1500 0.8790380 0.8853293 0.8884434
## Lag 2000 0.8346394 0.8438667 0.8470101
## Lag 2500 0.7901219 0.8003803 0.8065683
## Chain 13 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9546235 0.9657359 0.9726758
## Lag 1000 0.9107653 0.9321094 0.9453352
## Lag 1500 0.8683789 0.8972886 0.9136983
## Lag 2000 0.8225512 0.8612449 0.8842983
## Lag 2500 0.7728892 0.8240003 0.8531888
## Chain 14 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9488235 0.9582449 0.9629898
## Lag 1000 0.9042228 0.9161171 0.9338274
## Lag 1500 0.8559923 0.8741801 0.9053461
## Lag 2000 0.7960028 0.8356675 0.8757920
## Lag 2500 0.7400623 0.8010291 0.8464596
## Chain 15 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9429761 0.9631053 0.9404758
## Lag 1000 0.8885115 0.9294206 0.8859620
## Lag 1500 0.8352470 0.8953314 0.8359632
## Lag 2000 0.7886823 0.8611380 0.7816270
## Lag 2500 0.7404442 0.8258391 0.7162995
## Chain 16 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9620219 0.9403111 0.9611653
## Lag 1000 0.9200346 0.8827087 0.9185888
## Lag 1500 0.8732602 0.8312374 0.8809092
## Lag 2000 0.8268312 0.7811342 0.8407223
## Lag 2500 0.7765585 0.7288328 0.8033920
## Chain 17 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9263823 0.9689157 0.9564853
## Lag 1000 0.8564290 0.9363920 0.9047432
## Lag 1500 0.7879531 0.9044431 0.8545242
## Lag 2000 0.7261092 0.8724762 0.8047775
## Lag 2500 0.6641275 0.8399641 0.7579750
## Chain 18 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9702926 0.9690555 0.9736483
## Lag 1000 0.9330819 0.9387908 0.9417275
## Lag 1500 0.8953497 0.9070743 0.9092329
## Lag 2000 0.8570543 0.8731783 0.8740170
## Lag 2500 0.8180000 0.8417351 0.8393417
## Chain 19 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9187889 0.9558517 0.8891618
## Lag 1000 0.8412523 0.9094538 0.7941769
## Lag 1500 0.7819805 0.8686410 0.7264746
## Lag 2000 0.7226549 0.8292755 0.6592449
## Lag 2500 0.6663381 0.7877931 0.6049731
## Chain 20 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9704963 0.9642208 0.9755033
## Lag 1000 0.9394971 0.9277623 0.9516180
## Lag 1500 0.9060056 0.8940933 0.9268726
## Lag 2000 0.8742964 0.8616443 0.9024789
## Lag 2500 0.8407602 0.8320611 0.8783509
## Chain 21 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9561292 0.9612402 0.8822989
## Lag 1000 0.9162976 0.9238217 0.7890188
## Lag 1500 0.8773128 0.8854002 0.7158686
## Lag 2000 0.8296060 0.8474215 0.6524501
## Lag 2500 0.7894696 0.8152339 0.6069249
## Chain 22 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9366925 0.9681676 0.9424379
## Lag 1000 0.8724690 0.9385828 0.8934530
## Lag 1500 0.8180855 0.9101906 0.8491311
## Lag 2000 0.7647826 0.8805451 0.8123460
## Lag 2500 0.7041055 0.8490634 0.7709338
## Chain 23 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9539456 0.9474150 0.9445483
## Lag 1000 0.9149246 0.8941190 0.8917653
## Lag 1500 0.8826069 0.8470479 0.8441680
## Lag 2000 0.8438038 0.8001650 0.7996476
## Lag 2500 0.8065993 0.7562991 0.7576647
## Chain 24 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9224284 0.9620312 0.8850425
## Lag 1000 0.8583624 0.9263721 0.7847986
## Lag 1500 0.7914089 0.8911466 0.6721670
## Lag 2000 0.7184092 0.8555780 0.5568890
## Lag 2500 0.6475647 0.8187971 0.4438692
## Chain 25 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9490175 0.9681302 0.8821128
## Lag 1000 0.8951662 0.9364536 0.7531344
## Lag 1500 0.8537520 0.9064423 0.6595116
## Lag 2000 0.7999487 0.8755354 0.5778852
## Lag 2500 0.7489662 0.8430967 0.4719756
## Chain 26 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9484169 0.9591258 0.9466403
## Lag 1000 0.8935266 0.9164421 0.8964768
## Lag 1500 0.8457123 0.8714566 0.8577427
## Lag 2000 0.7969009 0.8287526 0.8241141
## Lag 2500 0.7607509 0.7929787 0.7936800
## Chain 27 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.8965449 0.9584912 0.9058383
## Lag 1000 0.7855672 0.9222252 0.8353042
## Lag 1500 0.7027027 0.8840063 0.7585613
## Lag 2000 0.6429195 0.8465731 0.6820658
## Lag 2500 0.6004471 0.8103405 0.6263977
## Chain 28 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.8639984 0.9640426 0.8971741
## Lag 1000 0.7617330 0.9292967 0.8147324
## Lag 1500 0.6809528 0.8928554 0.7445785
## Lag 2000 0.5976159 0.8545920 0.6667919
## Lag 2500 0.5453623 0.8180067 0.6212812
## Chain 29 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9275702 0.9663868 0.8597665
## Lag 1000 0.8761071 0.9352691 0.7746208
## Lag 1500 0.8330567 0.9028850 0.7133239
## Lag 2000 0.7895400 0.8687494 0.6652956
## Lag 2500 0.7516907 0.8359027 0.6198364
## Chain 30 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9626691 0.9668548 0.9493384
## Lag 1000 0.9295073 0.9323970 0.9000183
## Lag 1500 0.8975233 0.8992261 0.8496589
## Lag 2000 0.8634285 0.8649133 0.7966906
## Lag 2500 0.8247992 0.8291464 0.7378357
## Chain 31 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9712450 0.9545439 0.9410209
## Lag 1000 0.9424043 0.9146316 0.8911733
## Lag 1500 0.9143784 0.8783954 0.8524123
## Lag 2000 0.8842563 0.8399223 0.8165247
## Lag 2500 0.8543754 0.8020549 0.7713602
## Chain 32 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9510284 0.9660277 0.9272571
## Lag 1000 0.9062092 0.9326767 0.8697695
## Lag 1500 0.8686013 0.8994437 0.8160367
## Lag 2000 0.8298900 0.8636762 0.7529098
## Lag 2500 0.7920798 0.8272698 0.6904819
## Chain 33 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9608911 0.9710728 0.9537195
## Lag 1000 0.9202502 0.9440442 0.9057777
## Lag 1500 0.8794695 0.9159452 0.8584679
## Lag 2000 0.8438102 0.8865128 0.8225174
## Lag 2500 0.8098909 0.8563503 0.7920285
## Chain 34 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9680697 0.9690118 0.9581942
## Lag 1000 0.9367417 0.9370692 0.9180833
## Lag 1500 0.9015323 0.9050639 0.8738879
## Lag 2000 0.8708278 0.8714707 0.8339495
## Lag 2500 0.8410226 0.8371421 0.7988695
## Chain 35 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9728774 0.9633023 0.9573679
## Lag 1000 0.9407786 0.9251110 0.9152285
## Lag 1500 0.9018779 0.8857896 0.8590587
## Lag 2000 0.8632414 0.8445198 0.8080770
## Lag 2500 0.8221912 0.8033823 0.7597467
## Chain 36 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9202628 0.9660034 0.9619432
## Lag 1000 0.8389319 0.9313391 0.9136056
## Lag 1500 0.7580116 0.8996388 0.8625882
## Lag 2000 0.6655490 0.8656512 0.8067440
## Lag 2500 0.6043237 0.8308977 0.7540742
## Chain 37 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9655003 0.9629891 0.9586966
## Lag 1000 0.9240742 0.9261503 0.9165965
## Lag 1500 0.8822197 0.8878734 0.8775769
## Lag 2000 0.8403948 0.8498768 0.8400808
## Lag 2500 0.7988553 0.8152334 0.8094456
## Chain 38 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9455060 0.9636000 0.9217097
## Lag 1000 0.8962054 0.9270297 0.8465772
## Lag 1500 0.8396204 0.8907634 0.7545780
## Lag 2000 0.7905929 0.8556192 0.6763423
## Lag 2500 0.7487186 0.8194511 0.6202823
## Chain 39 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9438285 0.9610059 0.9091363
## Lag 1000 0.8977128 0.9234458 0.8425436
## Lag 1500 0.8450707 0.8839957 0.7714709
## Lag 2000 0.8046898 0.8431299 0.6985260
## Lag 2500 0.7565578 0.8019336 0.6312218
## Chain 40 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9668072 0.9534040 0.9334788
## Lag 1000 0.9355638 0.9088641 0.8741945
## Lag 1500 0.9063598 0.8693312 0.8253073
## Lag 2000 0.8741759 0.8272153 0.7556308
## Lag 2500 0.8408997 0.7889715 0.6764792
## Chain 41 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9121022 0.9575029 0.8571139
## Lag 1000 0.8165224 0.9136121 0.7044475
## Lag 1500 0.7196362 0.8712308 0.5871791
## Lag 2000 0.6383842 0.8336798 0.5089009
## Lag 2500 0.5836930 0.7975761 0.4800471
## Chain 42 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.7661200 0.9602736 0.9183231
## Lag 1000 0.5861403 0.9228751 0.8453760
## Lag 1500 0.4354245 0.8843249 0.7740007
## Lag 2000 0.3336248 0.8448061 0.7207077
## Lag 2500 0.2800213 0.8076420 0.6634333
## Chain 43 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9699669 0.9639922 0.9304108
## Lag 1000 0.9367055 0.9285499 0.8644598
## Lag 1500 0.9025531 0.8919597 0.8099055
## Lag 2000 0.8684131 0.8533722 0.7560748
## Lag 2500 0.8355929 0.8153867 0.7093863
## Chain 44 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.7858555 0.9568752 0.9014923
## Lag 1000 0.5554843 0.9106870 0.8011920
## Lag 1500 0.4412617 0.8650223 0.7414452
## Lag 2000 0.3485452 0.8228251 0.6868688
## Lag 2500 0.2499423 0.7787448 0.6326713
## Chain 45 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9524064 0.9632107 0.9295305
## Lag 1000 0.8993024 0.9270283 0.8491125
## Lag 1500 0.8512024 0.8895146 0.7824681
## Lag 2000 0.7992148 0.8531131 0.7166914
## Lag 2500 0.7544149 0.8161653 0.6569025
## Chain 46 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9281704 0.9618144 0.8713613
## Lag 1000 0.8550198 0.9248943 0.7776524
## Lag 1500 0.7867607 0.8921744 0.7225103
## Lag 2000 0.7216190 0.8585729 0.6626176
## Lag 2500 0.6613599 0.8237366 0.5848811
## Chain 47 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9554988 0.9567645 0.9562795
## Lag 1000 0.9065380 0.9113914 0.9185847
## Lag 1500 0.8601785 0.8602647 0.8748252
## Lag 2000 0.8216951 0.8126195 0.8335072
## Lag 2500 0.7888845 0.7660244 0.7809472
## Chain 48 
##              edges  triangle   twopath
## Lag 0    1.0000000 1.0000000 1.0000000
## Lag 500  0.9581234 0.9597004 0.9537961
## Lag 1000 0.9214989 0.9198634 0.9087878
## Lag 1500 0.8827211 0.8837994 0.8687999
## Lag 2000 0.8462808 0.8468103 0.8299526
## Lag 2500 0.8169768 0.8100600 0.7919665
## 
## Sample statistics burn-in diagnostic (Geweke):
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-10.5762492174221) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.19862198539223) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-46.912979533284) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.95472348869831) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.40105870299484) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-0.591930719077546) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (0.120754053891087) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-4.13416920392355) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-5.84183283087564) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-13.745152169446) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-11.8120192013421) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-5.73876387654383) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.5912609676227) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.22645166025596) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-20.6525377093871) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.19421676631191) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-3.10863748704727) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.12742243936345) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-7.42062578991247) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-0.931200042753602) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-3.51675554901897) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-5.91076784575181) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.91904210645507) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-6.14542494734144) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-0.453808613237977) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-12.325107316376) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-5.1378323552568) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-12.7350344701157) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.64502843282461) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.84351604607337) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.10017540830693) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-13.8574610540322) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-6.1840007444372) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-159.917169032231) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-3.0955388087167) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.16500846228769) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-8.61466304262925) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.75369849882616) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-3.90942232754214) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-6.0929711450913) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.75352200493721) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.18107035705276) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-2.4591481158038) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-5.76368079404408) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (0.658162896159733) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-1.11707645065572) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-7.88911410195877) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Warning in approx.hotelling.diff.test(x1, x2, var.equal = TRUE): Effective
## degrees of freedom (-10.0281804456948) must exceed the number of varying
## parameters (3). P-value will not be computed.
```

```
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -3.954   -4.244   -2.033 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 7.675341e-05 2.197204e-05 4.205752e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.050   -3.482   -2.804 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.446621e-09 4.978478e-04 5.049575e-03 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -10.279   -6.273  -14.973 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 8.765474e-25 3.533912e-10 1.109946e-50 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -3.005   -3.338   -1.701 
## 
## Individual P-values (lower = worse):
##       edges    triangle     twopath 
## 0.002655731 0.000844742 0.088989116 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.773   -4.947  -32.908 
## 
## Individual P-values (lower = worse):
##         edges      triangle       twopath 
##  1.257346e-11  7.554913e-07 1.672037e-237 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -8.148   -6.029   -6.531 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.709134e-16 1.650195e-09 6.542409e-11 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.617   -8.473  -10.031 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.673642e-11 2.384156e-17 1.109646e-23 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.442   -2.914   -2.127 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.176207e-10 3.566325e-03 3.337924e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -7.033   -4.819   -6.346 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.016850e-12 1.444495e-06 2.205230e-10 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.493   -8.219   -7.297 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.266626e-02 2.054653e-16 2.953381e-13 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.728   -4.675   -8.414 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.268058e-06 2.936740e-06 3.949775e-17 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -3.704   -4.492   -5.113 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.124790e-04 7.069640e-06 3.175773e-07 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.728   -5.761   -4.962 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.269629e-06 8.369959e-09 6.982393e-07 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.031   -5.646  -23.157 
## 
## Individual P-values (lower = worse):
##         edges      triangle       twopath 
##  5.554549e-05  1.647009e-08 1.231942e-118 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.384   -5.700   -5.802 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.164709e-05 1.200632e-08 6.541672e-09 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.129   -2.925  -10.161 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.644910e-05 3.440824e-03 2.973030e-24 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -7.301   -7.413  -25.653 
## 
## Individual P-values (lower = worse):
##         edges      triangle       twopath 
##  2.861352e-13  1.232006e-13 3.904285e-145 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -17.747   -6.488  -18.417 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.805102e-70 8.723578e-11 9.685058e-76 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -5.234   -6.839   -4.666 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.660261e-07 7.976981e-12 3.070455e-06 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -11.485   -7.511  -10.424 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.563203e-30 5.866687e-14 1.935588e-25 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.945   -5.813   -3.925 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.774644e-12 6.119282e-09 8.672320e-05 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##    -4.16   -16.94   -14.52 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.179068e-05 2.451107e-64 8.752284e-48 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -12.882   -7.866  -10.347 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 5.675382e-38 3.652142e-15 4.316421e-25 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -3.125   -4.486   -3.115 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.779901e-03 7.253460e-06 1.840640e-03 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.034   -4.379   -1.677 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.598182e-09 1.194604e-05 9.352540e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.174   -6.087   -3.348 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.971184e-02 1.153713e-09 8.147347e-04 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.769   -4.222   -3.534 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.848703e-06 2.425849e-05 4.096931e-04 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.555   -6.472    1.304 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.060743e-02 9.690812e-11 1.922301e-01 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -5.088   -4.950   -3.209 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.626811e-07 7.406282e-07 1.330843e-03 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.266   -5.413   -3.414 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.994378e-05 6.203994e-08 6.409490e-04 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -8.108   -5.513   -4.595 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 5.158000e-16 3.536484e-08 4.330377e-06 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.213   -4.311   -2.107 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 5.191313e-10 1.626081e-05 3.511823e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.301   -6.166   -3.359 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.696571e-05 6.991899e-10 7.821590e-04 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.923   -6.232   -4.743 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 8.502417e-07 4.609598e-10 2.109806e-06 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -5.932   -5.045   -4.748 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.989730e-09 4.525841e-07 2.049462e-06 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##    2.049   -7.151    1.911 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 4.043045e-02 8.636484e-13 5.602719e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -4.882   -5.041   -4.026 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.052032e-06 4.637386e-07 5.676644e-05 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -3.522   -5.139   -7.088 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 4.288787e-04 2.766823e-07 1.356384e-12 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -10.82    -5.27   -15.05 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.668164e-27 1.361063e-07 3.372094e-51 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -7.328   -4.749   -1.860 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.334275e-13 2.041305e-06 6.284779e-02 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -0.9751  -4.2348  -5.3967 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 3.295261e-01 2.287695e-05 6.787240e-08 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.459   -8.100   -2.918 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.393752e-02 5.481998e-16 3.519589e-03 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -16.199   -5.495  -17.327 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 5.114549e-59 3.917559e-08 2.946161e-67 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.550   -5.254   -6.419 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.077931e-02 1.490627e-07 1.372959e-10 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##  -10.178   -5.676   -2.582 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 2.475547e-24 1.378253e-08 9.808854e-03 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.606   -4.608   -7.006 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 9.167041e-03 4.065361e-06 2.454557e-12 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -2.827   -2.845   -1.518 
## 
## Individual P-values (lower = worse):
##       edges    triangle     twopath 
## 0.004693983 0.004439143 0.129108490 
## Joint P-value (lower = worse):  0 .
## Chain 48 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##    edges triangle  twopath 
##   -6.396   -5.370   -8.307 
## 
## Individual P-values (lower = worse):
##        edges     triangle      twopath 
## 1.592039e-10 7.889141e-08 9.815597e-17 
## Joint P-value (lower = worse):  0 .
```

```
## Loading required namespace: latticeExtra
```

![plot of chunk mod3-mcmc_diag](figure/mod3-mcmc_diag-1.png)

```
## 
## MCMC diagnostics shown here are from the last round of simulation, prior to computation of final parameter estimates. Because the final estimates are refinements of those used for this simulation run, these diagnostics may understate model performance. To directly assess the performance of the final model on in-model statistics, please use the GOF command: gof(ergmFitObject, GOF=~model).
```

```r
mcmc.model3
```

```
## $degeneracy.value
## NULL
## 
## $degeneracy.type
## NULL
```

```r
save(mcmc.model3, file = "modData/mcmc-smodel3.rda")
```




```r
save(model3, file = 'modData2/smodel3.rda')

#######################################################################################
#######################################################################################
#######################################################################################
```




```r
model4.f <- formula(mnet3.simple ~ edges + triangle + twopath + nodecov("timesCited") + nodecov('degree') 
                        + nodecov("numPub") + nodematch("communitySBM") + nodematch("collabType") 
                        + nodefactor("collabType"))
model4.f
```

```
## mnet3.simple ~ edges + triangle + twopath + nodecov("timesCited") + 
##     nodecov("degree") + nodecov("numPub") + nodematch("communitySBM") + 
##     nodematch("collabType") + nodefactor("collabType")
```

```r
summary.statistics(model4.f)
```

```
##                          edges                       triangle 
##                          95707                        9074443 
##                        twopath             nodecov.timesCited 
##                       28224895                       33146116 
##                 nodecov.degree                 nodecov.numPub 
##                       56641204                         495605 
##         nodematch.communitySBM           nodematch.collabType 
##                          75410                          58334 
## nodefactor.collabType.NATIONAL nodefactor.collabType.REGIONAL 
##                           6124                          30809
```




```r
t0 <- Sys.time()
model4 <- ergm(model4.f, iterations=1000, 
                         control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=1000,
                                              MCMC.interval=500, MCMLE.maxit = 1000, MCMC.burnin = 20000, seed=25))
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: algorithm did not convergeglm.fit: fitted probabilities
## numerically 0 or 1 occurred
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## Starting maximum likelihood estimation via MCMLE:
## Iteration 1 of at most 1000:
```

```
## Error in rbind(if (!is.null(x2)) t(gamma * t(x2crs) + (1 - margin * gamma) * : object 'm2crs' not found
```

```r
Sys.time() - t0
```

```
## Time difference of 4.646678 mins
```




```r
anova.ergm(model4)
```

```
## Error in anova.ergm(model4): object 'model4' not found
```

```r
summary(model4)
```

```
## Error in summary(model4): object 'model4' not found
```




```r
mcmc.model4 <- mcmc.diagnostics(model4)
```

```
## Error in mcmc.diagnostics(model4): object 'model4' not found
```

```r
mcmc.model4
```

```
## Error in eval(expr, envir, enclos): object 'mcmc.model4' not found
```

```r
save(mcmc.model4, file = "modData/mcmc-smodel4.rda")
```

```
## Error in save(mcmc.model4, file = "modData/mcmc-smodel4.rda"): object 'mcmc.model4' not found
```




```r
save(model4, file = 'modData2/smodel4.rda')
```

```
## Error in save(model4, file = "modData2/smodel4.rda"): object 'model4' not found
```




```r
t0 <- Sys.time()
diag.model2 <- gof(model2, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
```

```
## Warning: closing unused connection 53 (<-localhost:11785)
```

```
## Warning: closing unused connection 52 (<-localhost:11785)
```

```
## Warning: closing unused connection 51 (<-localhost:11785)
```

```
## Warning: closing unused connection 50 (<-localhost:11785)
```

```
## Warning: closing unused connection 49 (<-localhost:11785)
```

```
## Warning: closing unused connection 48 (<-localhost:11785)
```

```
## Warning: closing unused connection 47 (<-localhost:11785)
```

```
## Warning: closing unused connection 46 (<-localhost:11785)
```

```
## Warning: closing unused connection 45 (<-localhost:11785)
```

```
## Warning: closing unused connection 44 (<-localhost:11785)
```

```
## Warning: closing unused connection 43 (<-localhost:11785)
```

```
## Warning: closing unused connection 42 (<-localhost:11785)
```

```
## Warning: closing unused connection 41 (<-localhost:11785)
```

```
## Warning: closing unused connection 40 (<-localhost:11785)
```

```
## Warning: closing unused connection 39 (<-localhost:11785)
```

```
## Warning: closing unused connection 38 (<-localhost:11785)
```

```
## Warning: closing unused connection 37 (<-localhost:11785)
```

```
## Warning: closing unused connection 36 (<-localhost:11785)
```

```
## Warning: closing unused connection 35 (<-localhost:11785)
```

```
## Warning: closing unused connection 34 (<-localhost:11785)
```

```
## Warning: closing unused connection 33 (<-localhost:11785)
```

```
## Warning: closing unused connection 32 (<-localhost:11785)
```

```
## Warning: closing unused connection 31 (<-localhost:11785)
```

```
## Warning: closing unused connection 30 (<-localhost:11785)
```

```
## Warning: closing unused connection 29 (<-localhost:11785)
```

```
## Warning: closing unused connection 28 (<-localhost:11785)
```

```
## Warning: closing unused connection 27 (<-localhost:11785)
```

```
## Warning: closing unused connection 26 (<-localhost:11785)
```

```
## Warning: closing unused connection 25 (<-localhost:11785)
```

```
## Warning: closing unused connection 24 (<-localhost:11785)
```

```
## Warning: closing unused connection 23 (<-localhost:11785)
```

```
## Warning: closing unused connection 22 (<-localhost:11785)
```

```
## Warning: closing unused connection 21 (<-localhost:11785)
```

```
## Warning: closing unused connection 20 (<-localhost:11785)
```

```
## Warning: closing unused connection 19 (<-localhost:11785)
```

```
## Warning: closing unused connection 18 (<-localhost:11785)
```

```
## Warning: closing unused connection 17 (<-localhost:11785)
```

```
## Warning: closing unused connection 16 (<-localhost:11785)
```

```
## Warning: closing unused connection 15 (<-localhost:11785)
```

```
## Warning: closing unused connection 14 (<-localhost:11785)
```

```
## Warning: closing unused connection 13 (<-localhost:11785)
```

```
## Warning: closing unused connection 12 (<-localhost:11785)
```

```
## Warning: closing unused connection 11 (<-localhost:11785)
```

```
## Warning: closing unused connection 10 (<-localhost:11785)
```

```
## Warning: closing unused connection 9 (<-localhost:11785)
```

```
## Warning: closing unused connection 8 (<-localhost:11785)
```

```
## Warning: closing unused connection 7 (<-localhost:11785)
```

```
## Warning: closing unused connection 6 (<-localhost:11785)
```

```r
Sys.time() - t0
```

```
## Time difference of 14.82796 hours
```

```r
save(diag.model2, file = 'modData2/diag.smodel2.rda')
jpeg(filename = './modData2/gof-smodel2.jpg')
plot(diag.model2)
dev.off()
```

```
## png 
##   2
```




```r
t0 <- Sys.time()
diag.model4 <- gof(model4, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
```

```
## Error in gof(model4, control = control.gof.ergm(parallel = cores, parallel.type = "PSOCK")): object 'model4' not found
```

```r
Sys.time() - t0
```

```
## Time difference of 0.002240181 secs
```

```r
save(diag.model4, file = 'modData2/diag.smodel4.rda')
```

```
## Error in save(diag.model4, file = "modData2/diag.smodel4.rda"): object 'diag.model4' not found
```

```r
jpeg(filename = './modData2/gof-smodel4.jpg')
plot(diag.model4)
```

```
## Error in plot(diag.model4): object 'diag.model4' not found
```

```r
dev.off()
```

```
## png 
##   2
```

```r
#######################################################################################
#################################### END ##############################################
#######################################################################################
```

