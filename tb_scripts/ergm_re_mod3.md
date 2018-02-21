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

mnet3 <- read_graph('./graphs/CAnet_graph.graphml', format = 'graphml')
```

Since this processing is going to take a lot of time, we decide to parallelize all the computations:


```r
cores=detectCores() # get number of cores
```

```r
A <- get.adjacency(mnet3)
v.attrs <- get.data.frame(mnet3, what="vertices")
v.attrs$name <- V(mnet2)$name
v.attrs$timesCited <- V(mnet2)$timesCited
v.attrs$numPub <- V(mnet2)$numPub
v.attrs$community <- V(mnet2)$community
v.attrs$degree <- V(mnet2)$degree
v.attrs$affiliation <- V(mnet2)$place
v.attrs$city <- V(mnet2)$affil
v.attrs$country <- V(mnet2)$country

mnet3.simple <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(mnet3.simple, "timesCited", V(mnet2)$timesCited)
network::set.vertex.attribute(mnet3.simple, "numPub", V(mnet2)$numPub)
network::set.vertex.attribute(mnet3.simple, "community", V(mnet2)$community)
network::set.vertex.attribute(mnet3.simple, "degree", V(mnet2)$degree)
```

```r
mnet3.ergm <- formula(mnet3.simple ~ edges + nodecov("timesCited") + 
                        nodecov("numPub") + match("community"))
summary.statistics(mnet3.ergm)
```

```
##               edges  nodecov.timesCited      nodecov.numPub 
##               95707            33146116              495605 
## nodematch.community 
##               92011
```

```r
t0 <- Sys.time()
mnet3.ergm.fit <- ergm(mnet3.ergm, iterations=100,
                       control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=100,
                                            MCMC.interval=500, MCMC.burnin = 20000, seed=25, MCMLE.maxit = 100))
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init, : GLM model may be separable; restarting glm with zeros.

## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init, : GLM model may be separable; restarting glm with zeros.
```

```
## Evaluating log-likelihood at the estimate.
```

```r
Sys.time() - t0
```

```
## Time difference of 13.44294 secs
```

```r
save(mnet3.ergm.fit, file = './Rdata/mnet3.ergm.fit.rda')
```

```r
anova.ergm(mnet3.ergm.fit)
```

```
## Analysis of Variance Table
## 
## Model 1: mnet3.simple ~ edges + nodecov("timesCited") + nodecov("numPub") + 
##     match("community")
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                   1604736          0                 
## Model 1:  4  -309019   1604732     309019    < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(mnet3.ergm.fit)
```

```
## 
## ==========================
## Summary of model fit
## ==========================
## 
## Formula:   mnet3.simple ~ edges + nodecov("timesCited") + nodecov("numPub") + 
##     match("community")
## 
## Iterations:  9 out of 100 
## 
## Monte Carlo MLE Results:
##                       Estimate Std. Error MCMC % p-value    
## edges               -7.078e+00  1.917e-02      0  <1e-04 ***
## nodecov.timesCited   8.103e-03  2.497e-05      0  <1e-04 ***
## nodecov.numPub      -1.364e-01  6.851e-04      0  <1e-04 ***
## nodematch.community  5.258e+00  1.862e-02      0  <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2224636  on 1604736  degrees of freedom
##  Residual Deviance:  309019  on 1604732  degrees of freedom
##  
## AIC: 309027    BIC: 309076    (Smaller is better.)
```

```r
t0 <- Sys.time()
diag <- gof(mnet3.ergm.fit, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
save(diag, file = './Rdata/diag.mnet3.ergm.fit.rda')
Sys.time() - t0
```

```
## Time difference of 13.76587 hours
```

```r
mcmc.diagnostics(mnet3.ergm.fit)
```

```
## Error in mcmc.diagnostics.ergm(mnet3.ergm.fit): MCMC was not run or MCMC sample was not stored.
```

```r
par(mfrow=c(2, 2))
plot(diag)
```

![plot of chunk gof-viz](figure/gof-viz-1.png)

