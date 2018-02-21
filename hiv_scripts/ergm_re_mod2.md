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
my.ergm <- formula(mnet3.simple ~ edges + kstar(2) + kstar(3) + triangle)
summary.statistics(my.ergm)
```

```
##      edges     kstar2     kstar3   triangle 
##      95707   28224895 3207203947    9074443
```

```r
t0 <- Sys.time()
my.ergm.fit <- ergm(my.ergm, iterations=100,
                    control=control.ergm(parallel=cores, parallel.type="PSOCK", MCMC.samplesize=100,
                                         MCMC.interval=500, MCMLE.maxit = 100, MCMC.burnin = 20000, seed=25))
```

```
## Warning in ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init =
## init, : glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## Starting maximum likelihood estimation via MCMLE:
## Iteration 1 of at most 100: 
## The log-likelihood improved by 2.477 
## Iteration 2 of at most 100: 
## The log-likelihood improved by 1.848 
## Iteration 3 of at most 100: 
## The log-likelihood improved by 1.644 
## Iteration 4 of at most 100: 
## The log-likelihood improved by 2.14 
## Iteration 5 of at most 100: 
## The log-likelihood improved by 1.588 
## Iteration 6 of at most 100: 
## The log-likelihood improved by 1.854 
## Iteration 7 of at most 100: 
## The log-likelihood improved by 2.184 
## Iteration 8 of at most 100: 
## The log-likelihood improved by 1.296 
## Iteration 9 of at most 100: 
## The log-likelihood improved by 1.675 
## Iteration 10 of at most 100: 
## The log-likelihood improved by 1.786 
## Iteration 11 of at most 100: 
## The log-likelihood improved by 1.494 
## Iteration 12 of at most 100: 
## The log-likelihood improved by 1.232 
## Iteration 13 of at most 100: 
## The log-likelihood improved by 1.054 
## Iteration 14 of at most 100: 
## The log-likelihood improved by 1.596 
## Iteration 15 of at most 100: 
## The log-likelihood improved by 1.356 
## Iteration 16 of at most 100: 
## The log-likelihood improved by 1.397 
## Iteration 17 of at most 100: 
## The log-likelihood improved by 1.775 
## Iteration 18 of at most 100: 
## The log-likelihood improved by 1.019 
## Iteration 19 of at most 100: 
## The log-likelihood improved by 1.236 
## Iteration 20 of at most 100: 
## The log-likelihood improved by 1.212 
## Iteration 21 of at most 100: 
## The log-likelihood improved by 1.451 
## Iteration 22 of at most 100: 
## The log-likelihood improved by 1.766 
## Iteration 23 of at most 100: 
## The log-likelihood improved by 2.796 
## Iteration 24 of at most 100: 
## The log-likelihood improved by 2.064 
## Iteration 25 of at most 100: 
## The log-likelihood improved by 2.449 
## Iteration 26 of at most 100: 
## The log-likelihood improved by 1.709 
## Iteration 27 of at most 100: 
## The log-likelihood improved by 1.557 
## Iteration 28 of at most 100: 
## The log-likelihood improved by 1.341 
## Iteration 29 of at most 100: 
## The log-likelihood improved by 1.56 
## Iteration 30 of at most 100: 
## The log-likelihood improved by 2.551 
## Iteration 31 of at most 100: 
## The log-likelihood improved by 1.544 
## Iteration 32 of at most 100: 
## The log-likelihood improved by 1.771 
## Iteration 33 of at most 100: 
## The log-likelihood improved by 1.732 
## Iteration 34 of at most 100: 
## The log-likelihood improved by 1.186 
## Iteration 35 of at most 100: 
## The log-likelihood improved by 1.968 
## Iteration 36 of at most 100: 
## The log-likelihood improved by 1.887 
## Iteration 37 of at most 100: 
## The log-likelihood improved by 1.183 
## Iteration 38 of at most 100: 
## The log-likelihood improved by 1.378 
## Iteration 39 of at most 100: 
## The log-likelihood improved by 1.537 
## Iteration 40 of at most 100: 
## The log-likelihood improved by 1.608 
## Iteration 41 of at most 100: 
## The log-likelihood improved by 1.439 
## Iteration 42 of at most 100: 
## The log-likelihood improved by 1.35 
## Iteration 43 of at most 100: 
## The log-likelihood improved by 1.28 
## Iteration 44 of at most 100: 
## The log-likelihood improved by 1.285 
## Iteration 45 of at most 100: 
## The log-likelihood improved by 1.365 
## Iteration 46 of at most 100: 
## The log-likelihood improved by 1.145 
## Iteration 47 of at most 100: 
## The log-likelihood improved by 1.337 
## Iteration 48 of at most 100: 
## The log-likelihood improved by 1.132 
## Iteration 49 of at most 100: 
## The log-likelihood improved by 1.104 
## Iteration 50 of at most 100: 
## The log-likelihood improved by 1.105 
## Iteration 51 of at most 100: 
## The log-likelihood improved by 1.264 
## Iteration 52 of at most 100: 
## The log-likelihood improved by 1.471 
## Iteration 53 of at most 100: 
## The log-likelihood improved by 1.166 
## Iteration 54 of at most 100: 
## The log-likelihood improved by 1.106 
## Iteration 55 of at most 100: 
## The log-likelihood improved by 1.142 
## Iteration 56 of at most 100: 
## The log-likelihood improved by 1.114 
## Iteration 57 of at most 100: 
## The log-likelihood improved by 1.386 
## Iteration 58 of at most 100: 
## The log-likelihood improved by 1.357 
## Iteration 59 of at most 100: 
## The log-likelihood improved by 1.291 
## Iteration 60 of at most 100: 
## The log-likelihood improved by 1.23 
## Iteration 61 of at most 100: 
## The log-likelihood improved by 1.205 
## Iteration 62 of at most 100: 
## The log-likelihood improved by 1.3 
## Iteration 63 of at most 100: 
## The log-likelihood improved by 1.313 
## Iteration 64 of at most 100: 
## The log-likelihood improved by 1.562 
## Iteration 65 of at most 100: 
## The log-likelihood improved by 1.323 
## Iteration 66 of at most 100: 
## The log-likelihood improved by 1.23 
## Iteration 67 of at most 100: 
## The log-likelihood improved by 1.392 
## Iteration 68 of at most 100: 
## The log-likelihood improved by 1.47 
## Iteration 69 of at most 100: 
## The log-likelihood improved by 1.546 
## Iteration 70 of at most 100: 
## The log-likelihood improved by 1.587 
## Iteration 71 of at most 100: 
## The log-likelihood improved by 1.519 
## Iteration 72 of at most 100: 
## The log-likelihood improved by 1.612 
## Iteration 73 of at most 100: 
## The log-likelihood improved by 1.508 
## Iteration 74 of at most 100: 
## The log-likelihood improved by 1.588 
## Iteration 75 of at most 100: 
## The log-likelihood improved by 1.787 
## Iteration 76 of at most 100: 
## The log-likelihood improved by 1.604 
## Iteration 77 of at most 100: 
## The log-likelihood improved by 1.757 
## Iteration 78 of at most 100: 
## The log-likelihood improved by 1.69 
## Iteration 79 of at most 100: 
## The log-likelihood improved by 1.359 
## Iteration 80 of at most 100: 
## The log-likelihood improved by 1.29 
## Iteration 81 of at most 100: 
## The log-likelihood improved by 1.121 
## Iteration 82 of at most 100: 
## The log-likelihood improved by 0.926 
## Iteration 83 of at most 100: 
## The log-likelihood improved by 1.385 
## Iteration 84 of at most 100: 
## The log-likelihood improved by 1.282 
## Iteration 85 of at most 100: 
## The log-likelihood improved by 1.409 
## Iteration 86 of at most 100: 
## The log-likelihood improved by 1.505 
## Iteration 87 of at most 100: 
## The log-likelihood improved by 1.339 
## Iteration 88 of at most 100: 
## The log-likelihood improved by 1.285 
## Iteration 89 of at most 100: 
## The log-likelihood improved by 1.328 
## Iteration 90 of at most 100: 
## The log-likelihood improved by 1.537 
## Iteration 91 of at most 100: 
## The log-likelihood improved by 1.406 
## Iteration 92 of at most 100: 
## The log-likelihood improved by 1.204 
## Iteration 93 of at most 100: 
## The log-likelihood improved by 1.168 
## Iteration 94 of at most 100: 
## The log-likelihood improved by 1.258 
## Iteration 95 of at most 100: 
## The log-likelihood improved by 1.276 
## Iteration 96 of at most 100: 
## The log-likelihood improved by 1.076 
## Iteration 97 of at most 100: 
## The log-likelihood improved by 1.22 
## Iteration 98 of at most 100: 
## The log-likelihood improved by 1.255 
## Iteration 99 of at most 100: 
## The log-likelihood improved by 1.359 
## Iteration 100 of at most 100: 
## The log-likelihood improved by 1.252
```

```
## MCMLE estimation did not converge after 100 iterations. The estimated coefficients may not be accurate. Estimation may be resumed by passing the coefficients as initial values; see 'init' under ?control.ergm for details.
```

```
## Evaluating log-likelihood at the estimate. Using 20 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 .
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
Sys.time() - t0
```

```
## Time difference of 6.911152 hours
```

```r
save(my.ergm.fit, file = './Rdata/my.ergm.fit.rda')
```

```r
anova.ergm(my.ergm.fit)
```

```
## Analysis of Variance Table
## 
## Model 1: mnet3.simple ~ edges + kstar(2) + kstar(3) + triangle
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                   1604736          0                 
## Model 1:  4  -723636   1604732     723636    < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(my.ergm.fit)
```

```
## 
## ==========================
## Summary of model fit
## ==========================
## 
## Formula:   mnet3.simple ~ edges + kstar(2) + kstar(3) + triangle
## 
## Iterations:  100 out of 100 
## 
## Monte Carlo MLE Results:
##            Estimate Std. Error MCMC % p-value    
## edges    -4.202e+00  3.779e-03    100  <1e-04 ***
## kstar2   -7.106e-03  2.329e-04      3  <1e-04 ***
## kstar3    1.943e-06  1.600e-06      2   0.224    
## triangle  6.673e-02  2.672e-04      2  <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2224636  on 1604736  degrees of freedom
##  Residual Deviance:  723636  on 1604732  degrees of freedom
##  
## AIC: 723644    BIC: 723694    (Smaller is better.)
```

```r
t0 <- Sys.time()
diag <- gof(my.ergm.fit, control=control.gof.ergm(parallel=cores, parallel.type="PSOCK"))
save(diag, file = './Rdata/diag.my.erggm.fit.rda')
Sys.time() - t0
```

```
## Time difference of 16.864 hours
```

```r
mcmc.diagnostics(my.ergm.fit)
```

```
## Sample statistics summary:
## 
## Iterations = 20000:21000
## Thinning interval = 500 
## Number of chains = 48 
## Sample size per chain = 3 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##               Mean        SD  Naive SE Time-series SE
## edges    6.989e+02 9.318e+01 7.765e+00         0.5325
## kstar2   2.843e+05 1.958e+04 1.632e+03        57.2844
## kstar3   3.096e+07 2.483e+06 2.069e+05      7950.0680
## triangle 1.143e+05 6.080e+03 5.067e+02         4.7598
## 
## 2. Quantiles for each variable:
## 
##               2.5%       25%      50%       75%     97.5%
## edges    5.083e+02 6.525e+02      693 7.762e+02 8.577e+02
## kstar2   2.493e+05 2.679e+05   288960 2.990e+05 3.155e+05
## kstar3   2.661e+07 2.882e+07 31231074 3.273e+07 3.461e+07
## triangle 1.027e+05 1.101e+05   113597 1.170e+05 1.287e+05
## 
## 
## Sample statistics cross-correlations:
##              edges    kstar2    kstar3  triangle
## edges    1.0000000 0.6444450 0.5538330 0.3354276
## kstar2   0.6444450 1.0000000 0.9615471 0.7415277
## kstar3   0.5538330 0.9615471 1.0000000 0.6461732
## triangle 0.3354276 0.7415277 0.6461732 1.0000000
## 
## Sample statistics auto-correlation:
## Chain 1 
##               edges       kstar2      kstar3  triangle
## Lag 0     1.0000000  1.000000000  1.00000000  1.000000
## Lag 500  -0.1666667 -0.009483473 -0.03531413 -0.202722
## Lag 1000 -0.3333333 -0.490516527 -0.46468587 -0.297278
## Chain 2 
##               edges     kstar2     kstar3     triangle
## Lag 0     1.0000000  1.0000000  1.0000000  1.000000000
## Lag 500  -0.2772657 -0.2685914 -0.1128228 -0.004504505
## Lag 1000 -0.2227343 -0.2314086 -0.3871772 -0.495495495
## Chain 3 
##               edges      kstar2      kstar3   triangle
## Lag 0     1.0000000  1.00000000  1.00000000  1.0000000
## Lag 500  -0.3809524 -0.02932413 -0.00368935 -0.1630172
## Lag 1000 -0.1190476 -0.47067587 -0.49631065 -0.3369828
## Chain 4 
##                edges   kstar2    kstar3    triangle
## Lag 0     1.00000000  1.00000  1.000000  1.00000000
## Lag 500  -0.05128205 -0.05963 -0.113655 -0.09689922
## Lag 1000 -0.44871795 -0.44037 -0.386345 -0.40310078
## Chain 5 
##               edges      kstar2      kstar3    triangle
## Lag 0     1.0000000  1.00000000  1.00000000  1.00000000
## Lag 500  -0.6282051 -0.04724229 -0.01059809 -0.05128205
## Lag 1000  0.1282051 -0.45275771 -0.48940191 -0.44871795
## Chain 6 
##                 edges       kstar2       kstar3   triangle
## Lag 0     1.000000000  1.000000000  1.000000000  1.0000000
## Lag 500  -0.001539646 -0.000639927 -0.001997266 -0.1666667
## Lag 1000 -0.498460354 -0.499360073 -0.498002734 -0.3333333
## Chain 7 
##                edges      kstar2     kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.02380952 -0.01956081 -0.0145115 -0.3809524
## Lag 1000 -0.47619048 -0.48043919 -0.4854885 -0.1190476
## Chain 8 
##               edges       kstar2      kstar3   triangle
## Lag 0     1.0000000  1.000000000  1.00000000  1.0000000
## Lag 500  -0.2857143 -0.006667532 -0.01574013 -0.3205128
## Lag 1000 -0.2142857 -0.493332468 -0.48425987 -0.1794872
## Chain 9 
##          edges     kstar2     kstar3   triangle
## Lag 0      1.0  1.0000000  1.0000000  1.0000000
## Lag 500   -0.5 -0.2597705 -0.3073028 -0.1891792
## Lag 1000   0.0 -0.2402295 -0.1926972 -0.3108208
## Chain 10 
##                edges        kstar2      kstar3   triangle
## Lag 0     1.00000000  1.0000000000  1.00000000  1.0000000
## Lag 500  -0.04295533 -0.0003823233 -0.00939743 -0.1666667
## Lag 1000 -0.45704467 -0.4996176767 -0.49060257 -0.3333333
## Chain 11 
##                edges     kstar2      kstar3    triangle
## Lag 0     1.00000000  1.0000000  1.00000000  1.00000000
## Lag 500  -0.42982456 -0.1284844 -0.09895445 -0.00877193
## Lag 1000 -0.07017544 -0.3715156 -0.40104555 -0.49122807
## Chain 12 
##               edges     kstar2      kstar3    triangle
## Lag 0     1.0000000  1.0000000  1.00000000  1.00000000
## Lag 500  -0.1666667 -0.1482269 -0.03003165 -0.02380952
## Lag 1000 -0.3333333 -0.3517731 -0.46996835 -0.47619048
## Chain 13 
##                edges      kstar2     kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.05128205 -0.02806687 -0.1219269 -0.1580486
## Lag 1000 -0.44871795 -0.47193313 -0.3780731 -0.3419514
## Chain 14 
##               edges     kstar2      kstar3   triangle
## Lag 0     1.0000000  1.0000000  1.00000000  1.0000000
## Lag 500  -0.3522562 -0.6663531 -0.56507513 -0.6637037
## Lag 1000 -0.1477438  0.1663531  0.06507513  0.1637037
## Chain 15 
##               edges     kstar2      kstar3    triangle
## Lag 0     1.0000000  1.0000000  1.00000000  1.00000000
## Lag 500  -0.2857143 -0.6648968 -0.47773616 -0.03846154
## Lag 1000 -0.2142857  0.1648968 -0.02226384 -0.46153846
## Chain 16 
##               edges     kstar2       kstar3  triangle
## Lag 0     1.0000000  1.0000000  1.000000000  1.000000
## Lag 500  -0.3324125 -0.2409434 -0.008187184 -0.660853
## Lag 1000 -0.1675875 -0.2590566 -0.491812816  0.160853
## Chain 17 
##                edges      kstar2    kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.000000  1.0000000
## Lag 500  -0.56140351 -0.54426188 -0.629861 -0.2734628
## Lag 1000  0.06140351  0.04426188  0.129861 -0.2265372
## Chain 18 
##                edges      kstar2     kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.08164642 -0.07865448 -0.1403255 -0.1561451
## Lag 1000 -0.41835358 -0.42134552 -0.3596745 -0.3438549
## Chain 19 
##               edges      kstar2     kstar3   triangle
## Lag 0     1.0000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.1052632 -0.47818203 -0.6660067 -0.3809524
## Lag 1000 -0.3947368 -0.02181797  0.1660067 -0.1190476
## Chain 20 
##               edges      kstar2     kstar3   triangle
## Lag 0     1.0000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.1666667 -0.00246844 -0.6110414 -0.2734628
## Lag 1000 -0.3333333 -0.49753156  0.1110414 -0.2265372
## Chain 21 
##                edges       kstar2      kstar3   triangle
## Lag 0     1.00000000  1.000000000  1.00000000  1.0000000
## Lag 500  -0.06218905 -0.007680692 -0.07716527 -0.1769364
## Lag 1000 -0.43781095 -0.492319308 -0.42283473 -0.3230636
## Chain 22 
##               edges     kstar2      kstar3   triangle
## Lag 0     1.0000000  1.0000000  1.00000000  1.0000000
## Lag 500  -0.0356623 -0.3459887 -0.46178818 -0.1531633
## Lag 1000 -0.4643377 -0.1540113 -0.03821182 -0.3468367
## Chain 23 
##                edges     kstar2       kstar3   triangle
## Lag 0     1.00000000  1.0000000  1.000000000  1.0000000
## Lag 500  -0.48003072 -0.3868914 -0.009529565 -0.3809524
## Lag 1000 -0.01996928 -0.1131086 -0.490470435 -0.1190476
## Chain 24 
##                edges      kstar2      kstar3 triangle
## Lag 0     1.00000000  1.00000000  1.00000000      NaN
## Lag 500  -0.07142857 -0.04854956 -0.40459011      NaN
## Lag 1000 -0.42857143 -0.45145044 -0.09540989      NaN
## Chain 25 
##                edges      kstar2        kstar3    triangle
## Lag 0     1.00000000  1.00000000  1.0000000000  1.00000000
## Lag 500  -0.08998935 -0.01994856 -0.0002863956 -0.06218905
## Lag 1000 -0.41001065 -0.48005144 -0.4997136044 -0.43781095
## Chain 26 
##               edges      kstar2      kstar3   triangle
## Lag 0     1.0000000  1.00000000  1.00000000  1.0000000
## Lag 500  -0.1666667 -0.04484054 -0.08435006 -0.6666667
## Lag 1000 -0.3333333 -0.45515946 -0.41564994  0.1666667
## Chain 27 
##                edges     kstar2     kstar3   triangle
## Lag 0     1.00000000  1.0000000  1.0000000  1.0000000
## Lag 500  -0.55198777 -0.6538759 -0.6625243 -0.2073974
## Lag 1000  0.05198777  0.1538759  0.1625243 -0.2926026
## Chain 28 
##                edges     kstar2    kstar3    triangle
## Lag 0     1.00000000  1.0000000  1.000000  1.00000000
## Lag 500  -0.01079622 -0.2422754 -0.374068 -0.57482993
## Lag 1000 -0.48920378 -0.2577246 -0.125932  0.07482993
## Chain 29 
##               edges     kstar2     kstar3   triangle
## Lag 0     1.0000000  1.0000000  1.0000000  1.0000000
## Lag 500  -0.0356623 -0.2592007 -0.2677286 -0.2283105
## Lag 1000 -0.4643377 -0.2407993 -0.2322714 -0.2716895
## Chain 30 
##                edges     kstar2      kstar3   triangle
## Lag 0     1.00000000  1.0000000  1.00000000  1.0000000
## Lag 500  -0.05128205 -0.1626669 -0.02204467 -0.3205128
## Lag 1000 -0.44871795 -0.3373331 -0.47795533 -0.1794872
## Chain 31 
##                edges      kstar2        kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000000  1.0000000
## Lag 500  -0.07142857 -0.01463643 -0.0006680053 -0.1565358
## Lag 1000 -0.42857143 -0.48536357 -0.4993319947 -0.3434642
## Chain 32 
##                edges      kstar2     kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.05128205 -0.59426603 -0.6662811 -0.1666667
## Lag 1000 -0.44871795  0.09426603  0.1662811 -0.3333333
## Chain 33 
##                edges     kstar2      kstar3 triangle
## Lag 0     1.00000000  1.0000000  1.00000000      1.0
## Lag 500  -0.02885748 -0.1181639 -0.08291295      0.0
## Lag 1000 -0.47114252 -0.3818361 -0.41708705     -0.5
## Chain 34 
##                edges     kstar2      kstar3   triangle
## Lag 0     1.00000000  1.0000000  1.00000000  1.0000000
## Lag 500  -0.07142857 -0.6655455 -0.55181068 -0.1118721
## Lag 1000 -0.42857143  0.1655455  0.05181068 -0.3881279
## Chain 35 
##          edges     kstar2      kstar3    triangle
## Lag 0      1.0  1.0000000  1.00000000  1.00000000
## Lag 500   -0.5 -0.2022972 -0.09040076 -0.05128205
## Lag 1000   0.0 -0.2977028 -0.40959924 -0.44871795
## Chain 36 
##               edges     kstar2     kstar3   triangle
## Lag 0     1.0000000  1.0000000  1.0000000  1.0000000
## Lag 500  -0.3809524 -0.6230346 -0.3490611 -0.1666667
## Lag 1000 -0.1190476  0.1230346 -0.1509389 -0.3333333
## Chain 37 
##               edges      kstar2      kstar3   triangle
## Lag 0     1.0000000  1.00000000  1.00000000  1.0000000
## Lag 500  -0.6639344 -0.59929088 -0.48231784 -0.1429137
## Lag 1000  0.1639344  0.09929088 -0.01768216 -0.3570863
## Chain 38 
##                edges     kstar2     kstar3    triangle
## Lag 0     1.00000000  1.0000000  1.0000000  1.00000000
## Lag 500  -0.05128205 -0.2806025 -0.3308303 -0.08602151
## Lag 1000 -0.44871795 -0.2193975 -0.1691697 -0.41397849
## Chain 39 
##                edges     kstar2     kstar3 triangle
## Lag 0     1.00000000  1.0000000  1.0000000      1.0
## Lag 500  -0.42982456 -0.1445191 -0.1061602      0.0
## Lag 1000 -0.07017544 -0.3554809 -0.3938398     -0.5
## Chain 40 
##                edges      kstar2     kstar3   triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.0000000
## Lag 500  -0.07928803 -0.55765081 -0.6665222 -0.1626347
## Lag 1000 -0.42071197  0.05765081  0.1665222 -0.3373653
## Chain 41 
##                edges     kstar2      kstar3    triangle
## Lag 0     1.00000000  1.0000000  1.00000000  1.00000000
## Lag 500  -0.04295533 -0.2030979 -0.54357875 -0.07142857
## Lag 1000 -0.45704467 -0.2969021  0.04357875 -0.42857143
## Chain 42 
##               edges      kstar2     kstar3    triangle
## Lag 0     1.0000000  1.00000000  1.0000000  1.00000000
## Lag 500  -0.1118721 -0.08886146 -0.1136559 -0.00877193
## Lag 1000 -0.3881279 -0.41113854 -0.3863441 -0.49122807
## Chain 43 
##               edges     kstar2     kstar3   triangle
## Lag 0     1.0000000  1.0000000  1.0000000  1.0000000
## Lag 500  -0.6621622 -0.6665555 -0.6241077 -0.2634409
## Lag 1000  0.1621622  0.1665555  0.1241077 -0.2365591
## Chain 44 
##                edges       kstar2      kstar3    triangle
## Lag 0     1.00000000  1.000000000  1.00000000  1.00000000
## Lag 500  -0.01360544 -0.001675364 -0.05293682 -0.02985075
## Lag 1000 -0.48639456 -0.498324636 -0.44706318 -0.47014925
## Chain 45 
##                edges      kstar2    kstar3    triangle
## Lag 0     1.00000000  1.00000000  1.000000  1.00000000
## Lag 500  -0.45045045 -0.42599333 -0.123287 -0.01612903
## Lag 1000 -0.04954955 -0.07400667 -0.376713 -0.48387097
## Chain 46 
##                edges      kstar2     kstar3    triangle
## Lag 0     1.00000000  1.00000000  1.0000000  1.00000000
## Lag 500  -0.00877193 -0.05015295 -0.1577876 -0.09689922
## Lag 1000 -0.49122807 -0.44984705 -0.3422124 -0.40310078
## Chain 47 
##                 edges      kstar2     kstar3    triangle
## Lag 0     1.000000000  1.00000000  1.0000000  1.00000000
## Lag 500  -0.004504505 -0.04784212 -0.3010526 -0.05128205
## Lag 1000 -0.495495495 -0.45215788 -0.1989474 -0.44871795
## Chain 48 
##          edges    kstar2     kstar3   triangle
## Lag 0      1.0  1.000000  1.0000000  1.0000000
## Lag 500    0.0 -0.295326 -0.1556696 -0.1747291
## Lag 1000  -0.5 -0.204674 -0.3443304 -0.3252709
## 
## Sample statistics burn-in diagnostic (Geweke):
```

```
## Error in ar.yw.default(x, aic = aic, order.max = order.max, na.action = na.action, : 'order.max' must be >= 1
```

```r
par(mfrow=c(2, 2))
plot(diag)
```

![plot of chunk viz-diag](figure/viz-diag-1.png)

