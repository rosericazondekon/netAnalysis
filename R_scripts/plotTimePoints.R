# library(igraph)
# load("Rdata/netByTimePoints.rda")
par(mfrow=c(2,4), cex.main=3,xpd=NA,oma=c(1,0,1.8,0),mar=c(1,0,1.8,0),xpd=NA)
# for(i in length(net2)){
#   nt<-igraph::delete_vertices(net2[[1]], which(V(net2[[i]])$degree==0))
#   plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.01)
# }
#delete.vertices(simplify(g), degree(g)==0)

colr="steelblue"
s=10

nt<-igraph::delete_vertices(net2[[1]], which(V(net2[[1]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '1996-2006',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[2]], which(V(net2[[2]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2007-2009',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[3]], which(V(net2[[3]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2010-2011',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[4]], which(V(net2[[4]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2012-2013',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[5]], which(V(net2[[5]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2014',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[6]], which(V(net2[[6]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2015',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)

nt<-igraph::delete_vertices(net2[[7]], which(V(net2[[7]])$degree==0))
plot(nt,layout =layout.kamada.kawai(nt),vertex.label=NA,edge.width = 0.5,main = '2016',
     vertex.size=7+5 * sqrt(estimate_betweenness(nt,cutoff = 0)/1000),vertex.color=colr)