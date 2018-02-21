library(plotly)
library(igraph)
library(plotrix)

#data(karate, package="igraphdata")
#G <- upgrade_graph(karate)
G<-mnet3
# L <- layout.kamada.kawai(G)
L <- l.mnet3

vs <- V(G)
es <- as.data.frame(get.edgelist(G))

# Create vertices and edges
Nv <- length(vs)
Ne <- length(es[1]$V1)

# Create Nodes
Xn <- L[,1]
Yn <- L[,2]

nodeColors<-color.scale(V(mnet3)$community,c(0,1,1),c(1,1,0),0)

network <- plot_ly(x = ~Xn, y = ~Yn, color = nodeColors, size = 1+1 * sqrt(bw.mnet/1000), mode = "markers", text = vs$name, hoverinfo = "text")

# Create Edges
edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = 0.01),
    x0 = Xn[v0],
    y0 = Yn[v0],
    x1 = Xn[v1],
    y1 = Yn[v1]
  )
  
  edge_shapes[[i]] <- edge_shape
}

axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

p <- layout(
  network,
  title = 'Malaria Co-authorship Network in Benin',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis
)
hide_legend(p)

plotly('rosericazondekon','E3lmpmz5qKTNRn70CbkH')
# Create Network
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="coauthorship-network-r")
chart_link
