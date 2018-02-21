# Plotly Map
library(plotly)
library(plotrix)
packageVersion('plotly')

# Flight Paths Map
library(dplyr)
# airport locations
#air <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')
# flights between airports
# flights <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')
# flights$id <- seq_len(nrow(flights))
complete<-edges[which(is.na(edges$long_target) != T & is.na(edges$lat_target) != T & 
                        is.na(edges$long_source) != T & is.na(edges$lat_source) != T),]
auth_data2<-auth_data[which(is.na(auth_data$long)==F & is.na(auth_data$lat)==F),]

complete$id<-seq_len(nrow(complete))
infoAuth <- paste(c('<b>'),auth_data2$name,c('</b>'),
                  c('<br>contributions: '),auth_data2$numPub,c(' publications'),
                  c('<br>Times cited: '),auth_data2$timesCited,
                  c('<br>centre: '),auth_data2$affiliation,
                  c('<br>city: '),auth_data2$city,
                  c('<br>country: '),auth_data2$country, sep = '')
infoEdge <- paste(paste(complete$source, complete$target, sep = ' -- '),
                  c('<br>timesCited: '), complete$timesCited,sep = '')
complete$infoEdge <- infoEdge
auth_data2$infoAuth <- infoAuth

# map projection
geo <- list(
  # scope = 'north america',
  projection = list(type = 'Mercator'),
  showland = TRUE,
  showframe = FALSE,
  showcoastlines = TRUE,
  landcolor = toRGB("gray90"),
  # coastlinecolor = toRGB("black"),
  showcountries = F,
  bgcolor = toRGB("white", alpha = 0)
  ,countrycolor = toRGB("gray80")
)

p <- plot_geo(
    # locationmode = 'USA-states', 
    # locationmode = 'world', 
    color = I("blue")) %>%
  add_markers(
    # data = air, x = ~long, y = ~lat, text = ~airport,
    # size = ~cnt, hoverinfo = "text", alpha = 0.5
    data = auth_data2, x = ~long, y = ~lat, text = ~infoAuth,
    size = ~degree, hoverinfo = "text", alpha = 0.8
  ) %>%
  add_segments(
    # data = group_by(flights, id),
    # x = ~start_lon, xend = ~end_lon,
    # y = ~start_lat, yend = ~end_lat,
    # alpha = 0.3, size = I(1), hoverinfo = "none"
    data = group_by(complete, id),
    x = ~long_source, xend = ~long_target,
    y = ~lat_source, yend = ~lat_target, text = ~infoEdge,
    alpha = 0.07, size = I(0.3), color = I("pink"), hoverinfo = "text"
  ) %>%
  layout(
    title = 'Malaria Research Collaboration Map <br>(Hover for author names)',
    geo = geo, showlegend = FALSE, height=800
  )

plotly('rosericazondekon','E3lmpmz5qKTNRn70CbkH')
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="coauthorshipMalariaBenin")
chart_link
