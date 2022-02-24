library(sf)
library(ggplot2)
ex<-st_read("GSIM/GSIM_metadata/GSIM_catchments/au_0002256.SHP")

ggplot() + 
  geom_sf(data = ex, size = 3, color = "black", fill = "orange") + 
  ggtitle(" ") + 
  coord_sf()
