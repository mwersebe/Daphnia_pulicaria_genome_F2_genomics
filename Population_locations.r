#Map For Sample Locations: 
#Matthew Wersebe, University of Oklahoma
#June 7th, 2021 
###############################################################################
library(ggplot2)
theme_set(theme_light())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(rgeos)
library(ggspatial)
library(maps)
library(tools)

#Population locations and proportions of animals isolated. 
locations = read.table("/home/weider/LinkageMap_Genomics/population_locations.txt", header =T, sep = "\t")
head(locations)

#USA and Canada shapes with state boundaries:
Countries = ne_states(country = c('united states of america', 'canada'), returnclass = "sf")

ggplot(data = Countries)+
  geom_sf(fill = "antiquewhite1")+
  coord_sf(ylim = c(41, 56), xlim = c(-116.2, -78.5))+
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl", pad_y = unit(0.36, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(0.75, "cm"), height = unit(0.5, "cm"), width = unit(0.5, "cm"))+
  ggtitle("Sample Locations")+
  geom_point(data = locations, aes(x = Long, y = Lat, color= Type), size = 0.8)+
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
  