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

pop_locations <- ggplot(data = Countries)+
  geom_sf(fill = "antiquewhite1")+
  coord_sf(ylim = c(41, 56), xlim = c(-116.2, -78.5))+
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl", pad_y = unit(0.36, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(0.75, "cm"), height = unit(0.5, "cm"), width = unit(0.5, "cm"))+
  ggtitle("Sample Locations")+
  geom_point(data = locations, aes(x = Long, y = Lat, color= Type), size = 3.0)+
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))+
  theme(text = element_text(size = 20)) 
pop_locations


###############################################################################\
#pca plot:
# PCA:

daphnia.pca <- glPca(single_snp.genlight, nf = 3)

#PCA Plots:
#Eigenvalues
barplot(100*daphnia.pca$eig/sum(daphnia.pca$eig), col = heat.colors(180), main = "PCA Eigenvalues")
title(ylab="Percent Variance\nexplained", line =2)
title(xlab="Eignevalues", line = 1)


Variance <- 100*daphnia.pca$eig/sum(daphnia.pca$eig)
head(Variance)
#Scatter 
pca.scores <- as.data.frame(daphnia.pca$scores)
pca.scores$pop <- pop(single_snp.genlight)

daphnia.plot <- ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop))+
  geom_point(size = 2)+
  stat_ellipse(level = 0.95, size = 1)+
  scale_color_manual(labels = c("Pulex", "Hybrid", "Pulicaria"), values = cols)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  labs (y = "PC2 (5.17%)", x = "PC1 (26.18%)")+
  theme(text = element_text(size = 20), legend.title = element_blank())
daphnia.plot

###############################################################################
# faststructure barplot SVG:

library(grImport2)

readPicture("/home/weider/LinkageMap_Genomics/fastStructure/Pulex-Pulicaria-struct_edit.svg")


