#South Center Marey Map Plots
#Matthew Wersebe 3-14-2021
#South Center 13B LOD 26 Linkage Map.
######################################
library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd("/home/matt/Desktop/Linkage_plots")

#Chromosome One:
CHR1 <- read.table("SC_map_Chrom01-marey.txt", header=T)

#make it pretty with ggplot, Matthew.
#marey map
chrom1 <- ggplot(CHR1, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 1 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+
  theme_light()
chrom1
#run window recomb rate estimates
Window_size = 1000000
Step_size = 10000
recombination_1 <- Recomb_rate(Data = CHR1, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_1 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 1 Recombination")+
  theme_light()+ theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
rec_rate_1
#STITCH TOGETHER
#ggarrange(chrom1, rec_rate_1, labels = c("A) ", "B) "), common.legend =T, widths = c(5,5))

  ###########################################################################
#Chromosome #2
CHR2 <- read.table("SC_map_Chrom02-marey.txt", header=T)
chrom2 <- ggplot(CHR2, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 2 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+ theme_light()
chrom2
#run window recomb rate estimates

recombination_2 <- Recomb_rate(Data = CHR2, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_2 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 2Recombination")+
  theme_light()+ theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
#STITCH TOGETHER
#ggarrange(chrom2, rec_rate_2, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))

##################################################################################
#Chromosome 3
CHR03 <- read.table("SC_map_Chrom03-marey.txt", header=T)

chrom3 <- ggplot(CHR03, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 3 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+theme_light()

#run window recomb rate estimates

recombination_3 <- Recomb_rate(Data = CHR03, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_3 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 3 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light()+ theme(axis.title.x=element_blank(),
                                                                                            axis.text.x=element_blank(),
                                                                                            axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom3, rec_rate_3, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
#######################################################################################################
#Chromosome #4
CHR04 <- read.table("SC_map_Chrom04-marey.txt", header=T)

chrom4 <- ggplot(CHR04, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 4 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+theme_light()

#run window recomb rate estimates

recombination_4 <- Recomb_rate(Data = CHR04, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_4 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 4 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light()+ theme(axis.title.x=element_blank(),
                                                                                            axis.text.x=element_blank(),
                                                                                            axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom4, rec_rate_4, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
#####################################################################################
#Chromosome 5
CHR05 <- read.table("SC_map_Chrom05-marey.txt", header=T)

chrom5 <- ggplot(CHR05, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 5 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+theme_light()

#run window recomb rate estimates

recombination_5  <- Recomb_rate(Data = CHR05, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_5 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 5 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light() + theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom5, rec_rate_5, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
####################################################################################################
#Chromosome 6 
CHR06 <- read.table("SC_map_Chrom06-marey.txt", header=T)

chrom6 <- ggplot(CHR06, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 6 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+ theme_light()

#run window recomb rate estimates

recombination_6 <- Recomb_rate(Data = CHR06, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_6 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 6 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light() + theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom6, rec_rate_6, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))

#Chromosome 7
CHR07 <- read.table("SC_map_Chrom07-marey.txt", header=T)

chrom7 <- ggplot(CHR07, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 7 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)")+theme_light()

#run window recomb rate estimates

recombination_7 <- Recomb_rate(Data = CHR07, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_7 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 7 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+ theme_light() + theme(axis.title.x=element_blank(),
                                                                                              axis.text.x=element_blank(),
                                                                                              axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom7, rec_rate_7, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
#######################################################################################################
#Chromosomr 8
CHR08 <- read.table("SC_map_Chrom08-marey.txt", header=T)

chrom8 <- ggplot(CHR08, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 8 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)") + theme_light()

#run window recomb rate estimates

recombination_8 <- Recomb_rate(Data = CHR08, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_8 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 8 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male")) +theme_light() +theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom8, rec_rate_8, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
###################################################################################################################
#Chromosome 9
CHR09 <- read.table("SC_map_Chrom09-marey.txt", header=T)

chrom9 <- ggplot(CHR09, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 9 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)") +theme_light()

#run window recomb rate estimates

recombination_9 <- Recomb_rate(Data = CHR09, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_9 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 9 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+ theme_light() + theme(axis.title.x=element_blank(),
                                                                                              axis.text.x=element_blank(),
                                                                                              axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom9, rec_rate_9, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
####################################################################################################
#Chromosome 10
CHR10 <- read.table("SC_map_Chrom10-marey.txt", header=T)

chrom10 <- ggplot(CHR10, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 10 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)") +theme_light()

#run window recomb rate estimates

recombination_10 <- Recomb_rate(Data = CHR10, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_10 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 10 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light() + theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom10, rec_rate_10, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
###################################################################################################
#Chromosome 11
CHR11 <- read.table("SC_map_Chrom11-marey.txt", header=T)

chrom11 <- ggplot(CHR11, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 11 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)") + theme_light()

#run window recomb rate estimates

recombination_11 <- Recomb_rate(Data = CHR11, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_11 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 11 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male")) +theme_light()+ theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom11, rec_rate_11, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))

#Chromosome 12
CHR12 <- read.table("SC_map_Chrom12-marey.txt", header=T)

chrom12 <- ggplot(CHR12, aes(x=Physical_position)) + 
  geom_point(aes(y=cM_male, color= "steelblue"))+
  geom_point(aes(y=cM_female, color= "darkred"))+
  ggtitle("LG 12 Marey Map")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+
  xlab("Phyiscal Position (MB)")+ ylab("Genetic Position (cM)") +theme_light()

#run window recomb rate estimates

recombination_12 <- Recomb_rate(Data = CHR12, Window_Size = Window_size, Step_Size = Step_size)
#Plot
rec_rate_12 <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ggtitle("LG 12 Recombination")+
  xlab("Phyiscal Position (MB)")+ ylab("R (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))+theme_light() + theme(axis.title.x=element_blank(),
                                                                                             axis.text.x=element_blank(),
                                                                                             axis.ticks.x=element_blank())

#STITCH TOGETHER
#ggarrange(chrom12, rec_rate_12, labels = c("A)", "B)"), common.legend =T, widths = c(5,5))
##################################################################################################
# Marey Map Stack

ggarrange(chrom1, chrom2, chrom3, chrom4, chrom5, chrom6, chrom7, chrom8, chrom9, chrom10, chrom11, chrom12, ncol = 3, nrow = 4, legend = "bottom", common.legend = T)


############################################################################
#Linkage Map Plots:
# 3-17-2021
#SC 13B LOD26 Linkage Map
library(LinkageMapView)
setwd("/home/matt/Desktop/Linkage_plots")
male_map <- read.table("SC_map_13BLOD26-mapviewermale.txt", header= T, )
female_map <- read.table("SC_map_13BLOD26-mapviewerfemale.txt", header= T, )
#male map
outfile = file.path("/home/matt/Desktop/Linkage_plots", "SC_13B_LOD26_male.pdf")

maxpos <- floor(max(male_map$position[male_map$Group == "CHR01"]))
at.axis <- seq(0, maxpos)


axlab <- vector() 

for (lab in 0:maxpos) {if (!lab %% 10){axlab <-c(axlab, lab)} else {axlab <-c(axlab, NA)}}

lmv.linkage.plot(male_map, outfile, mapthese = c("CHR01", "CHR02", "CHR03", "CHR04", "CHR05", "CHR06", "CHR07", "CHR08", "CHR09", "CHR10", "CHR11", "CHR12"), dupnbr=TRUE, cex.axis = 1, at.axis = at.axis, labels.axis = axlab)

#female map
maxpos <- floor(max(female_map$position[male_map$Group == "CHR01"]))
at.axis <- seq(0, maxpos)


axlab <- vector() 

for (lab in 0:maxpos) {if (!lab %% 10){axlab <-c(axlab, lab)} else {axlab <-c(axlab, NA)}}


outfile_f = file.path("/home/matt/Desktop/Linkage_plots", "SC_13B_LOD26_female.pdf")
lmv.linkage.plot(female_map, outfile_f, mapthese = c("CHR01", "CHR02", "CHR03", "CHR04", "CHR05", "CHR06", "CHR07", "CHR08", "CHR09", "CHR10", "CHR11", "CHR12"), dupnbr =TRUE, cex.axis = 1, at.axis = at.axis, labels.axis = axlab)

###########################################################################################
#Relationship of estimated phyiscal and genetic lengths of chromosomes

correlation <- read.table("correlation.txt", header=F)
head(correlation)

sum(correlation$V4)

sum(correlation$V5)

corr.plot_m <-ggplot(correlation, aes(x=V2, y=V4))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("Phyiscal Length (MB)")+ ylab("Genetic Length (cM)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, data = correlation)+theme_light()
corr.plot_m

corr.plot_f <-ggplot(correlation, aes(x=V2, y=V5))+
  geom_point(color=  "darkred")+
  ggtitle("Female Map")+
  xlab("Phyiscal Length (MB)")+ ylab("Genetic Length (cM)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, data = correlation)+theme_light()
corr.plot_f


#Recombination Rates

male_recomb_rate <- (correlation$V4/(correlation$V2/1000000))
male_recomb_rate

female_recomb_rates <- (correlation$V5/(correlation$V2/1000000))
female_recomb_rates

correlation <- cbind.data.frame(correlation, male_recomb_rate)
head(correlation)

correlation <- cbind.data.frame(correlation, female_recomb_rates)
head(correlation)

mean(correlation$male_recomb_rate)
mean(correlation$female_recomb_rates)

genomewiderecombratemale <- (sum(correlation$V4)/(sum(correlation$V2)/1000000))
genomewiderecombratemale

genomewiderecombratefemale <- (sum(correlation$V5)/(sum(correlation$V2)/1000000))
genomewiderecombratefemale

rate.plot_m <-ggplot(correlation, aes(x=V2, y=male_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("Phyiscal Length (MB)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, data = correlation)+theme_light()
rate.plot_m

rate.plot_f <-ggplot(correlation, aes(x=V2, y=female_recomb_rates))+
  geom_point(color=  "darkred")+
  ggtitle("Female Map")+
  xlab("Phyiscal Length (MB)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, data = correlation)+theme_light()
rate.plot_f

ggarrange(corr.plot_m, corr.plot_f, rate.plot_m, rate.plot_f, labels = c("A)", "B)", "C)", "D)"), ncol = 2, nrow = 2 )
