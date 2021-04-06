

CHR01_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR01")

rep_content <- ggplot(CHR01_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 1 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


CHR02_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR02")

rep_content <- ggplot(CHR02_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 2 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR03_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR03")

rep_content <- ggplot(CHR03_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 3 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR04_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR04")

rep_content <- ggplot(CHR04_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 4 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR05_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR05")

rep_content <- ggplot(CHR05_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 5 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR06_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR06")

rep_content <- ggplot(CHR06_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 6 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


CHR07_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR07")

rep_content <- ggplot(CHR07_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 7 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR08_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR08")

rep_content <- ggplot(CHR08_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 8 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


CHR09_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR09")

rep_content <- ggplot(CHR09_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 9 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


CHR10_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR10")

rep_content <- ggplot(CHR10_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 10 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


CHR11_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR11")

rep_content <- ggplot(CHR11_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 11 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

CHR12_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR12")

rep_content <- ggplot(CHR12_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 12 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content


######################################################################################
#Repeat Landscape:

library(reshape)
library(ggplot2)
library(viridis)
library(tidyverse)
library(gridExtra)

landscape <- read.csv("Repeat_landscape.csv", sep = " ",header = T)
head(landscape)

genome_size = 185212202
kd_melt = melt(landscape,id="Div")
kd_melt$norm = kd_melt$value/genome_size * 100

TELandscape <- ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black") +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  labs(fill = "Repeat Type") +
  scale_fill_viridis(discrete = T)+
  coord_cartesian(xlim = c(0, 50))
TELandscape


