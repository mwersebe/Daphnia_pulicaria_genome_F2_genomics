# Windowed Repeat Content Calculations: 
# Custom R script by Matthew Wersebe
# Department of Biology, University of Oklahoma
##############################################################################
Step_size = 10000
Window_size = 1000000
#Full Gff file from repeat modeler
repeats <- read.table("SC_F1-1A_hap13B_scaffold_LOD26_purged.fa.out.gff", header = F, sep = "\t")
head(repeats)
#gff file with just LTR repeats. 
LTR <- read.table("LTR.gff", sep = "\t", header = F)
DNA <- read.table("DNA.gff", sep = "\t", header = F)
LINE <- read.table("LINE.gff", sep = "\t", header = F) 
UNKNOWN <- read.table("UNKNOWN.gff", sep = "\t", header = F)
coding <- read.table("cds.gff", sep = "\t", header = F)

repeat_content <- function(Data_total_repeat,Data_LTR, Data_DNA, Data_LINE, Data_UNKNOWN, Data_CODING, Window_size, Step_size, Chromosome){
  
  #Subset the total repeat gff file to the specificed chromosome
  
  chromosome_repeats = subset(Data_total_repeat, Data_total_repeat$V1 == Chromosome) 
  
  #Subset the the gff with the LTR and other data
  
  LTR_repeats <- subset(Data_LTR, Data_LTR$V1 == Chromosome)
  
  DNA_repeats <- subset(Data_DNA, Data_DNA$V1 == Chromosome)
  
  LINE_repeats <- subset(Data_LINE, Data_LINE$V1 == Chromosome)
  
  UNKNOWN_repeats <- subset(Data_UNKNOWN, Data_UNKNOWN$V1 == Chromosome)
  
  CODING_seq <- subset(Data_CODING, Data_CODING$V1 == Chromosome)
  
  #Determine the total length of the contig for estimating the number of windows
  Length <- max(chromosome_repeats$V5)
  Length
  Length_CDS <- max(CODING_seq$V5)
  Length_CDS
  #Number of windows Repeats
  
  number_windows <- (Length/Step_size)
  number_windows
  
  Windows = seq(from = 1, to = number_windows, by = 1)
  Windows
  
  #Initialize a matrix to hold the repeat content across the the chromosome.
  proportion = matrix(ncol = 7, nrow = length(Windows)) 
  #for loop fpr the gene content:
  
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(CODING_seq, CODING_seq$V4 > Step_size*(i - 1) & CODING_seq$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    cds_lengths <- (Window$V5 - Window$V4) #vector containing cds lengths
    
    total_cds <- sum(cds_lengths)
    
    prop_cds <- (total_cds/Window_size)
    
    #print the rate
    
    proportion[i,7] <- prop_cds}
  
  
  #for loop for the total repeat content
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(chromosome_repeats, chromosome_repeats$V4 > Step_size*(i - 1) & chromosome_repeats$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    repeat_lengths <- (Window$V5 - Window$V4) #vector containing repeat lengths
    
    total_repeat <- sum(repeat_lengths)
    
    prop <- (total_repeat/Window_size)
    
    #print the rate
    
    proportion[i,2] <- prop
    proportion[i,1] <- Step_size*(i)}
  
  #for loop for the LTR repeats
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(LTR_repeats, LTR_repeats$V4 > Step_size*(i - 1) & LTR_repeats$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    LTR_lengths <- (Window$V5 - Window$V4) #vector containing repeat lengths
    
    total_LTR <- sum(LTR_lengths)
    
    prop_LTR <- (total_LTR/Window_size)
    #print the rate
    proportion[i,3] <- prop_LTR}
  
  #DNA For loop:
  
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(DNA_repeats, DNA_repeats$V4 > Step_size*(i - 1) & DNA_repeats$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    DNA_lengths <- (Window$V5 - Window$V4) #vector containing repeat lengths
    
    total_DNA <- sum(DNA_lengths)
    
    prop_DNA <- (total_DNA/Window_size)
    #print the rate
    proportion[i,4] <- prop_DNA}
  
  #LINE For loop 
  
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(LINE_repeats, LINE_repeats$V4 > Step_size*(i - 1) & LINE_repeats$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    LINE_lengths <- (Window$V5 - Window$V4) #vector containing repeat lengths
    
    total_LINE <- sum(LINE_lengths)
    
    prop_LINE <- (total_LINE/Window_size)
    #print the rate
    proportion[i,5] <- prop_LINE}
  
  #Unknown for loop 
  for (i in 1:length(Windows)){
    
    #Determine the window
    Window <- subset(UNKNOWN_repeats, UNKNOWN_repeats$V4 > Step_size*(i - 1) & UNKNOWN_repeats$V5 < Window_size+(Step_size*(i - 1)))
    # calc repeat content  
    
    UNKNOWN_lengths <- (Window$V5 - Window$V4) #vector containing repeat lengths
    
    total_UNKNOWN <- sum(UNKNOWN_lengths)
    
    prop_UNKNOWN <- (total_UNKNOWN/Window_size)
    #print the rate
    proportion[i,6] <- prop_UNKNOWN}
  
  #Convert to data frame and return. 
  Repeats <- as.data.frame(proportion)
  
  return(Repeats)
}

#Example Use:

CHR01_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR01")

rep_content_1 <- ggplot(CHR01_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  



CHR02_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR02")


rep_content_2 <- ggplot(CHR02_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())


CHR03_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR03")


rep_content_3 <- ggplot(CHR03_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Centent")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


CHR04_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR04")


rep_content_4 <- ggplot(CHR04_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR05_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR05")


rep_content_5 <- ggplot(CHR05_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR06_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR06")


rep_content_6 <- ggplot(CHR06_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR07_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR07")


rep_content_7 <- ggplot(CHR07_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


CHR08_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR08")


rep_content_8 <- ggplot(CHR08_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR09_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR09")


rep_content_9 <- ggplot(CHR09_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR10_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR10")


rep_content_10 <- ggplot(CHR10_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR11_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR11")


rep_content_11 <- ggplot(CHR11_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CHR12_repeat <- repeat_content(Data_total_repeat = repeats, Data_LTR = LTR, Data_LINE = LINE, Data_DNA = DNA, Data_UNKNOWN = UNKNOWN, Data_CODING = coding, Window_size = Window_size, Step_size = Step_size, Chromosome = "CHR12")


rep_content_12 <- ggplot(CHR12_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

stack_1 <- ggarrange(rec_rate_1, rep_content_1, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_1 = stack_1$`1`

stack_2 <- ggarrange(rec_rate_2, rep_content_2, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_2 = stack_2$`1`

stack_3 <- ggarrange(rec_rate_3, rep_content_3, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_3 = stack_3$`1`

stack_4 <- ggarrange(rec_rate_4, rep_content_4, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_4 = stack_4$`1`

stack_5 = ggarrange(rec_rate_5, rep_content_5, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_5 = stack_5$`1`

stack_6 <- ggarrange(rec_rate_6, rep_content_6, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_6 = stack_6$`1`

stack_7 <- ggarrange(rec_rate_7, rep_content_7, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_7 = stack_7$`1`

stack_8 <- ggarrange(rec_rate_8, rep_content_8, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_8 = stack_8$`1`

stack_9 <- ggarrange(rec_rate_9, rep_content_9, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_9 = stack_9$`1`

stack_10 <- ggarrange(rec_rate_10, rep_content_10, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_10 = stack_10$`1`

stack_11 <- ggarrange(rec_rate_11, rep_content_11, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_11 = stack_11$`1`

stack_12 <- ggarrange(rec_rate_12, rep_content_12, ncol = 1, nrow = 2, align = "v", axis = "1", legend = "none", common.legend = F)
stack_12 = stack_12$`1`


ggarrange(stack_1, stack_2, stack_3, stack_4, stack_5, stack_6, stack_7, stack_8, stack_9, stack_10, stack_11, stack_12, nrow = 3, ncol = 4, legend = "bottom", common.legend = T)
