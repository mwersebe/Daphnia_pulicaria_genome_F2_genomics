# Windowed Repeat Content Calculations: 
# Custom R script by Matthew Wersebe
# Department of Biology, University of Oklahoma
# March 18,2021

#Repeat_content: Returns a dataframe with proportion of total repeat content, total LTR content, and position for plotting for a specific contig or chromosome.

setwd("/home/matt/Desktop/Linkage_plots/Gene_features")

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
head(CHR01_repeat)

rep_content <- ggplot(CHR01_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position (MB)")+ ylab("Window Coentent (Proportion)")+
  ggtitle("LG 1 Repeat Content Landscape")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))
rep_content

Phys_position <- CHR01_repeat$V1
total_repeat_content <- CHR01_repeat$V2
LTR_content <- CHR01_repeat$V3
DNA_content <- CHR01_repeat$V4
LINE_content <- CHR01_repeat$V5
UNKNOWN_content <- CHR01_repeat$V6
CODING_content<- CHR01_repeat$V7


male_window_recomb_rate <- recombination$V1
female_window_recomb_rate <- recombination$V2

sex_average_recomb_rate <- vector(mode = "numeric", length = length(recombination[,1]))
for (i in 1:length(sex_average_recomb_rate)){
  sex_average_recomb_rate[i] = ((recombination[i,1] + recombination[i,2])/2)
}

sex_average_recomb_rate

Rep_Recomb_correlation <-cbind.data.frame(Phys_position, total_repeat_content, LTR_content, DNA_content, LINE_content, CODING_content, male_window_recomb_rate, female_window_recomb_rate, sex_average_recomb_rate)
head(Rep_Recomb_correlation)

rate.plot_m_LTR <-ggplot(Rep_Recomb_correlation, aes(x=LTR_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("LTR Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_m_LTR

rate.plot_f_LTR <-ggplot(Rep_Recomb_correlation, aes(x=LTR_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("Female Map")+
  xlab("LTR Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_f_LTR

rate.plot_a_LTR <-ggplot(Rep_Recomb_correlation, aes(x=LTR_content, y=sex_average_recomb_rate))+
  geom_point(color= "black")+
  ggtitle("Sex Averaged Map")+
  xlab("LTR Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_a_LTR

rate.plot_m_LINE <-ggplot(Rep_Recomb_correlation, aes(x=LINE_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("LINE Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_m_LINE


rate.plot_f_LINE <-ggplot(Rep_Recomb_correlation, aes(x=LINE_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("female Map")+
  xlab("LINE Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_f_LINE

rate.plot_m_DNA <-ggplot(Rep_Recomb_correlation, aes(x=DNA_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("DNA Transposon Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_m_DNA


rate.plot_f_DNA <-ggplot(Rep_Recomb_correlation, aes(x=DNA_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("female Map")+
  xlab("DNA Transposon Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_f_DNA

rate.plot_m_CDS <-ggplot(Rep_Recomb_correlation, aes(x=CODING_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map")+
  xlab("Coding Sequence Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_m_CDS

rate.plot_f_CDS <-ggplot(Rep_Recomb_correlation, aes(x=CODING_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("female Map")+
  xlab("Coding Sequence Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_f_CDS

rate.plot_a_CDS <-ggplot(Rep_Recomb_correlation, aes(x=CODING_content, y=sex_average_recomb_rate))+
  geom_point(color= "black")+
  ggtitle("Sex Averaged Map")+
  xlab("Coding Sequence Content (Proportion)")+ ylab("Sex Averaged Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_a_CDS

rate.plot_a_TRep<-ggplot(Rep_Recomb_correlation, aes(x=total_repeat_content, y=female_window_recomb_rate))+
  geom_point(color= "black")+
  ggtitle("Sex Averaged Map")+
  xlab("Total Repeat Content (Proportion)")+ ylab("Sex Averaged Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
rate.plot_a_TRep