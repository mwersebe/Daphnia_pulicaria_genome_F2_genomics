## South Center Genome
## Populations Genomics Using the SRA Daphnia pulex and Daphnia pulicaria 
## ddRADseq.
## Matthew Wersebe (University of Oklahoma) 4-5-2021
#############################################################################
#set working directotry:

setwd("/home/weider/LinkageMap_Genomics/")

# Load Required Libraries

library(vcfR)
library(adegenet)
library(ape)
library(poppr)
library(ggplot2)
library(tidyverse)
library(qqman)
library(pcadapt)
library(qvalue)
###############################################################################

#read in the single snp vcf file to perform DAPC

single_snp.vcf <- read.vcfR("/home/weider/LinkageMap_Genomics/VCFs/Pulex_Pulicaria_onesnp.vcf.gz")

# 32949 variants in total 

# Read in the popmap to assign the population to the correct individuals

pop.map <- read.table("popmap.txt", sep = "\t", header = F)
head(pop.map)

#convert the vcfR object to a genlight for analysis:

single_snp.genlight <- vcfR2genlight(single_snp.vcf)

#add important meta data

ploidy(single_snp.genlight) <- 2
pop(single_snp.genlight) <- pop.map$V2

# no issues move on to genetic distances and population structure:
# Genetic Distance Tree:

daphnia.tree <- aboot(single_snp.genlight, tree = "upgma", distance = "bitwise.dist", sample = 100, showtree = F, cutoff = 50, quiet = T)

#plot the tree with ape:

cols <- c("aquamarine","orange","red")

plot.phylo(daphnia.tree, cex = 0.5, font = 0.3, adj = 0, tip.color = cols[pop(single_snp.genlight)])
nodelabels(daphnia.tree$node.lable, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
legend('topleft', legend = c("Daphnia Pulex","Hybrid", "Daphnia Pulicaria"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Bitwise Genetic distance (proportion of loci that are different)")


###############################################################################

# PCA:

daphnia.pca <- glPca(single_snp.genlight, nf = 3)

#PCA Plots:
#Eigenvalues
barplot(100*daphnia.pca$eig/sum(daphnia.pca$eig), col = heat.colors(180), main = "PCA Eigenvalues")
title(ylab="Percent Variance\nexplained", line =2)
title(xlab="Eignevalues", line = 1)

#Scatter 
pca.scores <- as.data.frame(daphnia.pca$scores)
pca.scores$pop <- pop(single_snp.genlight)

daphnia.plot <- ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop))+
  geom_point(size = 1)+
  stat_ellipse(level = 0.95, size = 0.5)+
  scale_color_manual(values = cols)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()
daphnia.plot

###############################################################################
###############################################################################


#Adaptive Genomics:
#Re set working directory:
setwd("/home/matt/Linkag_Genomics")

pe_pi_singlesnp_fst <- read.table("pulex_pulicaria_onesnp.weir.fst", sep = "\t", header = T)
head(pe_pi_singlesnp_fst)
quantile(pe_pi_singlesnp_fst$WEIR_AND_COCKERHAM_FST, c(0.975, 0.995), na.rm = T)

#determine the outlier space:

threshold.fst <- quantile(pe_pi_singlesnp_fst$WEIR_AND_COCKERHAM_FST, c(0.975), na.rm = T)
threshold.fst

#figure out whichones are outliers:
pe_pi_singlesnp_fst <- pe_pi_singlesnp_fst %>% mutate(outlier = ifelse(WEIR_AND_COCKERHAM_FST > threshold.fst, "outlier", "background"))

#tells you how many outliers

pe_pi_singlesnp_fst %>% group_by(outlier) %>% tally()

outliers <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$outlier == "outlier")
head(outliers)

# Manhattan plot:

pe_pi_subset <- pe_pi_singlesnp_fst[complete.cases(pe_pi_singlesnp_fst),]

PE_PI_SNPs <- c(1:(nrow(pe_pi_subset)))

SNP_df <- data.frame(PE_PI_SNPs, pe_pi_subset)

manhattan(SNP_df, chr = "CHROM", bp = "POS", p= "WEIR_AND_COCKERHAM_FST", snp = "PE_PI_SNPs", logp = F, ylab = "Weir and Cockerham Fst", ylim = c(0, 1.1), main = "Pulex-Pulicaria Per Site Fst")
abline(h=threshold.fst, col = "red")

################################################################################################################
# PCAdapt outlier analysis


# read in the vcf with just pulex and pulicaria:

pp_pcaadapt <- read.pcadapt("Pulex_Pulicaria_onesnp_pcadapt.vcf.gz", type = "vcf")

#Choose the best K of PCs for analysis:

scree <- pcadapt(input = pp_pcaadapt, K =20)
plot(scree, option = "screeplot")

  #Score plots:
pop.names <- c(rep("pulex", 82), rep("pulicaria", 41))
plot(scree, option = "scores", pop = pop.names)
plot(scree, option = "scores", i =3, j= 4, pop=pop.names)

#compute the test statistic:

Pca_test <- pcadapt(pp_pcaadapt, K=2)

summary(Pca_test)

#Manhattan Plots:

plot(Pca_test, option = "manhattan")

#QQ plots:

plot(Pca_test, option = "qqplot")

#The deviation away from the expected distribution confirms the presence of outliers:

hist(Pca_test$pvalues, xlab = "p-values", main = "Histogram of P-values", breaks = 50, col = "grey")

# Plot of the test statistic:

plot(Pca_test, option = "stat.distribution")

# Looks good: 

#Move on to outlier detection:

# Bonferroni Correction: (Conservative)

padj <- p.adjust(Pca_test$pvalues,method="bonferroni")
alpha <- 0.1
outliers_BF_pcadapt <- which(padj < alpha)
length(outliers_BF_pcadapt)


# To DO: 

# Graph the windowed (100 KB 50 % overlap) Fst, Pi and Dxy estimates from the from Genomic General Ouput (Results On OSCER).

# Finalize the graphing: Marey Maps, Recomb Landscapes, Differentiation v. Recomb (Boxplots), pop structure, Graph the linkage map (marker density).

# Perform/graph correlations between recomb rate, repeat content, diversity + divergence and differentiation. 

# Collate Results: 

# Done.

#####################################################################################
# Fst outliers and Recombination Rate Variation:

outliers <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$outlier == "outlier")
head(outliers)

#Function Outliers_recomb_rate produceds a dataframe with the outlier and the recomb rate for further analysis. 
Outlier_recomb_rate = function(Outliers, Marey_Map, Chromosome){

#get all the outliers on the chromosome:
chromosome_outliers = subset(Outliers, Outliers$CHROM == Chromosome)

#Initalize a vector with a length the same as chromosome_outliers:

recom_rate <- matrix(ncol = 2, nrow = length(chromosome_outliers$POS))

out_pos <- chromosome_outliers$POS

Size_scalar = (1000000/1000000)
#for loop to perfrom the the calcualtions:

for (i in 1:length(out_pos)){
  
  #Determine the start/end:
  
  start_window <- out_pos[i] - 500000
  end_window <- out_pos[i] + 500000
  
  # Extract the the Data from the marey map 
  
  window = subset(Marey_map, Marey_map$Physical_position > start_window & Marey_map$Physical_position < end_window)
  head(window)
  # Determine the rate (i.e., slope) by fitting a linear model:
  
  male = (((max(window$cM_male))-(min(window$cM_male)))/(Size_scalar))
  female = (((max(window$cM_female))-(min(window$cM_female)))/(Size_scalar))
  
  
  # Put the data in the matrix recomb_rate 
  
  recom_rate[i,1] <- male
  recom_rate[i,2] <- female
}

Outlier_recomb_rate <- cbind.data.frame(chromosome_outliers, recom_rate)

return(Outlier_recomb_rate)}

#Outliers from fst outlier test:
Outliers = outliers

#chromosome 1

Marey_map_01 = read.table("/home/matt/Downloads/SC_map_Chrom01-marey.txt", sep = "\t", header = T)

chr01 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_01, Chromosome = 1)

#chromosome 2 
Marey_map_02 <- read.table("/home/matt/Downloads/SC_map_Chrom02-marey.txt", sep = "\t", header = T)

chr02 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_02, Chromosome = 2)
head(chr02)

#chromosome 3 
Marey_map_03 <- read.table("/home/matt/Downloads/SC_map_Chrom03-marey.txt", sep = "\t", header = T)

chr03 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_03, Chromosome = 3)


#chromosome 4
Marey_map_04 <- read.table("/home/matt/Downloads/SC_map_Chrom04-marey.txt", sep = "\t", header = T)

chr04 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_04, Chromosome = 4)


#chromosome 5 
Marey_map_05 <- read.table("/home/matt/Downloads/SC_map_Chrom05-marey.txt", sep = "\t", header = T)

chr05 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_05, Chromosome = 5)


#chromosome 6 
Marey_map_06 <- read.table("/home/matt/Downloads/SC_map_Chrom06-marey.txt", sep = "\t", header = T)

chr06 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_06, Chromosome = 6)


#chromosome 7 
Marey_map_07 <- read.table("/home/matt/Downloads/SC_map_Chrom07-marey.txt", sep = "\t", header = T)

chr07 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_07, Chromosome = 7)

#chromosome 8 
Marey_map_08 <- read.table("/home/matt/Downloads/SC_map_Chrom08-marey.txt", sep = "\t", header = T)

chr08 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_08, Chromosome = 8)

#chromosome 9 
Marey_map_09 <- read.table("/home/matt/Downloads/SC_map_Chrom09-marey.txt", sep = "\t", header = T)

chr09 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_09, Chromosome = 9)

#chromosome 10
Marey_map_10 <- read.table("/home/matt/Downloads/SC_map_Chrom10-marey.txt", sep = "\t", header = T)

chr10 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_10, Chromosome = 10)
chr10
#chromosome 11 
Marey_map_11 <- read.table("/home/matt/Downloads/SC_map_Chrom11-marey.txt", sep = "\t", header = T)

chr11 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_11, Chromosome = 11)
chr11
#chromosome 12 
Marey_map_12 <- read.table("/home/matt/Downloads/SC_map_Chrom12-marey.txt", sep = "\t", header = T)

chr12 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_12, Chromosome = 12)

#make a data frame: 
recomb_outliers <- rbind.data.frame(chr01, chr02, chr03, chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12)
tail(recomb_outliers)
length(recomb_outliers[,1])



# Subsample a group of 762 "background" Loci and measure the recomb rate:

background <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$outlier == "background")
#set.seed(273)
#Background <- sample_n(background, 762)

Background = background
#chr's 

chr01 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_01, Chromosome = 1)

chr02 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_02, Chromosome = 2)

chr03 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_02, Chromosome = 2)

chr04 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_04, Chromosome = 4)

chr05 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_05, Chromosome = 5)

chr06 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_06, Chromosome = 6)

chr07 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_07, Chromosome = 7)

chr08 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_08, Chromosome = 8)

chr09 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_09, Chromosome = 9)

chr10 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_10, Chromosome = 10)

chr11 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_11, Chromosome = 11)

chr12 <- Outlier_recomb_rate(Outliers = Background, Marey_Map = Marey_map_12, Chromosome = 12)
#make a data frame from results:
recomb_background <- rbind.data.frame(chr01, chr02, chr03, chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12)


#combine for analysis and name columns:

analyze_recomb_fst <- rbind.data.frame(recomb_outliers, recomb_background)

analyze_recomb_fst <- analyze_recomb_fst %>% mutate(outlier = ifelse(WEIR_AND_COCKERHAM_FST > threshold.fst, "outlier", "background"))

head(analyze_recomb_fst)

tail(analyze_recomb_fst)

names(analyze_recomb_fst)[5] <- "male_rate"
names(analyze_recomb_fst)[6] <- "female_rate"

#male boxplot
ggplot(analyze_recomb_fst, aes(x= outlier, y = male_rate))+
  geom_boxplot()+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)
#female boxplot:
ggplot(analyze_recomb_fst, aes(x= outlier, y = female_rate))+
  geom_boxplot()+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)
#outlier sites do not fall more in low recombining regions: 
########################################################################################################################################
# Window Estimates: 

windows <- read.csv("Pulex_Pulicaria_allsnp_popgenwindows.csv", header =T)
head(windows)
Window_recomb_rate = function(Windows, Marey_Map, Chromosome, Window_size){
#constant
Size_scalar = (Window_size/1000000)

chromosome_windows = subset(Windows, Windows$scaffold == Chromosome)

window_start = chromosome_windows$start
window_end = chromosome_windows$end

recom_rate <- matrix(ncol = 2, nrow = length(window_start))

for (i in 1:length(window_start)){
  
  window = subset(Marey_map_01, Marey_map_01$Physical_position > window_start[i] & Marey_map_01$Physical_position < window_end[i])
  
  male = (((max(window$cM_male))-(min(window$cM_male)))/(Size_scalar))
  female = (((max(window$cM_female))-(min(window$cM_female)))/(Size_scalar))
  
  recom_rate[i,1] <- male
  recom_rate[i,2] <- female
}

return(as.data.frame(recom_rate))
  
}

Chr01 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_01, Chromosome = "CHR01", Window_size = 100000)

Chr02 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_02, Chromosome = "CHR02", Window_size = 100000)

Chr03 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_03, Chromosome = "CHR03", Window_size = 100000)

Chr04 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_04, Chromosome = "CHR04", Window_size = 100000)

Chr05 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_05, Chromosome = "CHR05", Window_size = 100000)

Chr06 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_06, Chromosome = "CHR06", Window_size = 100000)

Chr07 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_07, Chromosome = "CHR07", Window_size = 100000)

Chr08 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_08, Chromosome = "CHR08", Window_size = 100000)

Chr09 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_09, Chromosome = "CHR09", Window_size = 100000)

Chr10 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_10, Chromosome = "CHR10", Window_size = 100000)

Chr11 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_11, Chromosome = "CHR11", Window_size = 100000)

Chr12 <- Window_recomb_rate(Windows = windows, Marey_Map = Marey_map_12, Chromosome = "CHR12", Window_size = 100000)


recomb_window <- rbind.data.frame(Chr01, Chr02, Chr03, Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12)


windows = cbind.data.frame(windows, recomb_window)

names(windows)[15] <- "male_rate"
names(windows)[16] <- "female_rate"

head(windows)

windows %>% 
  mutate_if(is.numeric, list(~na_if(., "-Inf")))
windows %>% 
  mutate_if(is.numeric, list(~na_if(., "Inf")))

ggplot(windows, aes(x = male_rate,y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_point()+
  stat_smooth(method = "lm")
ggplot(windows, aes(x = female_rate,y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_point()+
  stat_smooth(method = "lm")

CHR01 <- subset(windows, windows$scaffold == "CHR01")
head(CHR01)
ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = pi_Daphnia_pulicaria, color = "black"))
ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = pi_Daphnia_pulex, color = "black"))

ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = Fst_Daphnia_pulex_Daphnia_pulicaria, color = "black"))
ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = female_rate, color = "dark red"))

ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = pi_Daphnia_pulex_x_Daphnia_pulicaria, color = "dark red"))

ggplot(CHR01, aes(x = start))+
  geom_line(aes(y = dxy_Daphnia_pulex_Daphnia_pulicaria, color = "dark red"))

CHR02 <- subset(windows, windows$scaffold == "CHR02")

ggplot(CHR02, aes(x = start))+
  geom_line(aes(y = pi_Daphnia_pulicaria, color = "black"))
ggplot(CHR02, aes(x = start))+
  geom_line(aes(y = Fst_Daphnia_pulex_Daphnia_pulicaria, color = "black"))
ggplot(CHR02, aes(x = start))+
  geom_line(aes(y = female_rate, color = "dark red"))


CHR03 <- subset(windows, windows$scaffold == "CHR03")

ggplot(CHR03, aes(x = start))+
  geom_line(aes(y = pi_Daphnia_pulicaria, color = "black"))
ggplot(CHR03, aes(x = start))+
  geom_line(aes(y = Fst_Daphnia_pulex_Daphnia_pulicaria, color = "black"))
ggplot(CHR03, aes(x = start))+
  geom_line(aes(y = female_rate, color = "dark red"))
ggplot(CHR03, aes(x = start))+
  geom_line(aes(y = dxy_Daphnia_pulex_Daphnia_pulicaria, color = "dark red"))



