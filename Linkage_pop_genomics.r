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
library(ggpubr)
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

cols <- c("blue","purple","red")

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
  geom_point(size = 2)+
  stat_ellipse(level = 0.95, size = 1)+
  scale_color_manual(values = cols)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()
daphnia.plot

###############################################################################
###############################################################################
# Adatpive Genomics: 
#Re set working directory:
#setwd("/home/matt/Linkag_Genomics")

pe_pi_singlesnp_fst <- read.table("/home/weider/LinkageMap_Genomics/SlidingWindows/pulex_pulicaria_onesnp.weir.fst", sep = "\t", header = T)
head(pe_pi_singlesnp_fst)
tail(pe_pi_singlesnp_fst)
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

# By Chromosome Manhattans:

CHR01 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 01)
ggplot(CHR01, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")
  
CHR02 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 02)
ggplot(CHR02, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst") 

CHR03 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 03)
ggplot(CHR03, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR04 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 04)
ggplot(CHR04, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR05 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 05)
ggplot(CHR05, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR06 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 06)
ggplot(CHR06, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR07 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 07)
ggplot(CHR07, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR08 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 08)
ggplot(CHR08, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR09 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 09)
ggplot(CHR09, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR10 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 10)
ggplot(CHR10, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR11 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 11)
ggplot(CHR11, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

CHR12 <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$CHROM == 12)
ggplot(CHR12, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = outlier))+
  geom_point()+
  xlab("Phyiscal Position along Chromosome")+
  ylab("Weir and Cockerham Fst")

#Shape of Fst distribution:

ggplot(pe_pi_singlesnp_fst, aes(x = WEIR_AND_COCKERHAM_FST))+
  geom_histogram(color="darkblue", fill="lightblue")+
  xlab("Weir and Cockerham Fst")+
  ylab("Count")+
  theme_classic()

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

outliers_BF_pcadapt

write(outliers_BF_pcadapt, "outliers_pcaadapt.txt",ncolumns = 1, sep = "\t")
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
  
  Size_scalar = (100000/1000000)
  #for loop to perfrom the the calcualtions:
  
  for (i in 1:length(out_pos)){
    
    #Determine the start/end:
    
    start_window <- out_pos[i] - 50000
    end_window <- out_pos[i] + 50000
    
    # Extract the the Data from the marey map 
    
    window = subset(Marey_Map, Marey_Map$Physical_position > start_window & Marey_Map$Physical_position < end_window)
    
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

Marey_map_01 = read.table("/home/weider/Downloads/SC_map_Chrom01-marey.txt", sep = "\t", header = T)

chr01 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_01, Chromosome = 1)

#chromosome 2 
Marey_map_02 <- read.table("/home/weider/Downloads/SC_map_Chrom02-marey.txt", sep = "\t", header = T)

chr02 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_02, Chromosome = 2)
head(chr02)

#chromosome 3 
Marey_map_03 <- read.table("/home/weider/Downloads/SC_map_Chrom03-marey.txt", sep = "\t", header = T)

chr03 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_03, Chromosome = 3)


#chromosome 4
Marey_map_04 <- read.table("/home/weider/Downloads/SC_map_Chrom04-marey.txt", sep = "\t", header = T)

chr04 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_04, Chromosome = 4)


#chromosome 5 
Marey_map_05 <- read.table("/home/weider/Downloads/SC_map_Chrom05-marey.txt", sep = "\t", header = T)

chr05 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_05, Chromosome = 5)


#chromosome 6 
Marey_map_06 <- read.table("/home/weider/Downloads/SC_map_Chrom06-marey.txt", sep = "\t", header = T)

chr06 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_06, Chromosome = 6)


#chromosome 7 
Marey_map_07 <- read.table("/home/weider/Downloads/SC_map_Chrom07-marey.txt", sep = "\t", header = T)

chr07 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_07, Chromosome = 7)

#chromosome 8 
Marey_map_08 <- read.table("/home/weider/Downloads/SC_map_Chrom08-marey.txt", sep = "\t", header = T)

chr08 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_08, Chromosome = 8)

#chromosome 9 
Marey_map_09 <- read.table("/home/weider/Downloads/SC_map_Chrom09-marey.txt", sep = "\t", header = T)

chr09 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_09, Chromosome = 9)

#chromosome 10
Marey_map_10 <- read.table("/home/weider/Downloads/SC_map_Chrom10-marey.txt", sep = "\t", header = T)

chr10 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_10, Chromosome = 10)
chr10
#chromosome 11 
Marey_map_11 <- read.table("/home/weider/Downloads/SC_map_Chrom11-marey.txt", sep = "\t", header = T)

chr11 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_11, Chromosome = 11)

#chromosome 12 
Marey_map_12 <- read.table("/home/weider/Downloads/SC_map_Chrom12-marey.txt", sep = "\t", header = T)

chr12 <- Outlier_recomb_rate(Outliers = Outliers, Marey_Map = Marey_map_12, Chromosome = 12)

#make a data frame: 
recomb_outliers <- rbind.data.frame(chr01, chr02, chr03, chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12)
tail(recomb_outliers)
length(recomb_outliers[,1])



# Subsample a group of 762 "background" Loci and measure the recomb rate:

background <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$outlier == "background")
set.seed(273)
Background <- sample_n(background, 762)

#Background = background
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
  ylab("Male Recombination Rate (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)
#female boxplot:
ggplot(analyze_recomb_fst, aes(x= outlier, y = female_rate))+
  geom_boxplot()+
  ylab("Female Recombination Rate (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)
#outlier sites do not fall more in low recombining regions: 


#################################################################################################################
# Coding Seq Density and Outlier status: 

Outlier_gene_content <- function(Outliers, CDS_gff, Chromosome, Chromosome_gff, Window_size){
  
  chromosome_outliers = subset(Outliers, Outliers$CHROM == Chromosome)
  
  chromosome_cds = subset(CDS_gff, CDS_gff$V1 == Chromosome_gff)
  
  cds_content <- matrix(ncol = 1, nrow = length(chromosome_outliers$POS))
  
  out_pos <- chromosome_outliers$POS
  
  Window_size <- Window_size
  
  for (i in 1:length(out_pos)){
    
    #Determine the start/end:
    
    start_window <- out_pos[i] - 50000
    end_window <- out_pos[i] + 50000
    
    Window <- subset(chromosome_cds, chromosome_cds$V4 > start_window & chromosome_cds$V5 < end_window)
    
    cds_lengths <- (Window$V5 - Window$V4)
    
    total_cds <- sum(cds_lengths)
    
    prop_cds <- (total_cds/Window_size)
    
    cds_content[i,1] <- prop_cds

  }
  
  Outlier_cds_content <- cbind.data.frame(chromosome_outliers, cds_content)
  
  return(Outlier_cds_content)
  }

CDS_gff <- read.table("~/Downloads/cds.gff", sep = "\t", header = F)

CHR01 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 01, Chromosome_gff = "CHR01", Window_size = 100000)

CHR02 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 02, Chromosome_gff = "CHR02", Window_size = 100000)

CHR03 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 03, Chromosome_gff = "CHR03", Window_size = 100000)

CHR04 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 04, Chromosome_gff = "CHR04", Window_size = 100000)

CHR05 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 05, Chromosome_gff = "CHR05", Window_size = 100000)

CHR06 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 06, Chromosome_gff = "CHR06", Window_size = 100000)

CHR07 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 07, Chromosome_gff = "CHR07", Window_size = 100000)

CHR08 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 08, Chromosome_gff = "CHR08", Window_size = 100000)

CHR09 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 09, Chromosome_gff = "CHR09", Window_size = 100000)

CHR10 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 10, Chromosome_gff = "CHR10", Window_size = 100000)

CHR11 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 11, Chromosome_gff = "CHR11", Window_size = 100000)

CHR12 <- Outlier_gene_content(Outlier= Outliers, CDS_gff = CDS_gff, Chromosome = 12, Chromosome_gff = "CHR12", Window_size = 100000)

cds_outliers <- rbind.data.frame(CHR01, CHR02, CHR03, CHR04, CHR05, CHR06, CHR07, CHR08, CHR09, CHR10, CHR11, CHR12)

CHR01 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 01, Chromosome_gff = "CHR01", Window_size = 100000)

CHR02 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 02, Chromosome_gff = "CHR02", Window_size = 100000)

CHR03 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 03, Chromosome_gff = "CHR03", Window_size = 100000)

CHR04 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 04, Chromosome_gff = "CHR04", Window_size = 100000)

CHR05 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 05, Chromosome_gff = "CHR05", Window_size = 100000)

CHR06 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 06, Chromosome_gff = "CHR06", Window_size = 100000)

CHR07 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 07, Chromosome_gff = "CHR07", Window_size = 100000)

CHR08 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 08, Chromosome_gff = "CHR08", Window_size = 100000)

CHR09 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 09, Chromosome_gff = "CHR09", Window_size = 100000)

CHR10 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 10, Chromosome_gff = "CHR10", Window_size = 100000)

CHR11 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 11, Chromosome_gff = "CHR11", Window_size = 100000)

CHR12 <- Outlier_gene_content(Outlier= Background, CDS_gff = CDS_gff, Chromosome = 12, Chromosome_gff = "CHR12", Window_size = 100000)

cds_background <- rbind.data.frame(CHR01, CHR02, CHR03, CHR04, CHR05, CHR06, CHR07, CHR08, CHR09, CHR10, CHR11, CHR12)


#combine for analysis and name columns:

analyze_cds_fst <- rbind.data.frame(cds_outliers, cds_background)


analyze_cds_fst <- analyze_cds_fst %>% mutate(outlier = ifelse(WEIR_AND_COCKERHAM_FST > threshold.fst, "outlier", "background"))


names(analyze_cds_fst)[5] <- "CDS_content"


#Cds boxplot
ggplot(analyze_cds_fst, aes(x= outlier, y = CDS_content))+
  geom_boxplot()+
  ylab("Proportion of Window as Coding Sequence (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)

#Significant T test:

t.test(cds_outliers$cds_content, cds_background$cds_content, alternative = "greater")


####################################################################################################################################################
# Sliding Window Fst pi and Dxy: 

windows <- read.csv("~/Downloads/Pulex_Pulicaria_allsnp_popgenwindows_1MB-10KB.csv")


windows_01 <- subset(windows, windows$scaffold == "CHR01")

names(windows_01)
#Pulicaria pi
Puli_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 1")+
  theme_light()
  
# Pulex-Pulicaria Fst
PA_PX_Fst <- ggplot(windows_01, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "pink")+
  ggtitle("Differentiation on LG 1")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  theme_light()
  
#DXY Divergence: 
PA_PX_Dxy <- ggplot(windows_01, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "light blue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 1")+
  xlab("Phyiscal Position")+
  theme_light()

# PI Pulex
Pulex_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 1")+
  theme_light()


#PI Hybrids

hybrid_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 1")+
  theme_light()

#Fst Pulicaria and Hybrids: 
ggplot(windows_01, aes(x=start, y=Fst_Daphnia_pulicaria_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position on Chromosome")+
  ylab("Fst Daphnia pulicaria Vs Daphnia Hybrids")

#Arrange the plots along the X axis:
library(cowplot)

plot <- ggarrange(rec_rate_1, PA_PX_Fst, PA_PX_Dxy, rep_content, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, Pulex_pi, hybrid_pi, heights = c(3,3,3), ncol = 1, nrow = 3)




windows_02 <- subset(windows, windows$scaffold == "CHR02")

#Pulicaria pi
Puli_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 2")+
  theme_light()
  

# Pulex- Puulicaria Fst
PA_PX_Fst <- ggplot(windows_02, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "pink")+
  ggtitle("Differentiation on LG 2")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  theme_light()

#DXY Divergence: 
PA_PX_Dxy <- ggplot(windows_02, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
 geom_line(color = "light blue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 2")+
  xlab("Phyiscal Position")+
  theme_light()


# PI Pulex
Pulex_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 2")+
  theme_light()

#PI Hybrids

hybrid_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrid Nucleotide Diversity on LG 2")+
  theme_light()

#Fst Pulicaria and Hybrids: 
ggplot(windows_02, aes(x=start, y=Fst_Daphnia_pulicaria_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position on Chromosome")+
  ylab("Fst Daphnia pulicaria Vs Daphnia Hybrids")

plot <- ggarrange(rec_rate_2, PA_PX_Fst, PA_PX_Dxy, rep_content, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, Pulex_pi, hybrid_pi, heights = c(3,3,3), ncol = 1, nrow = 3)
