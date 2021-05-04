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
# Not run
#daphnia.tree <- aboot(single_snp.genlight, tree = "upgma", distance = "bitwise.dist", sample = 100, showtree = F, cutoff = 50, quiet = T)

#plot the tree with ape:

cols <- c("blue","purple","red")

#plot.phylo(daphnia.tree, cex = 0.5, font = 0.3, adj = 0, tip.color = cols[pop(single_snp.genlight)])
#nodelabels(daphnia.tree$node.lable, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
#legend('topleft', legend = c("Daphnia Pulex","Hybrid", "Daphnia Pulicaria"), fill = cols, border = FALSE, bty = "n", cex = 2)
#axis(side = 1)
#title(xlab = "Bitwise Genetic distance (proportion of loci that are different)")


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

#Shape of Fst distribution:
fst <- pe_pi_singlesnp_fst$WEIR_AND_COCKERHAM_FST
mean(fst, na.rm = T)

Fst_distro <- ggplot(pe_pi_singlesnp_fst, aes(x = WEIR_AND_COCKERHAM_FST))+
  geom_histogram(color="darkblue", fill="lightblue")+
  xlab("Weir and Cockerham Fst")+
  ylab("Count")+
  theme_classic()
Fst_distro
################################################################################################################
# PCAdapt outlier analysis
setwd("~/Downloads")

# read in the vcf with just pulex and pulicaria:

pp_pcaadapt <- read.pcadapt("Pulex_Pulicaria_onesnp_pcadapt.vcf.gz", type = "vcf")

#Choose the best K of PCs for analysis:

scree <- pcadapt(input = pp_pcaadapt, K =20)
pcadapt_scree <- plot(scree, option = "screeplot")
pcadapt_scree <- pcadapt_scree + theme_light()

#Score plots:
pop.names <- c(rep("pulex", 82), rep("pulicaria", 41))

score_1_2 <-plot(scree, option = "scores", pop = pop.names)
score_1_2 <- score_1_2 +theme_light()

score_3_4 <- plot(scree, option = "scores", i =3, j= 4, pop=pop.names)
score_3_4 <- score_3_4 +theme_light()

#compute the test statistic:

Pca_test <- pcadapt(pp_pcaadapt, K=2)

summary(Pca_test)

#QQ plots:

pcadapt_qq <- plot(Pca_test, option = "qqplot")
pcadapt_qq <- pcadapt_qq +theme_light()

#The deviation away from the expected distribution confirms the presence of outliers:

pcadapt_values <- as.data.frame(Pca_test$pvalues)
pcadapt_values$`Pca_test$pvalues`
pcadapt_hist <- ggplot(pcadapt_values, aes(x=pcadapt_values$`Pca_test$pvalues`))+
  geom_histogram(color="darkblue", fill="lightblue")+
  theme_classic()+
  ylab("p-values")+
  xlab("Frequency")+
  ggtitle("PCAdapt P-values")
pcadapt_hist

# Plot of the test statistic:

test_stat <- plot(Pca_test, option = "stat.distribution")
test_stat <- test_stat + theme_light()

# Looks good: 

ggarrange(pcadapt_scree, score_1_2, score_3_4, pcadapt_qq, pcadapt_hist, test_stat, ncol = 2, nrow = 3, labels = c("A)", "B)", "C)", "D)", "E)", "F)"))


#Move on to outlier detection:

# Bonferroni Correction: (Conservative)

padj <- p.adjust(Pca_test$pvalues,method="bonferroni")
alpha <- 0.1
log(alpha)
outliers_BF_pcadapt <- which(padj < alpha)
length(outliers_BF_pcadapt)

write(outliers_BF_pcadapt, "outliers_pcaadapt.txt",ncolumns = 1, sep = "\t")
#Write out a Manhattan Plot:

pca.snps <- read.table("pcadapt.snps", sep = "\t", header = F)

SNP <- c(1:(nrow(pca.snps)))
pca.snps <- cbind.data.frame(pca.snps, log(Pca_test$pvalues), SNP)
head(pca.snps)

manhattan(pca.snps, chr = "V1", bp = "V2", p= "log(Pca_test$pvalues)", snp = "SNP", logp = F, ylim = c(0, -80),ylab = "log(P-value)", main = "Pulex-Pulicaria PCAdapt Outliers")

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
male_BP <- ggplot(analyze_recomb_fst, aes(x= outlier, y = male_rate))+
  geom_boxplot()+
  ylab("Male Recombination Rate (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  theme_light()
#female boxplot:
female_BP <- ggplot(analyze_recomb_fst, aes(x= outlier, y = female_rate))+
  geom_boxplot()+
  ylab("Female Recombination Rate (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  theme_light()
#outlier sites do not fall more in low recombining regions: 
#male rate

recomb_background[mapply(is.infinite, recomb_background)] <- NA
recomb_outliers[mapply(is.infinite, recomb_outliers)] <- NA

t.test(recomb_outliers$`1`, recomb_background$`1`, alternative = "two.sided")

t.test(recomb_outliers$`2`, recomb_background$`2`, alternative = "two.sided")


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
cds_BP <- ggplot(analyze_cds_fst, aes(x= outlier, y = CDS_content))+
  geom_boxplot()+
  ylab("Window CDS Content (100 KB Window Centered on SNP)")+
  xlab("Outlier Status")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  theme_light()

#Significant T test:

t.test(cds_outliers$cds_content, cds_background$cds_content, alternative = "greater")


ggarrange(male_BP, female_BP, cds_BP, Fst_distro, ncol = 2, nrow = 2, labels = c("A)", "B)", "C)", "D)"))

####################################################################################################################################################
# Sliding Window Fst pi and Dxy: 
setwd("~/Downloads")
windows <- read.csv("~/Downloads/Pulex_Pulicaria_allsnp_popgenwindows_1MB-10KB.csv")

#chromsome #1:
windows_01 <- subset(windows, windows$scaffold == "CHR01")

names(windows_01)

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_01$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_01, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 1")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 44000000)
 

#DXY Divergence: 
mean_Dxy <- mean(windows_01$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_01, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 1")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 44000000)
 

#Pulicaria pi

Puli_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 1")+
  theme_light()+ expand_limits(x = 44000000)

# PI Pulex
Pulex_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 1")+
  theme_light()+ expand_limits(x = 44000000)


#PI Hybrids

hybrid_pi <- ggplot(windows_01, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 1")+
  theme_light()+ expand_limits(x = 44000000)

#Plot Tajima's D on Chromosomes:

puli_D <- read.table("/home/weider/LinkageMap_Genomics/VCFs/tajimad_pulicaria.tsv", header = T, sep = "\t")
mean(puli_D$TajimaD)

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR01")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 1")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 44000000)

pulex_D <- read.table("/home/weider/LinkageMap_Genomics/VCFs/tajimad_pulex.tsv", header = T, sep = "\t")
mean(pulex_D$TajimaD)

Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR01")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 1")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 44000000)

rec_rate_1 <- ggplot(recombination_1, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 1 Recombination")+
  theme_light()+ expand_limits(x = 44000000)

rep_content_1 <- ggplot(CHR01_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 44000000)
#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_1, PA_PX_Fst, PA_PX_Dxy, rep_content, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


ggarrange(Puli_pi, Pulex_pi, heights = c(6,6), ncol = 1, nrow = 2)


#Chromosome 2: 
windows_02 <- subset(windows, windows$scaffold == "CHR02")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_02$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_02, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 2")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 13000000)


#DXY Divergence: 
mean_Dxy <- mean(windows_02$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_02, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 2")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 13000000)


#Pulicaria pi

Puli_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 2")+
  theme_light()+ expand_limits(x = 13000000)

# PI Pulex
Pulex_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 2")+
  theme_light()+ expand_limits(x = 13000000)


#PI Hybrids

hybrid_pi <- ggplot(windows_02, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 2")+
  theme_light()+ expand_limits(x = 13000000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR02")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 2")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 13000000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR02")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 2")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 13000000)

rec_rate_2 <- ggplot(recombination_2, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 2 Recombination")+
  theme_light()+ expand_limits(x = 13000000)
rep_content_2 <- ggplot(CHR02_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 13000000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_2, PA_PX_Fst, PA_PX_Dxy, rep_content_2, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 3: 
windows_03 <- subset(windows, windows$scaffold == "CHR03")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_03$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_03, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 3")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 13600000)


#DXY Divergence: 
mean_Dxy <- mean(windows_03$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_03, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 3")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 13600000)


#Pulicaria pi

Puli_pi <- ggplot(windows_03, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 3")+
  theme_light()+ expand_limits(x = 13600000)

# PI Pulex
Pulex_pi <- ggplot(windows_03, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 3")+
  theme_light()+ expand_limits(x = 13600000)


#PI Hybrids

hybrid_pi <- ggplot(windows_03, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 3")+
  theme_light()+ expand_limits(x = 13600000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR03")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 3")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 13600000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR03")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 3")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 13600000)

rec_rate_3 <- ggplot(recombination_3, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 3 Recombination")+
  theme_light()+ expand_limits(x = 13600000)
rep_content_3 <- ggplot(CHR03_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 13600000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_3, PA_PX_Fst, PA_PX_Dxy, rep_content_3, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)



#Chromosome 4: 
windows_04 <- subset(windows, windows$scaffold == "CHR04")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_04$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_04, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 4")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 22600000)


#DXY Divergence: 
mean_Dxy <- mean(windows_04$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_04, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 4")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 22600000)


#Pulicaria pi

Puli_pi <- ggplot(windows_04, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 4")+
  theme_light()+ expand_limits(x = 22600000)

# PI Pulex
Pulex_pi <- ggplot(windows_04, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 4")+
  theme_light()+ expand_limits(x = 22600000)


#PI Hybrids

hybrid_pi <- ggplot(windows_04, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 4")+
  theme_light()+ expand_limits(x = 22600000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR04")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 4")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 22600000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR04")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 4")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 22600000)

rec_rate_4 <- ggplot(recombination_4, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 4 Recombination")+
  theme_light()+ expand_limits(x = 22600000)
rep_content_4 <- ggplot(CHR04_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 22600000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_4, PA_PX_Fst, PA_PX_Dxy, rep_content_4, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 5: 
windows_05 <- subset(windows, windows$scaffold == "CHR05")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_05$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_05, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 5")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 7300000)


#DXY Divergence: 
mean_Dxy <- mean(windows_05$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_05, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 5")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 7300000)


#Pulicaria pi

Puli_pi <- ggplot(windows_05, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 5")+
  theme_light()+ expand_limits(x = 7300000)

# PI Pulex
Pulex_pi <- ggplot(windows_05, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 5")+
  theme_light()+ expand_limits(x = 7300000)


#PI Hybrids

hybrid_pi <- ggplot(windows_05, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 5")+
  theme_light()+ expand_limits(x = 7300000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR05")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 5")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 7300000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR05")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 5")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 7300000)

rec_rate_5 <- ggplot(recombination_5, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 5 Recombination")+
  theme_light()+ expand_limits(x = 7300000)
rep_content_5 <- ggplot(CHR05_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 7300000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_5, PA_PX_Fst, PA_PX_Dxy, rep_content_5, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 6: 
windows_06 <- subset(windows, windows$scaffold == "CHR06")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_06$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_06, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 6")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 24300000)


#DXY Divergence: 
mean_Dxy <- mean(windows_06$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_06, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 6")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 24300000)


#Pulicaria pi

Puli_pi <- ggplot(windows_06, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 6")+
  theme_light()+ expand_limits(x = 24300000)

# PI Pulex
Pulex_pi <- ggplot(windows_06, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 6")+
  theme_light()+ expand_limits(x = 24300000)


#PI Hybrids

hybrid_pi <- ggplot(windows_06, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 6")+
  theme_light()+ expand_limits(x = 24300000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR06")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 6")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 24300000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR06")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 6")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 24300000)

rec_rate_6 <- ggplot(recombination_6, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 6 Recombination")+
  theme_light()+ expand_limits(x = 24300000)
rep_content_6 <- ggplot(CHR06_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 24300000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_6, PA_PX_Fst, PA_PX_Dxy, rep_content_6, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 7: 
windows_07 <- subset(windows, windows$scaffold == "CHR07")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_07$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_07, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 7")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 18900000)


#DXY Divergence: 
mean_Dxy <- mean(windows_07$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_07, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 7")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 18900000)


#Pulicaria pi

Puli_pi <- ggplot(windows_07, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 7")+
  theme_light()+ expand_limits(x = 18900000)

# PI Pulex
Pulex_pi <- ggplot(windows_07, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 7")+
  theme_light()+ expand_limits(x = 18900000)


#PI Hybrids

hybrid_pi <- ggplot(windows_07, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 7")+
  theme_light()+ expand_limits(x = 18900000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR07")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 7")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 18900000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR07")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 7")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 18900000)

rec_rate_7 <- ggplot(recombination_7, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 7 Recombination")+
  theme_light()+ expand_limits(x = 18900000)
rep_content_7 <- ggplot(CHR07_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 18900000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_7, PA_PX_Fst, PA_PX_Dxy, rep_content_7, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 8: 
windows_08 <- subset(windows, windows$scaffold == "CHR08")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_08$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_08, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 8")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 16000000)


#DXY Divergence: 
mean_Dxy <- mean(windows_08$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_08, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 8")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 16000000)


#Pulicaria pi

Puli_pi <- ggplot(windows_08, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 8")+
  theme_light()+ expand_limits(x = 16000000)

# PI Pulex
Pulex_pi <- ggplot(windows_08, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 8")+
  theme_light()+ expand_limits(x = 16000000)


#PI Hybrids

hybrid_pi <- ggplot(windows_08, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 8")+
  theme_light()+ expand_limits(x = 16000000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR08")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 8")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 16000000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR08")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 8")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 16000000)

rec_rate_8 <- ggplot(recombination_8, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 8 Recombination")+
  theme_light()+ expand_limits(x = 16000000)
rep_content_8 <- ggplot(CHR08_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 16000000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_8, PA_PX_Fst, PA_PX_Dxy, rep_content_8, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)


#Chromosome 9: 
windows_09 <- subset(windows, windows$scaffold == "CHR09")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_09$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_09, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 9")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light() + expand_limits(x = 4500000)


#DXY Divergence: 
mean_Dxy <- mean(windows_09$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_09, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 9")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 4500000)


#Pulicaria pi

Puli_pi <- ggplot(windows_09, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 9")+
  theme_light()+ expand_limits(x = 4500000)

# PI Pulex
Pulex_pi <- ggplot(windows_09, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 9")+
  theme_light()+ expand_limits(x = 4500000)


#PI Hybrids

hybrid_pi <- ggplot(windows_09, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 9")+
  theme_light()+ expand_limits(x = 4500000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR09")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 9")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 4500000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR09")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 9")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 4500000)

rec_rate_9 <- ggplot(recombination_9, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 9 Recombination")+
  theme_light()+ expand_limits(x = 4500000)
rep_content_9 <- ggplot(CHR09_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 4500000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_9, PA_PX_Fst, PA_PX_Dxy, rep_content_9, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)



#Chromosome 10: 
windows_10<- subset(windows, windows$scaffold == "CHR10")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_10$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_10, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 10")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 7600000)


#DXY Divergence: 
mean_Dxy <- mean(windows_10$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_10, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 10")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 7600000)


#Pulicaria pi

Puli_pi <- ggplot(windows_10, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 10")+
  theme_light()+ expand_limits(x = 7600000)

# PI Pulex
Pulex_pi <- ggplot(windows_10, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 10")+
  theme_light()+ expand_limits(x = 7600000)


#PI Hybrids

hybrid_pi <- ggplot(windows_10, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 10")+
  theme_light()+ expand_limits(x = 7600000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR10")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 10")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 7600000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR10")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 10")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 7600000)

rec_rate_10 <- ggplot(recombination_10, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 10 Recombination")+
  theme_light()+ expand_limits(x = 7600000)
rep_content_10 <- ggplot(CHR10_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 7600000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_10, PA_PX_Fst, PA_PX_Dxy, rep_content_10, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)



#Chromosome 11: 
windows_11<- subset(windows, windows$scaffold == "CHR11")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_11$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_11, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 11")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light()+ expand_limits(x = 4700000)


#DXY Divergence: 
mean_Dxy <- mean(windows_11$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_11, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 11")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light()+ expand_limits(x = 4700000)


#Pulicaria pi

Puli_pi <- ggplot(windows_11, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 11")+
  theme_light()+ expand_limits(x = 4700000)

# PI Pulex
Pulex_pi <- ggplot(windows_11, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 11")+
  theme_light()+ expand_limits(x = 4700000)


#PI Hybrids

hybrid_pi <- ggplot(windows_11, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 11")+
  theme_light()+ expand_limits(x = 4700000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR11")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 11")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 4700000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR11")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 11")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 4700000)

rec_rate_11 <- ggplot(recombination_11, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 11 Recombination")+
  theme_light()+ expand_limits(x = 4700000)
rep_content_11 <- ggplot(CHR11_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 4700000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_11, PA_PX_Fst, PA_PX_Dxy, rep_content_11, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)



#Chromosome 12: 
windows_12<- subset(windows, windows$scaffold == "CHR12")

# Pulex-Pulicaria Fst
mean_fst <- mean(windows_12$Fst_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Fst <- ggplot(windows_12, aes(x=start, y= Fst_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "deeppink3")+
  ggtitle("Differentiation on LG 12")+
  ylab("Fst")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_fst)+
  theme_light() + expand_limits(x = 3500000)


#DXY Divergence: 
mean_Dxy <- mean(windows_12$dxy_Daphnia_pulex_Daphnia_pulicaria)

PA_PX_Dxy <- ggplot(windows_12, aes(x=start, y= dxy_Daphnia_pulex_Daphnia_pulicaria))+
  geom_line(color = "steelblue")+
  ylab("Dxy")+
  ggtitle("Divergence on LG 12")+
  xlab("Phyiscal Position")+
  geom_hline(yintercept = mean_Dxy)+
  theme_light() + expand_limits(x = 3500000)


#Pulicaria pi

Puli_pi <- ggplot(windows_12, aes(x= start, y = pi_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulicaria Nucleotide Diversity on LG 12")+
  theme_light()+ expand_limits(x = 3500000)

# PI Pulex
Pulex_pi <- ggplot(windows_12, aes(x= start, y = pi_Daphnia_pulex))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia pulex Nucleotide Diversity on LG 12")+
  theme_light()+ expand_limits(x = 3500000)


#PI Hybrids

hybrid_pi <- ggplot(windows_12, aes(x= start, y = pi_Daphnia_pulex_x_Daphnia_pulicaria))+
  geom_line()+
  xlab("Phyiscal Position")+
  ylab("Pi")+
  ggtitle("Daphnia hybrids Nucleotide Diversity on LG 12")+
  theme_light()+ expand_limits(x = 3500000)

#Plot Tajima's D on Chromosomes:

Puli_D_chr01 <- subset(puli_D, puli_D$CHROM == "CHR12")
mean_D <- mean(Puli_D_chr01$TajimaD)

TD <- ggplot(Puli_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulicaria Tajima's D on LG 12")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 3500000)



Pulex_D_chr01 <- subset(pulex_D, pulex_D$CHROM == "CHR12")

mean_D <- mean(Pulex_D_chr01$TajimaD)

TD_px <- ggplot(Pulex_D_chr01, aes(x=BIN_START, y=TajimaD))+
  geom_line()+
  xlab("Physical Position")+
  ylab("Tajima's D")+
  ggtitle("Daphnia pulex Tajima's D on LG 12")+
  geom_hline(yintercept = mean_D)+
  theme_light()+ expand_limits(x = 3500000)

rec_rate_12 <- ggplot(recombination_12, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  ylab("R (cM/MB)")+
  xlab("Physical Position")+
  scale_color_identity(guide = "legend", labels = c("Female Map", "Male Map"))+
  ggtitle("LG 12 Recombination")+
  theme_light()+ expand_limits(x = 3500000)
rep_content_12 <- ggplot(CHR12_repeat, aes(x=V1)) + 
  geom_line(aes(y = V3, color= "red"))+
  geom_line(aes(y = V4, color= "blue"))+
  geom_line(aes(y = V5, color= "darkgreen"))+
  geom_line(aes(y = V7, color = "black"))+
  xlab("Phyiscal Position")+ ylab("Window Content")+
  ggtitle("Sequence Content")+
  scale_color_identity(guide = "legend", labels = c("Coding Sequence", "DNA", "LINE", "LTR"))+
  theme_light()+ expand_limits(x = 3500000)

#Arrange the plots along the X axis:


plot <- ggarrange(rec_rate_12, PA_PX_Fst, PA_PX_Dxy, rep_content_12, heights = c(4,3,3,4), ncol = 1, nrow = 4, align = "v", axis = "2", legend = "bottom", common.legend = F)
plot$`1`

ggarrange(Puli_pi, TD, Pulex_pi, TD_px, heights = c(3,3,3,3), ncol = 1, nrow = 4)

