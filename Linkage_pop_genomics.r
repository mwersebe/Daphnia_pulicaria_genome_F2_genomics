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
library(qqman)

#Adaptive Genomics:
# Site-wise Fst outliers (calculated by VCFtools)
# Pulex to Pulicaria
fst_pe_pi <- read.table("/home/weider/LinkageMap_Genomics/SlidingWindows/pulex_pulicaria.weir.fst", header = T)
head(fst_pe_pi)
#Using qqman 
pe_pi_subset <- fst_pe_pi[complete.cases(fst_pe_pi),]

PE_PI_SNPs <- c(1:(nrow(pe_pi_subset)))

SNP_df <- data.frame(PE_PI_SNPs, pe_pi_subset)

manhattan(SNP_df, chr = "CHROM", bp = "POS", p= "WEIR_AND_COCKERHAM_FST", snp = "PE_PI_SNPs", logp = F, ylab = "Weir and Cockerham Fst", ylim = c(0, 1.1), main = "Pulex-Pulicaria Per Site Fst")

#Pulex to Hybrid

fst_pe_hy <- read.table("/home/weider/LinkageMap_Genomics/SlidingWindows/pulex_hybrid.weir.fst", header = T)
head(fst_pe_hy)
#Using qqman 
pe_hy_subset <- fst_pe_hy[complete.cases(fst_pe_hy),]

PE_HY_SNPs <- c(1:(nrow(pe_hy_subset)))

SNP_df <- data.frame(PE_HY_SNPs, pe_pi_subset)

manhattan(SNP_df, chr = "CHROM", bp = "POS", p= "WEIR_AND_COCKERHAM_FST", snp = "PE_HY_SNPs", logp = F, ylab = "Weir and Cockerham Fst", ylim = c(0, 1.1), main = "Pulex-Hybrids Per Site Fst")

#Pulicaria to Hybrid:

fst_pi_hy <- read.table("/home/weider/LinkageMap_Genomics/SlidingWindows/pulicaria_hybrid.weir.fst", header = T)
head(fst_pi_hy)
#Using qqman 
pi_hy_subset <- fst_pi_hy[complete.cases(fst_pi_hy),]

PI_HY_SNPs <- c(1:(nrow(pi_hy_subset)))

SNP_df <- data.frame(PI_HY_SNPs, pi_hy_subset)

manhattan(SNP_df, chr = "CHROM", bp = "POS", p= "WEIR_AND_COCKERHAM_FST", snp = "PI_HY_SNPs", logp = F, ylab = "Weir and Cockerham Fst", ylim = c(0, 1.1), main = "Pulicaria-Hybrid Per Site Fst")

#Window Fst, DXY and Pi measurements:

big_window<- read.table("/home/weider/LinkageMap_Genomics/SlidingWindows/Pulex_Pulicaria_allsnp.slidingWindow.csv", sep = ",", header = T)
head(big_window)  

#Pulex-Pulicaria:

CHR01 <- subset(big_window, big_window$scaffold == "CHR01")

ggplot(CHR01, aes(x=mid, y=pi_Daphnia_pulicaria))+
  geom_line(color = "black", linetype = 1)+
  ylab("Average Pairwise Nucleotide Differences (pi)")+
  xlab("Window Midpoint (BP)")+
  theme_light()
  

ggplot(CHR01, aes(x=mid, y=dxy_Daphnia_pulicaria_Daphnia_pulex))+
  geom_line(color = "black", linetype = 1)+
  ylab("Absolute Nucleotide Differences (Dxy)")+
  xlab("Window Midpoint (BP)")+
  theme_light()

ggplot(CHR01, aes(x=mid, y=Fst_Daphnia_pulicaria_Daphnia_pulex))+
  geom_line(color = "black", linetype = 1)+
  ylab("Fixation Index (Fst)")+
  xlab("Window Midpoint (BP)")+
  theme_light()
