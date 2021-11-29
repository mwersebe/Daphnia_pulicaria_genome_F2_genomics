## South Center Genome
## Permutation of background SNP recombination rates:
## Matthew Wersebe, Dept. of Biology, University of Oklahoma
## November 19, 2021
## ############################################################################

# Store each set of background SNPs in a list containing the matrix of data

permutations <- vector(mode = "list", length = 1000)

# All the background SNPs from the persite Fst scan:

background <- subset(pe_pi_singlesnp_fst, pe_pi_singlesnp_fst$outlier == "background")

# For loop to fill the loop with new permutations:

for (i in 1:1000){
  permutations[[i]] <- sample_n(background, 762)
}

# Store the results in a new list:

perm_recomb <- vector(mode = "list", length = 1000)

# for loop to fill this new list: 

for (i in 1:1000){
  
  chr01 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_01, Chromosome = 1)
  
  chr02 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_02, Chromosome = 2)
  
  chr03 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_02, Chromosome = 2)
  
  chr04 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_04, Chromosome = 4)
  
  chr05 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_05, Chromosome = 5)
  
  chr06 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_06, Chromosome = 6)
  
  chr07 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_07, Chromosome = 7)
  
  chr08 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_08, Chromosome = 8)
  
  chr09 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_09, Chromosome = 9)
  
  chr10 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_10, Chromosome = 10)
  
  chr11 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_11, Chromosome = 11)
  
  chr12 <- Outlier_recomb_rate(Outliers = permutations[[i]], Marey_Map = Marey_map_12, Chromosome = 12)
  #make a data frame from results:
  perm_recomb[[i]] <- rbind.data.frame(chr01, chr02, chr03, chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12)
  
  names(perm_recomb[[i]])[5] <- "male_rate"
  names(perm_recomb[[i]])[6] <- "female_rate"
  
}



male_means <- vector(mode = "numeric", length = 1000)
female_means <- vector(mode = "numeric", length = 1000)
for (i in 1:1000){
  perm_recomb[[i]][mapply(is.infinite, perm_recomb[[i]])] <- NA
  male <- perm_recomb[[i]][,5]
  female <- perm_recomb[[i]][,6]
  
  male_means[i] <- mean(male, na.rm = T)
  female_means[i] <- mean(female, na.rm = T)
}

male_means <- as.data.frame(male_means)
female_means <- as.data.frame(female_means)


male_test <- ggplot(male_means, aes(x=male_means))+
  geom_histogram(color="darkblue", fill="lightblue")+
  xlab("Mean Recombination Rate \n (100 KB windows Centered on SNP)")+
  ylab("Count")+
  theme_classic()+
  geom_vline(xintercept = males_maps, color = "red")+
  ggtitle("Male Map")


female_test <- ggplot(female_means, aes(x=female_means))+
  geom_histogram(color="darkblue", fill="lightblue")+
  xlab("Mean Recombination Rate \n (100 KB windows Centered on SNP)")+
  ylab("Count")+
  theme_classic()+
  geom_vline(xintercept = female_maps, color = "red")+
  ggtitle("Female Map")

outlier_male <- recomb_outliers[,5]
outlier_male

outlier_female <- recomb_outliers[,6]
outlier_female

males_maps <- mean(outlier_male, na.rm = T)
female_maps <- mean(outlier_female, na.rm = T)

quantile(male_means$male_means, c(0.975, 0.995), na.rm = T)

length(male_means[male_means>7.523678])

length(female_means[female_means>7.52986])


### CDS Content:
set.seed(NULL)
perm_cds <- vector(mode = "list", length = 1000)
for (i in 1:1000){
  
  CHR01 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 01, Chromosome_gff = "CHR01", Window_size = 100000)
  
  CHR02 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 02, Chromosome_gff = "CHR02", Window_size = 100000)
  
  CHR03 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 03, Chromosome_gff = "CHR03", Window_size = 100000)
  
  CHR04 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 04, Chromosome_gff = "CHR04", Window_size = 100000)
  
  CHR05 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 05, Chromosome_gff = "CHR05", Window_size = 100000)
  
  CHR06 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 06, Chromosome_gff = "CHR06", Window_size = 100000)
  
  CHR07 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 07, Chromosome_gff = "CHR07", Window_size = 100000)
  
  CHR08 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 08, Chromosome_gff = "CHR08", Window_size = 100000)
  
  CHR09 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 09, Chromosome_gff = "CHR09", Window_size = 100000)
  
  CHR10 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 10, Chromosome_gff = "CHR10", Window_size = 100000)
  
  CHR11 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 11, Chromosome_gff = "CHR11", Window_size = 100000)
  
  CHR12 <- Outlier_gene_content(Outlier= permutations[[i]], CDS_gff = CDS_gff, Chromosome = 12, Chromosome_gff = "CHR12", Window_size = 100000)
  
 perm_cds[[i]] <- rbind.data.frame(CHR01, CHR02, CHR03, CHR04, CHR05, CHR06, CHR07, CHR08, CHR09, CHR10, CHR11, CHR12)
 
}
for (i in 1:1000){
  names(perm_cds[[i]])[5] <- "CDS_content"
}

cds_means <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  perm_cds[[i]][mapply(is.infinite, perm_cds[[i]])] <- NA
  cds <- perm_cds[[i]][,5]
  cds_means[i] <- mean(cds, na.rm = T)
}

cds_means <- as.data.frame(cds_means)


out_cds <-cds_outliers[,5]
out_cds
mean_out_cds <- mean(out_cds)
mean_out_cds
cds_test <- ggplot(cds_means, aes(x=cds_means))+
  geom_histogram(color="darkblue", fill="lightblue")+
  xlab("Mean CDS Content \n (100 KB windows Centered on SNP)")+
  ylab("Count")+
  theme_classic()+
  geom_vline(xintercept = mean_out_cds , color = "red")+
  ggtitle("Coding Sequence Content")
cds_test
