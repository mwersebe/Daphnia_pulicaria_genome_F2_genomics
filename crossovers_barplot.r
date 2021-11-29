## South Center Genome
## Crossovers per chromosome plots.
## Matthew Wersebe, University of Oklahoma, Nov. 18, 2021

#############################################################################
library(ggplot2)
library(plyr)

setwd("/home/weider/Linkage/LepAnchor_LOD26-13B/Recombination_stderr")
crossovers <- read.csv("crossovers.csv", header = T)
head(crossovers)

chr01 <- count(crossovers, "chr01")
chr12

chr02 <- count(crossovers, "chr02")
chr03 <- count(crossovers, "chr03")
chr04 <- count(crossovers, "chr04")
chr05 <- count(crossovers, "chr05")
chr06 <- count(crossovers, "chr06")
chr07 <- count(crossovers, "chr07")
chr08 <- count(crossovers, "chr08")
chr09 <- count(crossovers, "chr09")
chr10 <- count(crossovers, "chr10")
chr11 <- count(crossovers, "chr11")
chr12 <- count(crossovers, "chr12")


plot_data <- read.csv("crossover_freq.csv", header = T)

head(plot_data)

ggplot(plot_data, aes(fill=Crossovers, y=Frequency, x=Chromosome)) +
  geom_bar(position = "stack", stat = "identity")+
  ylab("Count of Individuals")
