## South Center Genome
## Gene/ Annotation Averages based on augustus gff3 files.
## Matthew Wersebe, University of Oklahoma, Nov. 17, 2021
############################################################################
library(plotrix)
setwd("/home/weider/Genome_assemblies/Hap13B")
## Haplotype 13B:

## in shell: grep "transcript" <gff3 file> |grep -v "#" > transcript.gff3 
## grep -v "h1tg" transcript.gff3 > placed_transcripts.gff3

transcripts <- read.table("transcript.gff3", header = F)
introns <- read.table("intron.gff3", header =F)
exons <- read.table("exon.gff3", header = F)
placed_transcripts <- read.table("placed_transcripts.gff3", header = F)

#gene lengths 
transcripts_lengths <- transcripts$V5 - transcripts$V4

mean(transcripts_lengths)
std.error(transcripts_lengths)

#introns
intron_lengths <- introns$V5 - introns$V4

mean(intron_lengths)
std.error(intron_lengths)

#exons
exon_lengths <- exons$V5 - exons$V4

mean(exon_lengths)
std.error(exon_lengths)



#Inter gene length:

head(transcripts)

transcripts_end <- placed_transcripts$V5
transcripts_end

inter_gene_length = vector(mode = "numeric", length = length(transcripts_end))

for(i in 1:length(transcripts_end)){
  
  inter_gene_length[i] = placed_transcripts[(i+1), 4] - transcripts_end[i]
}

#check locations of chromosome splits:
inter_gene_length[22520]

#mean/se without chromosome splits:

mean(inter_gene_length[-c(5183, 7038, 8245, 11106, 12144, 15086, 17683, 19743, 20582, 21626, 22520)], na.rm = T)
std.error(inter_gene_length[-c(5183, 7038, 8245, 11106, 12144, 15086, 17683, 19743, 20582, 21626, 22520)])

###############################################################################

##haplotype 2X 

setwd("/home/weider/Genome_assemblies/Hap2X")
## Haplotype 13B:

## in shell: grep "transcript" <gff3 file> |grep -v "#" > transcript.gff3 
## grep -v "h1tg" transcript.gff3 > placed_transcripts.gff3

transcripts <- read.table("transcript.gff3", header = F)
introns <- read.table("intron.gff3", header =F)
exons <- read.table("exon.gff3", header = F)
placed_transcripts <- read.table("placed_transcript.gff3", header = F)

#gene lengths 
transcripts_lengths <- transcripts$V5 - transcripts$V4

mean(transcripts_lengths)
std.error(transcripts_lengths)

#introns
intron_lengths <- introns$V5 - introns$V4

mean(intron_lengths)
std.error(intron_lengths)

#exons
exon_lengths <- exons$V5 - exons$V4

mean(exon_lengths)
std.error(exon_lengths)



#Inter gene length:


transcripts_end <- placed_transcripts$V5
transcripts_end

inter_gene_length = vector(mode = "numeric", length = length(transcripts_end))

for(i in 1:length(transcripts_end)){
  
  inter_gene_length[i] = placed_transcripts[(i+1), 4] - transcripts_end[i]
}

#check locations of chromosome splits:
inter_gene_length[23375]

#mean/se without chromosome splits:

mean(inter_gene_length[-c(4893, 6519, 7849, 10621, 11603, 14856, 17023, 19135, 19808, 20903, 21795)], na.rm = T)
std.error(inter_gene_length[-c(4893, 6519, 7849, 10621, 11603, 14856, 17023, 19135, 19808, 20903, 21795)])

