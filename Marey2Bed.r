# Script to make a bed file from marey map data:

# read in the marey map: 

map <- read.table("~/Downloads/SC_map_Chrom01-marey.txt", header = T, sep = "\t")

Start <- vector(mode = "numeric", length = length(map$Physical_position))

End = map$Physical_position

for (i in 1:length(End)){
  Start[i] = End[i]-1
}


output <- cbind.data.frame(map$Chromosome, Start, End, map$Marker_ID, map$cM_male, map$cM_female)

write.table(output, file = "Marey2bed.txt",col.names = F, sep = "\t", quote = F, row.names = F)
