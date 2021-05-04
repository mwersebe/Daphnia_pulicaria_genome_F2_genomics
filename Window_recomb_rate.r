# Windowed Recombination Rates: 
# Custom R script by Matthew Wersebe
# Department of Biology, University of Oklahoma
# March 18,2021

#Function Recomb_rate: Returns a data-frame with physical position and recombination rate in cM/MB for plotting. 

Recomb_rate <- function(Data, Window_Size, Step_Size){
#1: Get the biggest value in the Data$physical position:

physical_length <- max(Data$Physical_position)
physical_length

#2: Get the number of windows to the nearest half length
number_windows = (physical_length/Step_size)

Windows = seq(from = 1, to = number_windows, by = 1)

rates = matrix(ncol = 3, nrow = length(Windows))

Size_scalar = (Window_size/1000000)

for (i in 1:length(Windows)){
  #determine the window
  Window <- subset(Data, Data$Physical_position > Step_size*(i - 1) & Data$Physical_position < Window_size+(Step_size*(i - 1)))
  #calc rates 
  rate_male <- (((max(Window$cM_male))-(min(Window$cM_male)))/(Size_scalar))
  rate_female <- (((max(Window$cM_female))-(min(Window$cM_female)))/(Size_scalar))
  #print the rate
  rates[i,1] <- rate_male
  rates[i,2] <- rate_female
  rates[i,3] <- Step_size*(i)
}
recombination <- as.data.frame(rates)

return(recombination)
}
#Example Use

Data <- read.table("SC_map_Chrom01-marey.txt", header=T)

Window_size = 1000000
Step_size = 10000

recombination <- Recomb_rate(Data = Data, Window_Size = Window_size, Step_Size = Step_size)
head(recombination)

recombination %>% 
  mutate_if(is.numeric, list(~na_if(., "-Inf")))
recombination %>% 
  mutate_if(is.numeric, list(~na_if(., "Inf")))

#Plot
plot <- ggplot(recombination, aes(x=V3)) + 
  geom_line(aes(y = V1, color= "steelblue"))+
  geom_line(aes(y = V2, color = "darkred"))+
  xlab("Phyiscal Position (MB)")+ ylab("Recombination Rate (cM/MB)")+
  scale_color_identity(guide = "legend", labels = c("Female", "Male"))
plot







