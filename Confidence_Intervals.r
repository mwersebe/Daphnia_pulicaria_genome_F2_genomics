#Demographic Analysis: 95% Confidence Intervals for Parameter Estimates
#Matthew Wersebe, University of Oklahoma
#05/25/2021
###############################################################################

estimates <- read.table("Demographic_bootstrap_estimates.txt", header = F, sep = "\ ")
head(estimates)
#Extract 95% CI from output:

#Nanc
t.test(estimates$V3, y = NULL, alternative = "two.sided")

#NPa
t.test(estimates$V4, y = NULL, alternative = "two.sided")

#NPx
t.test(estimates$V5, y = NULL, alternative = "two.sided")

#TDiv
t.test(estimates$V6, y = NULL, alternative = "two.sided")

#Mig21
t.test(estimates$V7, y = NULL, alternative = "two.sided")

#Mig12 
t.test(estimates$V8, y = NULL, alternative = "two.sided")
