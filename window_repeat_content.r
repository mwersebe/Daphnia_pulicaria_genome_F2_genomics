
# March 18,2021
##############################################################################
#Plots and Tests:

Phys_position <- CHR01_repeat$V1
total_repeat_content <- CHR01_repeat$V2
LTR_content <- CHR01_repeat$V3
CDS_content <- CHR01_repeat$V7


male_window_recomb_rate <- recombination_1$V1
female_window_recomb_rate <- recombination_1$V2

Rep_Recomb_correlation_1 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate, female_window_recomb_rate)
head(Rep_Recomb_correlation_1)


rate.plot_m_rep <-ggplot(Rep_Recomb_correlation, aes(x=total_repeat_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map Repeat Relationship")+
  xlab("Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)+theme_light()
rate.plot_m_rep

rate.plot_f_rep <-ggplot(Rep_Recomb_correlation, aes(x=total_repeat_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("Female Map Repeat Relationship")+
  xlab("Repeat Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)+theme_light()
rate.plot_f_rep

rate.plot_f_cds <-ggplot(Rep_Recomb_correlation, aes(x=CDS_content, y=female_window_recomb_rate))+
  geom_point(color= "darkred")+
  ggtitle("Female Map CDS Relationship")+
  xlab("CDS Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)+theme_light()
rate.plot_f_cds

rate.plot_m_cds <-ggplot(Rep_Recomb_correlation, aes(x=CDS_content, y=male_window_recomb_rate))+
  geom_point(color= "steelblue")+
  ggtitle("Male Map CDS Relationship")+
  xlab("CDS Content (Proportion)")+ ylab("Recombination Rate (cM/MB)")+
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)+theme_light()
rate.plot_m_cds


ggarrange(rate.plot_m_rep, rate.plot_f_rep, rate.plot_m_cds, rate.plot_f_cds, nrow = 2, ncol = 2, labels =c("A)", "B)", "C)", "D)"))

#Chromosome 1:

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


#################################################################################

#Chromosome 2:

Phys_position <- CHR02_repeat$V1
total_repeat_content <- CHR02_repeat$V2
CDS_content <- CHR02_repeat$V7

male_window_recomb_rate <- recombination_2$V1
female_window_recomb_rate <- recombination_2$V2

Rep_Recomb_correlation_2 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate, female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


#################################################################################

#Chromosome 3:

Phys_position <- CHR03_repeat$V1
total_repeat_content <- CHR03_repeat$V2
CDS_content <- CHR03_repeat$V7

male_window_recomb_rate <- recombination_3$V1
female_window_recomb_rate <- recombination_3$V2
dim(recombination_3)

Rep_Recomb_correlation_3 <-cbind.data.frame(Phys_position[1:1332], total_repeat_content[1:1332], CDS_content[1:1332], male_window_recomb_rate, female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)

################################################################################
#Chromsome 4: 

Phys_position <- CHR04_repeat$V1
total_repeat_content <- CHR04_repeat$V2
CDS_content <- CHR04_repeat$V7
dim(CHR04_repeat)
male_window_recomb_rate <- recombination_4$V1
female_window_recomb_rate <- recombination_4$V2
dim(recombination_4)

getOption("na.action")

Rep_Recomb_correlation_4 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate[1:2248], female_window_recomb_rate[1:2248])

Rep_Recomb_correlation_4[mapply(is.infinite, Rep_Recomb_correlation_4)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$`male_window_recomb_rate[1:2248]`, use="na.omit")

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$`female_window_recomb_rate[1:2248]`)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$`male_window_recomb_rate[1:2248]`)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$`female_window_recomb_rate[1:2248]`)

################################################################################

#chromosome 5:

Phys_position <- CHR05_repeat$V1
total_repeat_content <- CHR05_repeat$V2
CDS_content <- CHR05_repeat$V7
dim(CHR05_repeat)
male_window_recomb_rate <- recombination_5$V1
female_window_recomb_rate <- recombination_5$V2
dim(recombination_5)


Rep_Recomb_correlation_5 <-cbind.data.frame(Phys_position[1:719], total_repeat_content[1:719], CDS_content[1:719], male_window_recomb_rate, female_window_recomb_rate)

Rep_Recomb_correlation_5[mapply(is.infinite, Rep_Recomb_correlation_5)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)

################################################################################

#Chromsome 6:

Phys_position <- CHR06_repeat$V1
total_repeat_content <- CHR06_repeat$V2
CDS_content <- CHR06_repeat$V7
dim(CHR06_repeat)
male_window_recomb_rate <- recombination_6$V1
female_window_recomb_rate <- recombination_6$V2
dim(recombination_6)


Rep_Recomb_correlation_6 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate[1:2418], female_window_recomb_rate[1:2418])

Rep_Recomb_correlation_6[mapply(is.infinite, Rep_Recomb_correlation_6)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


#################################################################################

#Chromosome 7:

Phys_position <- CHR07_repeat$V1
total_repeat_content <- CHR07_repeat$V2
CDS_content <- CHR07_repeat$V7
dim(CHR07_repeat)
male_window_recomb_rate <- recombination_7$V1
female_window_recomb_rate <- recombination_7$V2
dim(recombination_7)


Rep_Recomb_correlation_7 <-cbind.data.frame(Phys_position[1:1884], total_repeat_content[1:1884], CDS_content[1:1884], male_window_recomb_rate, female_window_recomb_rate)

Rep_Recomb_correlation_7[mapply(is.infinite, Rep_Recomb_correlation_7)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


################################################################################

#Chromosome 8:

Phys_position <- CHR08_repeat$V1
total_repeat_content <- CHR08_repeat$V2
CDS_content <- CHR08_repeat$V7
dim(CHR08_repeat)
male_window_recomb_rate <- recombination_8$V1
female_window_recomb_rate <- recombination_8$V2
dim(recombination_8)


Rep_Recomb_correlation_8 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate[1:1585], female_window_recomb_rate[1:1585])

Rep_Recomb_correlation_8[mapply(is.infinite, Rep_Recomb_correlation_8)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


################################################################################

#Chromosome 9:

Phys_position <- CHR09_repeat$V1
total_repeat_content <- CHR09_repeat$V2
CDS_content <- CHR09_repeat$V7
dim(CHR09_repeat)
male_window_recomb_rate <- recombination_9$V1
female_window_recomb_rate <- recombination_9$V2
dim(recombination_9)


Rep_Recomb_correlation_9 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate, female_window_recomb_rate)

Rep_Recomb_correlation_9[mapply(is.infinite, Rep_Recomb_correlation_9)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)

################################################################################

#Chromosome 10:

Phys_position <- CHR10_repeat$V1
total_repeat_content <- CHR10_repeat$V2
CDS_content <- CHR10_repeat$V7
dim(CHR10_repeat)
male_window_recomb_rate <- recombination_10$V1
female_window_recomb_rate <- recombination_10$V2
dim(recombination_10)


Rep_Recomb_correlation_10 <-cbind.data.frame(Phys_position[1:717], total_repeat_content[1:717], CDS_content[1:717], male_window_recomb_rate, female_window_recomb_rate)

Rep_Recomb_correlation_10[mapply(is.infinite, Rep_Recomb_correlation_10)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)

################################################################################

#Chromsome 11:

Phys_position <- CHR11_repeat$V1
total_repeat_content <- CHR11_repeat$V2
CDS_content <- CHR11_repeat$V7
dim(CHR11_repeat)
male_window_recomb_rate <- recombination_11$V1
female_window_recomb_rate <- recombination_11$V2
dim(recombination_11)


Rep_Recomb_correlation_11 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate[1:461], female_window_recomb_rate[1:461])

Rep_Recomb_correlation_11[mapply(is.infinite, Rep_Recomb_correlation_11)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)


#Chromosome 12:

Phys_position <- CHR12_repeat$V1
total_repeat_content <- CHR12_repeat$V2
CDS_content <- CHR12_repeat$V7
dim(CHR12_repeat)
male_window_recomb_rate <- recombination_12$V1
female_window_recomb_rate <- recombination_12$V2
dim(recombination_12)


Rep_Recomb_correlation_12 <-cbind.data.frame(Phys_position, total_repeat_content, CDS_content, male_window_recomb_rate[1:322], female_window_recomb_rate[1:322])

Rep_Recomb_correlation_11[mapply(is.infinite, Rep_Recomb_correlation_11)] <- NA


cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$CDS_content, Rep_Recomb_correlation$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation$total_repeat_content, Rep_Recomb_correlation$female_window_recomb_rate)

######################################################################################
# Overall Correlations:
names(Rep_Recomb_correlation_2) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_3) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_4) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_5) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_6) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_7) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_8) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_9) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_10) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_11) <- names(Rep_Recomb_correlation_1)
names(Rep_Recomb_correlation_12) <- names(Rep_Recomb_correlation_1)
Rep_Recomb_correlation_all <- rbind(Rep_Recomb_correlation_1, Rep_Recomb_correlation_2, Rep_Recomb_correlation_3, Rep_Recomb_correlation_4, Rep_Recomb_correlation_5, Rep_Recomb_correlation_6, Rep_Recomb_correlation_7, Rep_Recomb_correlation_8, Rep_Recomb_correlation_9, Rep_Recomb_correlation_10, Rep_Recomb_correlation_11, Rep_Recomb_correlation_12)

cor.test(Rep_Recomb_correlation_all$CDS_content, Rep_Recomb_correlation_all$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation_all$CDS_content, Rep_Recomb_correlation_all$female_window_recomb_rate)

cor.test(Rep_Recomb_correlation_all$total_repeat_content, Rep_Recomb_correlation_all$male_window_recomb_rate)

cor.test(Rep_Recomb_correlation_all$total_repeat_content, Rep_Recomb_correlation_all$female_window_recomb_rate)
