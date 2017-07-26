# Katie Dio
# created: 7/20/2017
# last updated: 7/24/2017

# Samsi 2017
# Rho Project

#Set working directory to where the data is located
setwd("C:/Users/student/Google Drive/URI Grad/SAMSI/Data/SAMSI_2017")


#############
# Packages: #
#############

library(dplyr)
library(ggplot2)
library(GGally)

#############
# Read Data #
#############

dust<-readRDS("MB_SAMSI_dust_OTU.rds")
nasal<-readRDS("MB_SAMSI_nasal_OTU.rds")

master<-readRDS("MB_SAMSI_master.rds")
master_cr<-readRDS("master_cr.rds")

##############
# Data Fixes #
##############
#Need to remove all zero OTUs for the bar plot of Nasal and Dust OTU Abundance and Instances
OTU<-rbind(nasal, dust)
OTU_nz<-OTU[,which(colSums(OTU[,2:ncol(OTU)])!=0)]

##########################
# Bar Plots of Dust/Nasal#
##########################

# Nasal and Dust Abundance and Instances of OTUs
counts_nas<-colSums(OTU_nz[1:74,2:ncol(OTU_nz)])
counts_dust<-colSums(OTU_nz[75:148,2:ncol(OTU_nz)])

inst_nas<-colSums(OTU_nz[1:74,2:ncol(OTU_nz)]!=0)
inst_dust<-colSums(OTU_nz[75:148,2:ncol(OTU_nz)]!=0)

jpeg("Unpolished_OTU.jpg")
par(mfrow=c(2,2))
barplot(sqrt(counts_nas), main="Sqrt Abundance in Nasal Samples", xaxt='n', ylim=c(0,750), xlab="OTU", ylab="Sum of Counts")
barplot(inst_nas, main="Instances in Nasal Samples", xaxt='n', ylim=c(0,74), xlab="OTU", ylab="Count of Appearances")

barplot(sqrt(counts_dust), main="Sqrt Abundance in Dust Samples", xaxt='n', ylim=c(0,300), xlab="OTU", ylab="Sum of Counts")
barplot(inst_dust, main="Instances in Dust Samples", xaxt='n', ylim=c(0,74), xlab="OTU", ylab="Count of Appearances")
dev.off()


par(mfrow=c(1,1))

##############
# Data Fixes #
##############
#Using master data file, need studyID, wheeze_atopy_grp, data_source, Chao1, Pielou, Site
dataBP<-data.frame("StudyID"=master$StudyID, "wheeze_atopy_grp"=master$wheeze_atopy_grp, "data_source"=master$data_source, "Chao1"=master$Chao1, "Pielou"=master$Pielou, "Faith"=master$Faith, "site"=master$site)
dustBP<-dataBP[which(dataBP$data_source=="Dust"),]
nasBP<-dataBP[which(dataBP$data_source=="Nasal"),]

#medians
medCHdu<-median(dustBP$Chao1)
medCHna<-median(nasBP$Chao1)

#new variable for high or low richness compared to median
dustBP$Chao1HL<-cut(dustBP$Chao1,
                    breaks = c(-Inf, medCHdu, Inf),
                    labels = c("dust.low", "dust.high"),
                    right = FALSE)

nasBP$Chao1HL<-cut(nasBP$Chao1,
                   breaks = c(-Inf, medCHna, Inf),
                   labels = c("nas.low", "nas.high"),
                   right = FALSE)
trial<-inner_join(nasBP, dustBP, by="StudyID")
trial1<-mutate(trial, chao1HLBoth=interaction(Chao1HL.x,Chao1HL.y))

tab3<-prop.table(table(trial1$chao1HLBoth, trial1$wheeze_atopy_grp.x),2)

tab4<-prop.table(table(dustBP$wheeze_atopy_grp, dustBP$site),2)

###########
# GGPlots #
###########

#Boxplot of Chao1 by Wheeze Atopy Group
bp_1<-ggplot(dustBP, aes(x=factor(wheeze_atopy_grp), y=Chao1, fill=wheeze_atopy_grp))+geom_boxplot()+
  ggtitle("Boxplot of Dust Richness by Atopy Group")+labs(x="Wheeze Atopy Group", y="Chao1")+
  theme(legend.position="none")+theme(plot.title = element_text(hjust = 0.5))

bp_2<-ggplot(nasBP, aes(x=factor(wheeze_atopy_grp), y=Chao1, fill=wheeze_atopy_grp))+geom_boxplot()+
  ggtitle("Boxplot of Nasal Richness by Atopy Group")+labs(x="Wheeze Atopy Group", y="Chao1")+
  theme(legend.position="none")+theme(plot.title = element_text(hjust = 0.5))

#Boxplot by Site
bp_3<-ggplot(dustBP, aes(x=factor(site), y=Chao1, fill=site))+geom_boxplot()+
  ggtitle("Boxplot of Dust Richness by Site")+labs(x="Site", y="Chao1")+
  theme(legend.position="none")+theme(plot.title = element_text(hjust = 0.5))

bp_4<-ggplot(nasBP, aes(x=factor(site), y=Chao1, fill=site))+geom_boxplot()+
  ggtitle("Boxplot of Nasal Richness by Site")+labs(x="Site", y="Chao1")+
  theme(legend.position="none")+theme(plot.title = element_text(hjust = 0.5))

#Boxplot by Nasal Dust
bp_5<-ggplot(master, aes(x=factor(data_source), y=sqrt(Chao1), fill=data_source))+geom_boxplot()+
  ggtitle("Boxplot of Sq. Root Richness by Source")+labs(x="Source", y="Sq. Root of Chao1")+
  theme(legend.position="none")+theme(plot.title = element_text(hjust = 0.5))

jpeg("Box_Dust_AtGRP.jpeg")
print(bp_1)
dev.off()

jpeg("Box_Nas_AtGRP.jpeg")
print(bp_2)
dev.off()

jpeg("Box_Dust_SITE.jpeg")
print(bp_3)
dev.off()

jpeg("Box_Nas_SITE.jpeg")
print(bp_4)
dev.off()

jpeg("Box_SOURCE.jpeg")
print(bp_5)
dev.off()


#Barplot of all vars together by Atopy Group
trial2<-as.data.frame(tab3)
colnames(trial2)<-c("Richness_Level", "Atopy_Status", "Rel_Freq")

Bp<-ggplot(data=trial2, aes(x=factor(Atopy_Status), y=Rel_Freq, fill=Richness_Level)) +
  geom_bar(stat='identity') + ggtitle("Distribution of Richness Level by Atopy Status") +
  labs(x="Atopy Status", y="Relative Frequency") +theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(labels=c("Nasal Low Dust Low", "Nasal High Dust Low", "Nasal Low Dust High",
                             "Nasal High Dust High"), values=c("#F8766D", "#7CAE00","#00BFC4", "#C77CFF"))+
  geom_text(aes(label=round(Rel_Freq,3)),position=position_stack(vjust=0.5))

jpeg("Bar_Atopy.jpeg")
print(Bp)
dev.off()

#Barplot of site and atopy
trial3<-as.data.frame(tab4)
colnames(trial3)<-c("Atopy_Status","Site", "Rel_Freq")

BP2<-ggplot(data=trial3, aes(x=factor(Site), y=Rel_Freq, fill=Atopy_Status)) +
  geom_bar(stat='identity') + ggtitle("Distribution of Atopy Status by Site") +
  labs(x="Site", y="Relative Frequency") +theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(labels=c("Atopy Only", "Both", "Neither",
                             "Wheeze Only"), values=c("#F8766D", "#7CAE00","#00BFC4", "#C77CFF"))+
  geom_text(aes(label=round(Rel_Freq,3)),position=position_stack(vjust=0.5))

jpeg("Bar_Site.jpeg")
print(BP2)
dev.off()


# Diversity Measures Correlation Plots
a<-ggpairs(master[which(master$data_source=="Dust"),13:15], title="Diversity Measures for Dust Samples")+
  theme(plot.title = element_text(hjust = 0.5)) #diversity measures
b<-ggpairs(master[which(master$data_source=="Nasal"),13:15], title="Diversity Measures for Nasal Samples")+
  theme(plot.title = element_text(hjust = 0.5))


#Correlation Plot between Dust and Nasal
trial_pairs<-trial
trial_pairs$"Dust.Chao1" <-trial$Chao1.y
trial_pairs$"Dust.Pielou" <-trial$Pielou.y
trial_pairs$"Dust.Faith" <-trial$Faith.y

trial_pairs$"Nasal.Chao1" <-trial$Chao1.x
trial_pairs$"Nasal.Pielou" <-trial$Pielou.x
trial_pairs$"Nasal.Faith" <-trial$Faith.x

corplot<-ggcorr(trial_pairs[,16:21], label=TRUE,  method=c("everything","spearman"),
       label_round=2, color="grey50", hjust=0.50, label_color = "black")+
ggtitle("Spearman Correlation of Diversity Measures")+
  theme(plot.title = element_text(hjust = 0.50))


jpeg("Dust_Div.jpeg")
print(a)
dev.off()

jpeg("Nasal_Div.jpeg")
print(b)
dev.off()

jpeg("Div_Cor.jpeg")
print(corplot)
dev.off()

# Barplot of Signficant OTUs

#Using master_cr data file, need studyID, wheeze_atopy_grp, data_source, Chao1, Pielou, Site
dataBPcr<-data.frame("StudyID"=master_cr$StudyID, "wheeze_atopy_grp"=master_cr$wheeze_atopy_grp, "data_source"=master_cr$data_source, "Chao1"=master_cr$Chao1, "Pielou"=master_cr$Pielou, "site"=master_cr$site)
dustBPcr<-dataBP[which(dataBPcr$data_source=="Dust"),]
nasBPcr<-dataBP[which(dataBPcr$data_source=="Nasal"),]

#medians
medCHducr<-median(dustBPcr$Chao1)
medCHnacr<-median(nasBPcr$Chao1)

#new variable for high or low richness compared to median
dustBPcr$Chao1HL<-cut(dustBPcr$Chao1,
                    breaks = c(-Inf, medCHducr, Inf),
                    labels = c("dust.low", "dust.high"),
                    right = FALSE)

nasBPcr$Chao1HL<-cut(nasBPcr$Chao1,
                   breaks = c(-Inf, medCHnacr, Inf),
                   labels = c("nas.low", "nas.high"),
                   right = FALSE)
trialcr<-inner_join(nasBPcr, dustBPcr, by="StudyID")
trial1cr<-mutate(trialcr, chao1HLBoth=interaction(Chao1HL.x,Chao1HL.y))

tab3cr<-prop.table(table(trial1cr$chao1HLBoth, trial1cr$wheeze_atopy_grp.x),2)

trial2cr<-as.data.frame(tab3cr)
colnames(trial2cr)<-c("Richness_Level", "Atopy_Status", "Rel_Freq")

Bpcr<-ggplot(data=trial2cr, aes(x=factor(Atopy_Status), y=Rel_Freq, fill=Richness_Level)) +
  geom_bar(stat='identity') + ggtitle("Distribution of Richness Level by Atopy Status") +
  labs(x="Atopy Status", y="Relative Frequency") +theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(labels=c("Nasal Low Dust Low", "Nasal High Dust Low", "Nasal Low Dust High",
                             "Nasal High Dust High"), values=c("#F8766D", "#7CAE00","#00BFC4", "#C77CFF"))+
  geom_text(aes(label=round(Rel_Freq,3)),position=position_stack(vjust=0.5))

jpeg("Bar_Atopycr.jpeg")
print(Bpcr)
dev.off()
