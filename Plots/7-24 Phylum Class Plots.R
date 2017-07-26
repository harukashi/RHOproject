

library(ggplot2)

setwd("C:/Users/student/Google Drive/URI Grad/SAMSI/Data/SAMSI_2017")

#read data sets
#create relative frequency variable

dust_class<-readRDS("df_dust_class.rds")
dust_class$rel_freq<-dust_class$Freq/sum(dust_class$Freq)
dust_class$source<-"original"

sup_dust_class<-readRDS("df_sup_learn_dust_class.rds")
sup_dust_class$rel_freq<-sup_dust_class$Freq/sum(sup_dust_class$Freq)
sup_dust_class$source<-"supervised"

nasal_class<-readRDS("df_nasal_class.rds")
nasal_class$rel_freq<-nasal_class$Freq/sum(nasal_class$Freq)
nasal_class$source<-"original"

sup_nasal_class<-readRDS("df_sup_learn_nasal_class.rds")
sup_nasal_class$rel_freq<-sup_nasal_class$Freq/sum(sup_nasal_class$Freq)
sup_nasal_class$source<-"supervised"

combined_class_dust<-rbind(dust_class, sup_dust_class)
combined_class_nasal<-rbind(nasal_class, sup_nasal_class)


dust_phylum<-readRDS("df_dust_phylum.rds")
dust_phylum$rel_freq<-dust_phylum$Freq/sum(dust_phylum$Freq)
dust_phylum$source<-"original"

sup_dust_phylum<-readRDS("df_sup_learn_dust_phylum.rds")
sup_dust_phylum$rel_freq<-sup_dust_phylum$Freq/sum(sup_dust_phylum$Freq)
sup_dust_phylum$source<-"supervised"

nasal_phylum<-readRDS("df_nasal_phylum.rds")
nasal_phylum$rel_freq<-nasal_phylum$Freq/sum(nasal_phylum$Freq)
nasal_phylum$source<-"original"

sup_nasal_phylum<-readRDS("df_sup_learn_nasal_phylum.rds")
sup_nasal_phylum$rel_freq<-sup_nasal_phylum$Freq/sum(sup_nasal_phylum$Freq)
sup_nasal_phylum$source<-"supervised"

combined_phylum_dust<-rbind(dust_phylum, sup_dust_phylum)
combined_phylum_nasal<-rbind(nasal_phylum, sup_nasal_phylum)

combined_phylum_nasal_dust<-rbind(nasal_phylum, dust_phylum)
combined_class_nasal_dust<-rbind(nasal_class, dust_class)

#barplots

#dust phylum
dp<-ggplot(data=combined_phylum_dust[which(combined_phylum_dust$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle("Dust Phylum Distribution", subtitle = "Dust OTUs vs. Supervised Learning Subset") +
  labs(x="Phylum (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

#dust class
dc<-ggplot(data=combined_class_dust[which(combined_class_dust$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle("Dust Class Distribution", subtitle = "Dust OTUs vs. Supervised Learning Subset") +
  labs(x="Class (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))


#nasal phylum
np<-ggplot(data=combined_phylum_nasal[which(combined_phylum_nasal$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle("Nasal Phylum Distribution", subtitle = "Nasal OTUs vs. Supervised Learning Subset") +
  labs(x="Phylum (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

#nasal class
nc<-ggplot(data=combined_class_nasal[which(combined_class_nasal$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle(label="Nasal Class Distribution", subtitle = "Nasal OTUs vs. Supervised Learning Subset") + 
  labs(x="Class (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

#nasal vs dust
#class
ndc<-ggplot(data=combined_class_nasal_dust[which(combined_class_nasal_dust$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle(label="Class Distribution", subtitle = "Nasal vs. Dust") + 
  labs(x="Class (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

ndp<-ggplot(data=combined_phylum_nasal_dust[which(combined_phylum_nasal_dust$rel_freq>0.01),], aes(x=reorder(Var1, -rel_freq), y=rel_freq, fill=source)) +
  geom_bar(stat='identity', position='dodge') + ggtitle(label="Phylum Distribution", subtitle = "Nasal vs. Dust") + 
  labs(x="Phylum (with relative frequency > 0.01)", y="Relative Frequency (within source)") +theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

png("dust_phylum.png", width = 4, height = 4, units = 'in', res = 300)
print(dp)
dev.off()

png("dust_class.png", width = 4, height = 4, units = 'in', res = 300)
print(dc)
dev.off()

png("nasal_phylum.png", width = 4, height = 4, units = 'in', res = 300)
print(np)
dev.off()

png("nasal_class.png", width = 4, height = 4, units = 'in', res = 300)
print(nc)
dev.off()

png("nasal_dust_class.png", width = 4, height = 4, units = 'in', res = 300)
print(ndc)
dev.off()

png("nasal_dust_phylum.png", width = 4, height = 4, units = 'in', res = 300)
print(ndp)
dev.off()
