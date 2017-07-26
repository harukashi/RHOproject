# Katie Dio
# created: 7/22/2017
# last updated: 7/23/2017

# Samsi 2017
# Rho Asthma Project

#Set working directory to where the data is located
setwd("C:/Users/student/Google Drive/URI Grad/SAMSI/Data/SAMSI_2017")


#############
# Packages: #
#############

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2"); library(DESeq2) # significance testing
biocLite("ggplot2"); library(ggplot2) # best plotting package
biocLite("ggrepel"); library(ggrepel) # ggplot2 addon for better labels
biocLite("dplyr"); library(dplyr) # ggplot2 addon for better labels

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr) 

#############
# Read Data #
#############

master<-readRDS("MB_SAMSI_master.rds")

#create dataset for just the dust or nasal samples
master_dust<-master[which(master$data_source=="Dust"),]
master_nasal<-master[which(master$data_source=="Nasal"),]

# create a variable for neither= no symptoms and atleastone= wheeze, atopy or both
master_dust$neither_other<-ifelse(master_dust$wheeze_atopy_grp=="neither", c("neither"),c("atleastone"))
master_nasal$neither_other<-ifelse(master_nasal$wheeze_atopy_grp=="neither", c("neither"),c("atleastone"))


#################
# Volcano Plots # Response (none vs at least one)
#################

#function for geo means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

########
# DUST #
########

# matrix of OTU values
mat_dust=t(master_dust[,16:10133])
dim(mat_dust)
colnames(mat_dust)=paste0(factor(master_dust$neither_other),gsub("^subject_", "_", master_dust$StudyID))
mat_dust.coldata=matrix(factor(master_dust$neither_other))
colnames(mat_dust.coldata)="condition"
rownames(mat_dust.coldata)=paste0(master_dust$neither_other,gsub("^subject_", "_", master_dust$StudyID))
all(rownames(mat_dust.coldata) %in% colnames(mat_dust))
all(rownames(mat_dust.coldata) == colnames(mat_dust))


# DEseq2 package
dds_dust = DESeqDataSetFromMatrix(countData = mat_dust,
                              colData = mat_dust.coldata,
                              design = ~ condition)


dds_dust = dds_dust[ rowSums(counts(dds_dust)) > 0, ]


geoMeans_dust = apply(counts(dds_dust), 1, gm_mean)
dds_dust = estimateSizeFactors(dds_dust, geoMeans = geoMeans_dust)
dds_dust = DESeq(dds_dust)

# dds = DESeq(dds)
results_dust = results(dds_dust, contrast=c("condition", "atleastone", "neither"))


summary(results_dust$padj)
as.data.frame(results_dust@listData)



data_dust = data.frame(gene = row.names(results_dust),
                   pvalue = -log10(results_dust$padj), 
                   lfc = results_dust$log2FoldChange)

data_dust = na.omit(data_dust)

data_dust = data_dust %>%
  mutate(color = ifelse(data_dust$lfc > 0 & data_dust$pvalue > 1.3, 
                        yes = "atleastone", 
                        no = ifelse(data_dust$lfc < 0 & data_dust$pvalue > 1.3, 
                                    yes = "neither", 
                                    no = "none")))
#### Plot ####
colored_dust = ggplot(data_dust, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle = "Dust OTU Samples") +  
  xlab(expression(log[2]("Response" / "Neither"))) + 
  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = 0, colour = "black") + 
  geom_hline(yintercept = 1.3, colour = "black") + 
  scale_color_manual(values = c("atleastone" = "#E64B35", 
                                 "neither" = "#3182bd", 
                                 "none" = "#636363"),
                      name=NULL,
                      breaks=c("atleastone", "neither"),
                      labels=c("Response", "Neither"))+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Subset table to only show certain gene labels
top_labelled_dust <- top_n(data_dust, n = 10, wt = pvalue)

# Add layer of text annotation to volcano plot.
dust_volcano<-colored_dust + geom_text_repel(data = top_labelled_dust, 
                           mapping = aes(label = gene), 
                           size = 3,
                           fontface = 'bold', 
                           color = 'black',
                           box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.5, "lines")) + 
  scale_y_continuous(trans = "log1p") 

#### Save Plot ####
jpeg("Dust_Super_Volcano.jpeg")
print(dust_volcano)
dev.off()


#### Save Significant OTU ####
OTU_volc_dust<-data_dust[data_dust$color=="neither",]
OTU_volc_dust<-rbind(OTU_volc_dust,data_dust[data_dust$color=="atleastone",])


saveRDS(OTU_volc_dust, "OTU_volc_dust_super.rds")

###############################################################################################

#########
# nasal #
#########

# matrix of nasal OTU
mat_nasal=t(master_nasal[,16:10133])

dim(mat_nasal)
colnames(mat_nasal)=paste0(factor(master_nasal$neither_other),gsub("^subject_", "_", master_nasal$StudyID))
mat_nasal.coldata=matrix(factor(master_nasal$neither_other))
colnames(mat_nasal.coldata)="condition"
rownames(mat_nasal.coldata)=paste0(master_nasal$neither_other,gsub("^subject_", "_", master_nasal$StudyID))
all(rownames(mat_nasal.coldata) %in% colnames(mat_nasal))
all(rownames(mat_nasal.coldata) == colnames(mat_nasal))

# Package DEseq2
dds_nasal = DESeqDataSetFromMatrix(countData = mat_nasal,
                              colData = mat_nasal.coldata,
                              design = ~ condition)


dds_nasal = dds_nasal[ rowSums(counts(dds_nasal)) > 0, ]


geoMeans_nasal = apply(counts(dds_nasal), 1, gm_mean)
dds_nasal = estimateSizeFactors(dds_nasal, geoMeans = geoMeans_nasal)
dds_nasal = DESeq(dds_nasal)

results_nasal = results(dds_nasal, contrast=c("condition", "atleastone", "neither"))


summary(results_nasal$padj)
as.data.frame(results_nasal@listData)



data_nasal = data.frame(gene = row.names(results_nasal),
                   pvalue = -log10(results_nasal$padj), 
                   lfc = results_nasal$log2FoldChange)

data_nasal = na.omit(data_nasal)

data_nasal = data_nasal %>%
  mutate(color = ifelse(data_nasal$lfc > 0 & data_nasal$pvalue > 1.3, 
                        yes = "atleastone", 
                        no = ifelse(data_nasal$lfc < 0 & data_nasal$pvalue > 1.3, 
                                    yes = "neither", 
                                    no = "none")))

#### Plot ####
colored_nasal = ggplot(data_nasal, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle = "Nasal OTU Samples") +  
  xlab(expression(log[2]("Response" / "Neither"))) + 
  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = 0, colour = "black") + 
  geom_hline(yintercept = 1.3, colour = "black") + 
  scale_color_manual(values = c("atleastone" = "#E64B35", 
                                "neither" = "#3182bd", 
                                "none" = "#636363"),
                     name=NULL,
                     breaks=c("atleastone", "neither"),
                     labels=c("Response", "Neither")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Plot

# Subset table to only show certain gene labels
top_labelled_nasal <- top_n(data_nasal, n = 2, wt = pvalue)

# Add layer of text annotation to volcano plot.
nasal_volcano<-colored_nasal + geom_text_repel(data = top_labelled_nasal, 
                           mapping = aes(label = gene), 
                           size = 3,
                           fontface = 'bold', 
                           color = 'black',
                           box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.5, "lines")) + 
  scale_y_continuous(trans = "log1p") 
  
#### Save Plot ####
jpeg("Nasal_Super_Volcano.jpeg")
print(nasal_volcano)
dev.off()

#### Save Significant OTU ####
OTU_volc_nasal<-data_nasal[data_nasal$color=="neither",]
OTU_volc_nasal<-rbind(OTU_volc_nasal,data_nasal[data_nasal$color=="atleastone",])


saveRDS(OTU_volc_nasal, "OTU_volc_nasal_super.rds")
