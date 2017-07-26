
# Step 1: Load Packages Needed for Analysis

library(plyr, quietly = T, warn.conflicts=F, verbose=F) 
library(dplyr, quietly = T, warn.conflicts=F, verbose=F) 
library(ggplot2, quietly = T, warn.conflicts=F, verbose=F)
library(lattice) 
library(stats) 
library(utils) 
library(gridExtra) 
library(knitr)
library(car, quietly = T, warn.conflicts=F, verbose=F)
library(tidyr)
library(stats)
library(xtable)
library(reshape2)
library(graphics)
library(grid)
library(caret, quietly = T, warn.conflicts=F, verbose=F)
library(Hmisc, quietly = T, warn.conflicts=F, verbose=F)
library(vegan, quietly = T, warn.conflicts=F, verbose=F)
library(vegan3d, quietly = T, warn.conflicts=F, verbose=F)

# Step 1: Change WD
setwd("/Volumes/Lexar/SAMSI_2017_Project/Data")

# Step 2: Check files names in WD
list.files()

# Step 3: Read in Data to descriptive objects
master <- read.csv("master.csv", header=T)
taxonomy <- read.csv("MB_SAMSI_taxonomy.csv", header=T)
volcano_OTU <- readRDS("volcano_OTU.rds")

# Step 4: Make edited data frame merging frames
# Here we have subseted the OTU with phylogentic information based
# on significant OTU from volcano plot
#
# We need to match corresponding OTU beteen the two data frames
# and then merge togther 

a <- taxonomy[taxonomy$OTU_ID %in% volcano_OTU$OTU_ID,]
b <- volcano_OTU[volcano_OTU$OTU_ID %in% taxonomy$OTU_ID,]
taxonomy2 <- cbind(a,b) # Used to merge together
taxonomy2 <- taxonomy2[, !duplicated(colnames(taxonomy2))]

# Step 5: Remove unwanted objects for space/speed
rm(list=c("a", "b", "volcano_OTU"))


# Step 6: Write Summary Statistics function
summary.statistics <- function(x) c(Mean = mean(x), Median=median(x), St.Dev = sd(x), 
                                    Minimum=range(x)[1], Maximum=range(x)[2] )




# Step 1: Break up pylogenetic data by different groups
# Here we are going to change a few column names for more 
# descriptive purposes
names(master)[10:13] <- c("Nasal_Score", "Nasal_Score_Cutoff", "Location", "Sample_Source")

names(master)[9] <- "Characterisation"

# Step 2: Group by characteristics 

master2 <- master %>% 
  group_by(Sample_Source, Characterisation) %>% 
  summarise(
    Amount = n(), Avg.faith =mean(Faith), sd(Faith), Avg.Chaol=mean(Chao1), sd(Chao1), Avg.Pielo=mean(Pielou), sd(Pielou) )

# Step 3: Group by Taxonomy general counts
# Remove data frames after table setup
names(taxonomy) <- c("OTU_ID", "Kingdom", "Phylum", "Class")
kingdom_grp <- taxonomy %>% group_by(Kingdom, Phylum) %>% mutate(Count = n())
kingdom_grp2 <- taxonomy %>% group_by(Kingdom, Phylum) %>% tally()
kingdom_grp3 <- taxonomy %>% group_by(Kingdom) %>% tally()
kingdom_grp4 <- taxonomy %>% group_by(Kingdom, Phylum, Class) %>% tally()


# Step 3: Group by Taxonomy general percentage
# Remove data frames after ggplot setup

tax_dust <- taxonomy2[taxonomy2$prevalence=="Dust",]
tax_dust <- tax_dust[complete.cases(tax_dust), ]
names(tax_dust) <- c("OTU_ID", "Kindgom", "Phylum", "Class", "Log_Change", "P_Value", "Prevalence")
tax_nasal <- taxonomy2[taxonomy2$prevalence=="Nasal",]
tax_nasal<- tax_nasal[complete.cases(tax_nasal), ]
names(tax_nasal) <- c("OTU_ID", "Kindgom", "Phylum", "Class", "Log_Change", "P_Value", "Prevalence")

group <- tax_dust %>% group_by(Phylum) %>% tally()
group <- ddply(group, "Phylum", mutate, Phylum_Percent_Weight = n/sum(group$n))
Phylum_Percent_Weight <- rowsum(group$Phylum_Percent_Weight, group$Phylum_Percent_Weight <= 0.010)[2]
n <- rowsum(group$n, group$n < 20)[2]
group <- group[group$Phylum_Percent_Weight > 0.010,]
p <- "Other_Phylum_Combined"
p <- data.frame(cbind(p,n, Phylum_Percent_Weight))
names(p)[1] <- "Phylum"
group1 <- rbind(group, p)
group1$Phylum_Percent_Weight <- as.numeric(group1$Phylum_Percent_Weight)
group1$Phylum_Percent_Weight <- round(group1$Phylum_Percent_Weight, digits=3)
rm(list=c("group", "n", "p", "Phylum_Percent_Weight" ))


group2 <- tax_nasal %>% group_by(Phylum) %>% tally()
group2 <- ddply(group2, "Phylum", mutate, Phylum_Percent_Weight = n/sum(group2$n))

# Step 4: Set up latex table
Summary1 <- kable(group1, format = 'latex', digits=4, row.names =T)
Summary1 <- sub('l|r}', 'l|r|}', Summary1, fixed=T)
Summary1 <- sub('{l|l', '{|l|l', Summary1, fixed=T)

Summary2 <- kable(group2, format = 'latex', digits=4, row.names =T)
Summary2 <- sub('r|r}', 'r|r|}', Summary2, fixed=T)
Summary2 <- sub('{l|l', '{|l|l', Summary2, fixed=T)


# Step 5: Print latex output


cat(c("\\begin{table}[h] \\centering ",
      Summary1,
      "\\hspace{1cm} \\centering ",
      "\\caption{Phylogenitic Summary} \\end{table}"))


cat(c("\\begin{table}[h] \\centering ",
      Summary2,
      "\\hspace{1cm} \\centering ",
      "\\caption{Phylogenitic Summary} \\end{table}"))

rm(list=c("tax_dust","tax_nasal","group1", "group2", "Summary1", "Summary2"))

# Step 1: Make kniter and latex table for descriptive statistics

Summary.Faith <- tapply(master$Faith, master$Sample_Source, FUN=summary.statistics)
Summary.Faith <- do.call(cbind, Summary.Faith)
Summary.Faith <- kable(Summary.Faith, format = 'latex', digits=3, row.names =T)
Summary.Faith <- sub('r|r}', 'r|r|}', Summary.Faith, fixed=T)
Summary.Faith <- sub('{l|r', '{|l|r', Summary.Faith, fixed=T)


Summary.Chao1 <- tapply(master$Chao1, master$Sample_Source, FUN=summary.statistics)
Summary.Chao1 <- do.call(cbind, Summary.Chao1)
Summary.Chao1 <- kable(Summary.Chao1, format = 'latex', digits=3, row.names =T)
Summary.Chao1 <- sub('r|r}', 'r|r|}', Summary.Chao1, fixed=T)
Summary.Chao1 <- sub('{l|r', '{|l|r', Summary.Chao1, fixed=T)

Summary.Pielou <- tapply(master$Pielou, master$Sample_Source, FUN=summary.statistics)
Summary.Pielou <- do.call(cbind, Summary.Pielou)
Summary.Pielou <- kable(Summary.Pielou, format = 'latex', digits=3, row.names =T)
Summary.Pielou <- sub('r|r}', 'r|r|}', Summary.Pielou, fixed=T)
Summary.Pielou <- sub('{l|r', '{|l|r', Summary.Pielou, fixed=T)

rm(summary.statistics)

# Step 2: Print latex output

cat(c("\\begin{table}[h] \\centering ",
      Summary.Faith,
      "\\hspace{1cm} \\centering ",
      "\\caption{Faith Scores} \\end{table}"))


cat(c("\\begin{table}[h] \\centering ",
      Summary.Chao1,
      "\\hspace{1cm} \\centering ",
      "\\caption{Chao1 Scores} \\end{table}"))


cat(c("\\begin{table}[h] \\centering ",
      Summary.Pielou,
      "\\hspace{1cm} \\centering ",
      "\\caption{Pielou Scores} \\end{table}"))


## Taxanomy Information:

# Step 1: Make kniter and latex table for descriptive statistics
Taxa1 <- kable(kingdom_grp3, format = 'latex', digits=3, row.names =T)
Taxa1 <- sub('l|r}', 'l|r|}', Taxa1, fixed=T)
Taxa1 <- sub('{l|l', '{|l|r', Taxa1, fixed=T)

kingdom_grp2 <- kingdom_grp2[complete.cases(kingdom_grp2),]

Taxa2.1 <- kable(kingdom_grp2[1:12,], format = 'latex', digits=3, row.names =T)
Taxa2.2 <- kable(kingdom_grp2[13:24,], format = 'latex', digits=3, row.names =T)
Taxa2.3 <- kable(kingdom_grp2[25:36,], format = 'latex', digits=3, row.names =T)

Taxa2.1 <- sub('l|r}', 'l|r|}', Taxa2.1, fixed=T)
Taxa2.1 <- sub('{l|l', '{|l|l', Taxa2.1, fixed=T)

Taxa2.2 <- sub('l|r}', 'l|r|}', Taxa2.2, fixed=T)
Taxa2.2 <- sub('{l|l', '{|l|l', Taxa2.2, fixed=T)

Taxa2.3 <- sub('l|r}', 'l|r|}', Taxa2.3, fixed=T)
Taxa2.3 <- sub('{l|l', '{|l|l', Taxa2.3, fixed=T)

# Step 2: Print latex output
cat(c("\\begin{table}[h] \\centering ", 
      Taxa2.1,
      "\\hspace{1cm} \\centering ",
      Taxa2.2,
      "\\hspace{1cm} \\centering ",
      Taxa2.3,
      "\\caption{Taxonomy Table} \\end{table}"))  


# Step 3: Remove unwanted objects for space/speed
rm(list=c("kingdom_grp","kingdom_grp2", "kingdom_grp3", "kingdom_grp4") )

# Step 1: Histogram
histo <- ggplot(master, aes(x = Faith)) +
  geom_histogram(aes(fill = ..count..), binwidth = 5) +
  scale_x_continuous(name = "Mean Faith Scores",
                     breaks = seq(0, 175, 25),
                     limits=c(0, 175)) +
  scale_y_continuous(name = "Count") +
  ggtitle("Frequency Histogram of Mean Faith Diversity Index") +
  scale_fill_gradient("Count", low = "blue", high = "red") + 
  theme(plot.title = element_text(hjust = 0.5, family="Times"))

histo1 <- ggplot(master, aes(x = Faith)) +
  geom_histogram(aes(fill = ..count..), binwidth = 5) +
  scale_x_continuous(name = "Mean Faith Scores",
                     breaks = seq(0, 175, 25),
                     limits=c(0, 175)) +
  scale_y_continuous(name = "Count") +
  ggtitle("Frequency Histogram of Mean Faith Diversity Index") +
  scale_fill_gradient("Count", low = "blue", high = "red")+ facet_grid(~Sample_Source)+ 
  theme(plot.title = element_text(hjust = 0.5, family="Times"))

ggsave("histo.png", width=11, height=7, units = "in")
ggsave("histo1.png", width=11, height=7, units = "in")

# Step 1: Sepearte Nasal from Dust
# We want to create a matrix which realates count of
# OTU to Sample_Source

class(master$Sample_Source)
master$Sample_Source <- as.character(master$Sample_Source)
master_dust <- master[master$Sample_Source=="Dust",]
master_nasal <- master[master$Sample_Source=="Nasal",]

# Step 2: Remove unwanted objects for space/speed
rm(list=c("master_dust", "master_nasal") )


# Step 1: Rearrange Data frame; for OTU related to sample type
# Methods 1: Other method found in extra chunks

New.testSet <- master[,grep("^OTU_|data_source", colnames(master))]
data <- data.frame(colSums(New.testSet[New.testSet$data_source=="Nasal",-1]))
colnames(data) <- "Nasal"
data2 <- data.frame(colSums(New.testSet[New.testSet$data_source=="Dust",-1]))
colnames(data2) <- "Dust"

# Step 2: Combine data set and turn into a matrix
data3 <- cbind(data2, data)
data3 <- as.matrix(data3)
data3 <- data3[complete.cases(data3), ]

# Step 3: Remove unwanted objects for space/speed
rm(list=c("data", "data2") )

# Step 1: Rearrange Data frame; for taxonomy relation 

allKingdoms <- dcast(taxonomy, Phylum~Kingdom, fun.aggregate=length)
allKingdoms <- allKingdoms[complete.cases(allKingdoms),]
allPhylum <- as.vector(allKingdoms[,"Phylum"])
allKingdomsMat <- as.matrix(allKingdoms[,-1])
rownames(allKingdomsMat) <- allPhylum

# Step 2: Remove unwanted objects for space/speed
rm(list=c("allKingdoms") )
