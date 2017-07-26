source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2"); library(DESeq2) # significance testing
biocLite("ggplot2"); library(ggplot2) # best plotting package
biocLite("ggrepel"); library(ggrepel) # ggplot2 addon for better labels
biocLite("dplyr"); library(dplyr) # ggplot2 addon for better labels

# biocLite("pasilla"); library("pasilla")
# 
# data("pasillaGenes")
# 
# pasCts <- system.file("extdata",
#                       "pasilla_gene_counts.tsv",
#                       package="pasilla", mustWork=TRUE)
# pasAnno <- system.file("extdata",
#                        "pasilla_sample_annotation.csv",
#                        package="pasilla", mustWork=TRUE)
# cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
# coldata <- read.csv(pasAnno, row.names=1)
# coldata <- coldata[,c("condition","type")]
# 
# head(cts)
# head(coldata)
# 
# rownames(coldata) <- sub("fb","",rownames(coldata))
# all(rownames(coldata) %in% colnames(cts))
# 
# cts <- cts[, rownames(coldata)]
# all(rownames(coldata) == colnames(cts))
# 
# library("DESeq2")
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ condition)
# dds

mat=t(master[,16:10133])
dim(mat)
colnames(mat)=paste0(master$data_source,gsub("^subject_", "_", master$StudyID))
mat.coldata=matrix(master$data_source)
colnames(mat.coldata)="condition"
rownames(mat.coldata)=paste0(master$data_source,gsub("^subject_", "_", master$StudyID))
all(rownames(mat.coldata) %in% colnames(mat))
all(rownames(mat.coldata) == colnames(mat))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = mat.coldata,
                              design = ~ condition)
dds

dds <- dds[ rowSums(counts(dds)) > 0, ]

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds)

# dds <- DESeq(dds)
results <- results(dds)
results

summary(results$padj)
as.data.frame(results@listData)



# Add gene names
# Add the -log10 pvalue 
# Add the pre-calculated log2 fold change
data <- data.frame(gene = row.names(results),
                   pvalue = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove rows that have NA values
data <- na.omit(data)

# View table
head(data)

simple <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(size = 3, alpha = 0.7, na.rm = T) + # Make dots bigger
  theme_bw(base_size = 16) + # change theme
  ggtitle(label = "Volcano Plot", subtitle = "Simple black & white") + # Add a title
  xlab(expression(log[2]("Nasal" / "Dust"))) + # x-axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
  geom_vline(xintercept = c(-2,2), colour = "darkgrey") + # Add cutoffs
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  scale_colour_gradient(low = "black", high = "black", guide = FALSE) + # Color black
  scale_x_continuous(limits = c(-4, 4)) # min/max of lfc

# Plot figure
simple

# Modify dataset to add new coloumn of colors
data <- data %>%
  mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                        yes = "Nasal", 
                        no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                    yes = "Dust", 
                                    no = "none")))

# Color corresponds to fold change directionality
colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
  theme_bw(base_size = 16) + # clean up theme
  ggtitle(label = "Volcano Plot", subtitle = "All OTUs Sampled") +  # add title
  xlab(expression(log[2]("Nasal" / "Dust"))) + # x-axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
  geom_vline(xintercept = 0, colour = "black") + # add line at 0
  geom_hline(yintercept = 1.3, colour = "black") + # p(0.05) = 1.3
  scale_color_manual(values = c("Nasal" = "#E64B35", 
                                "Dust" = "#3182bd", 
                                "none" = "#636363"), # change colors
                     name=NULL,
                     breaks=c("Nasal", "Dust"),
                     labels=c("Nasal", "Dust"))

# Plot figure
colored

colored + scale_y_continuous(trans = "log1p")

# Subset table to only show certain gene labels
top_labelled <- top_n(data, n = 4, wt = lfc)

# Add layer of text annotation to volcano plot.
colored + geom_text_repel(data = top_labelled, 
                          mapping = aes(label = gene), 
                          size = 3,
                          fontface = 'bold', 
                          color = 'black',
                          box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.5, "lines")) + 
  scale_y_continuous(trans = "log1p") 

volcano_OTU=data.frame(OTU_ID=results@rownames[which(results$padj<=0.05)],
         log2FoldChange=results$log2FoldChange[which(results$padj<=0.05)],
         padj=results$padj[which(results$padj<=0.05)])
volcano_OTU$prevalence="Dust"
volcano_OTU$prevalence[which(volcano_OTU$padj<=0.05 & volcano_OTU$log2FoldChange>0)]="Nasal"

saveRDS(volcano_OTU,"volcano_new_OTU.rds")



