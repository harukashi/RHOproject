vignette("eulerr")
set.seed(123)
df_venn=data.frame(OTU_ID=names.cr)
df_venn$volcano=df_venn$OTU_ID %in% volcano_OTU$OTU_ID
df_venn$cor=df_venn$OTU_ID %in% v.cor.OTU
df_venn$pcad=df_venn$OTU_ID %in% v.dfPCA1.OTU
df_venn$pcan=df_venn$OTU_ID %in% v.dfPCA2.OTU
df_venn$supvold=df_venn$OTU_ID %in% v.dsupvol.OTU
df_venn$supvoln=df_venn$OTU_ID %in% v.nsupvol.OTU
df_venn$rfd=df_venn$OTU_ID %in% v.rfd.OTU$Feature
df_venn$rfn=df_venn$OTU_ID %in% v.rfn.OTU$Feature
df_venn$svmd=df_venn$OTU_ID %in% v.df8.OTU
df_venn$svmn=df_venn$OTU_ID %in% v.df9.OTU
df_venn$vsd=df_venn$OTU_ID %in% v.dvs.OTU
df_venn$vsn=df_venn$OTU_ID %in% v.nvs.OTU
fit1 = euler(df_venn[,2:13])
par(mar = c(0, 0, 0, 0))
plot(fit1)



names.cr.d=unique(c(v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% volcano_OTU$OTU_ID],
                    v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.cor.OTU],
                    v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.dfPCA1.OTU],
                    
                    v.df8.OTU[v.df8.OTU %in% volcano_OTU$OTU_ID],
                    v.df8.OTU[v.df8.OTU %in% v.cor.OTU],
                    v.df8.OTU[v.df8.OTU %in% v.dfPCA1.OTU],
                    
                    v.dvs.OTU[v.dvs.OTU %in% volcano_OTU$OTU_ID],
                    v.dvs.OTU[v.dvs.OTU %in% v.cor.OTU],
                    v.dvs.OTU[v.dvs.OTU %in% v.dfPCA1.OTU],
                    
                    v.dsupvol.OTU[v.dsupvol.OTU %in% volcano_OTU$OTU_ID],
                    v.dsupvol.OTU[v.dsupvol.OTU %in% v.cor.OTU],
                    v.dsupvol.OTU[v.dsupvol.OTU %in% v.dfPCA1.OTU]))



set.seed(123)
df_venn_d=data.frame(OTU_ID=names.cr)
df_venn_d$Volcano=df_venn_d$OTU_ID %in% volcano_OTU$OTU_ID
df_venn_d$Correlation=df_venn_d$OTU_ID %in% v.cor.OTU
df_venn_d$PCA=df_venn_d$OTU_ID %in% v.dfPCA1.OTU
df_venn_d$Super_Volcano=df_venn_d$OTU_ID %in% v.dsupvol.OTU
df_venn_d$Random_Forest=df_venn_d$OTU_ID %in% v.rfd.OTU$Feature
df_venn_d$SVM=df_venn_d$OTU_ID %in% v.df8.OTU
df_venn_d$Elastic_Net=df_venn_d$OTU_ID %in% v.dvs.OTU
fit2 = euler(df_venn_d[,2:8])
par(mar = c(0, 0, 0, 0))
plot(fit2, counts = TRUE, labels = c("Volcano", "Correlation", "Sparse PCA", "Super Volcano", "Random Forest", "Sparse SVM", "Elastic Net"), main="Important OTUs from Dust Samples")



names.cr.n=unique(c(v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% volcano_OTU$OTU_ID],
                    v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.cor.OTU],
                    v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.dfPCA2.OTU],
                    
                    v.df9.OTU[v.df9.OTU %in% volcano_OTU$OTU_ID],
                    v.df9.OTU[v.df9.OTU %in% v.cor.OTU],
                    v.df9.OTU[v.df9.OTU %in% v.dfPCA2.OTU],
                    
                    v.nvs.OTU[v.nvs.OTU %in% volcano_OTU$OTU_ID],
                    v.nvs.OTU[v.nvs.OTU %in% v.cor.OTU],
                    v.nvs.OTU[v.nvs.OTU %in% v.dfPCA2.OTU],
                    
                    v.nsupvol.OTU[v.nsupvol.OTU %in% volcano_OTU$OTU_ID],
                    v.nsupvol.OTU[v.nsupvol.OTU %in% v.cor.OTU],
                    v.nsupvol.OTU[v.nsupvol.OTU %in% v.dfPCA2.OTU]))



set.seed(123)
df_venn_n=data.frame(OTU_ID=names.cr)
df_venn_n$Volcano=df_venn_n$OTU_ID %in% volcano_OTU$OTU_ID
df_venn_n$Correlation=df_venn_n$OTU_ID %in% v.cor.OTU
df_venn_n$PCA=df_venn_n$OTU_ID %in% v.dfPCA2.OTU
df_venn_n$Super_Volcano=df_venn_n$OTU_ID %in% v.nsupvol.OTU
df_venn_n$Random_Forest=df_venn_n$OTU_ID %in% v.rfn.OTU$Feature
df_venn_n$SVM=df_venn_n$OTU_ID %in% v.df9.OTU
df_venn_n$Elastic_Net=df_venn_n$OTU_ID %in% v.nvs.OTU
fit3 = euler(df_venn_n[,2:8])
par(mar = c(0, 0, 0, 0))
plot(fit3, counts = TRUE, labels = c("Volcano", "Correlation", "Sparse PCA", "Super Volcano", "Random Forest", "Sparse SVM", "Elastic Net"), main="Important OTUs from Nasal Samples")


# x <- euler(c("A" = 10, "B" = 5, "A&B" = 3))
# plot(x, labels = c("foo", "bar"), fill_alpha = 0.7)
