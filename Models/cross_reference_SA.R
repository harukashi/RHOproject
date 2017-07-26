tbl.cr=sort(table(c(v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% volcano_OTU$OTU_ID],
v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.cor.OTU],
v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.dfPCA1.OTU],
v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.dfPCA2.OTU],

v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% volcano_OTU$OTU_ID],
v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.cor.OTU],
v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.dfPCA1.OTU],
v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.dfPCA2.OTU],

v.df8.OTU[v.df8.OTU %in% volcano_OTU$OTU_ID],
v.df8.OTU[v.df8.OTU %in% v.cor.OTU],
v.df8.OTU[v.df8.OTU %in% v.dfPCA1.OTU],
v.df8.OTU[v.df8.OTU %in% v.dfPCA2.OTU],

v.df9.OTU[v.df9.OTU %in% volcano_OTU$OTU_ID],
v.df9.OTU[v.df9.OTU %in% v.cor.OTU],
v.df9.OTU[v.df9.OTU %in% v.dfPCA1.OTU],
v.df9.OTU[v.df9.OTU %in% v.dfPCA2.OTU],

v.dvs.OTU[v.dvs.OTU %in% volcano_OTU$OTU_ID],
v.dvs.OTU[v.dvs.OTU %in% v.cor.OTU],
v.dvs.OTU[v.dvs.OTU %in% v.dfPCA1.OTU],
v.dvs.OTU[v.dvs.OTU %in% v.dfPCA2.OTU],

v.nvs.OTU[v.nvs.OTU %in% volcano_OTU$OTU_ID],
v.nvs.OTU[v.nvs.OTU %in% v.cor.OTU],
v.nvs.OTU[v.nvs.OTU %in% v.dfPCA1.OTU],
v.nvs.OTU[v.nvs.OTU %in% v.dfPCA2.OTU],

v.dsupvol.OTU[v.dsupvol.OTU %in% volcano_OTU$OTU_ID],
v.dsupvol.OTU[v.dsupvol.OTU %in% v.cor.OTU],
v.dsupvol.OTU[v.dsupvol.OTU %in% v.dfPCA1.OTU],
v.dsupvol.OTU[v.dsupvol.OTU %in% v.dfPCA2.OTU],

v.nsupvol.OTU[v.nsupvol.OTU %in% volcano_OTU$OTU_ID],
v.nsupvol.OTU[v.nsupvol.OTU %in% v.cor.OTU],
v.nsupvol.OTU[v.nsupvol.OTU %in% v.dfPCA1.OTU],
v.nsupvol.OTU[v.nsupvol.OTU %in% v.dfPCA2.OTU])))

sort(table(c(v.cor.OTU[v.cor.OTU %in% volcano_OTU$OTU_ID],
             v.dfPCA1.OTU[v.dfPCA1.OTU %in% volcano_OTU$OTU_ID],
             v.cor.OTU[v.cor.OTU %in% v.dfPCA1.OTU],
             v.dfPCA1.OTU[v.dfPCA1.OTU %in% v.dfPCA2.OTU],
             v.dfPCA2.OTU[v.dfPCA2.OTU %in% volcano_OTU$OTU_ID],
             v.cor.OTU[v.cor.OTU %in% v.dfPCA2.OTU])))
             
sort(table(c(v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.rfn.OTU$Feature],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.df8.OTU],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.df9.OTU],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.dvs.OTU],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.nvs.OTU],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.dsupvol.OTU],
             v.rfd.OTU$Feature[v.rfd.OTU$Feature %in% v.nsupvol.OTU],
             
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.df8.OTU],
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.df9.OTU],
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.dvs.OTU],
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.nvs.OTU],
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.dsupvol.OTU],
             v.rfn.OTU$Feature[v.rfn.OTU$Feature %in% v.nsupvol.OTU],
             
             v.df8.OTU[v.df8.OTU %in% v.df9.OTU],
             v.df8.OTU[v.df8.OTU %in% v.dvs.OTU],
             v.df8.OTU[v.df8.OTU %in% v.nvs.OTU],
             v.df8.OTU[v.df8.OTU %in% v.dsupvol.OTU],
             v.df8.OTU[v.df8.OTU %in% v.nsupvol.OTU],

             v.df9.OTU[v.df9.OTU %in% v.dvs.OTU],
             v.df9.OTU[v.df9.OTU %in% v.nvs.OTU],
             v.df9.OTU[v.df9.OTU %in% v.dsupvol.OTU],
             v.df9.OTU[v.df9.OTU %in% v.nsupvol.OTU],
             
             v.dvs.OTU[v.dvs.OTU %in% v.nvs.OTU],
             v.dvs.OTU[v.dvs.OTU %in% v.dsupvol.OTU],
             v.dvs.OTU[v.dvs.OTU %in% v.nsupvol.OTU],
             
             v.nvs.OTU[v.nvs.OTU %in% v.dsupvol.OTU],
             v.nvs.OTU[v.nvs.OTU %in% v.nsupvol.OTU],
             
             v.dsupvol.OTU[v.dsupvol.OTU %in% v.nsupvol.OTU]
             )))

v.sup.learn=unique(c(v.rfd.OTU$Feature,v.rfn.OTU$Feature,v.df8.OTU,v.df9.OTU,v.dvs.OTU,v.nvs.OTU,v.dsupvol.OTU,v.nsupvol.OTU))
v.sup.learn.d=unique(c(v.rfd.OTU$Feature,v.df8.OTU,v.dvs.OTU,v.dsupvol.OTU))
v.sup.learn.n=unique(c(v.rfn.OTU$Feature,v.df9.OTU,v.nvs.OTU,v.nsupvol.OTU))

names.cr=names(tbl.cr)
master_cr=cbind(master[,1:12],master[,c("allergic_wheeze",names.cr)])

df.names.cr=data.frame(names.cr)
saveRDS(df.names.cr,"names_cr.rds")



df_sup_learn_dust_phylum=as.data.frame(table(taxonomy$phylum[taxonomy$OTU_ID %in% v.sup.learn.d]))
df_sup_learn_dust_class=as.data.frame(table(taxonomy$class[taxonomy$OTU_ID %in% v.sup.learn.d]))

df_sup_learn_nasal_phylum=as.data.frame(table(taxonomy$phylum[taxonomy$OTU_ID %in% v.sup.learn.n]))
df_sup_learn_nasal_class=as.data.frame(table(taxonomy$class[taxonomy$OTU_ID %in% v.sup.learn.n]))

df_all_phylum=as.data.frame(table(taxonomy$phylum))
df_all_class=as.data.frame(table(taxonomy$class))

master_dust=master[master$data_source=="Dust",]
keep.col=colnames(master_dust[,which(colSums(master_dust[,16:10133])>0)+15])
df_dust_phylum=as.data.frame(table(taxonomy$phylum[taxonomy$OTU_ID %in% keep.col]))
df_dust_class=as.data.frame(table(taxonomy$class[taxonomy$OTU_ID %in% keep.col]))

master_nasal=master[master$data_source=="Nasal",]
keep.col=colnames(master_nasal[,which(colSums(master_nasal[,16:10133])>0)+15])
df_nasal_phylum=as.data.frame(table(taxonomy$phylum[taxonomy$OTU_ID %in% keep.col]))
df_nasal_class=as.data.frame(table(taxonomy$class[taxonomy$OTU_ID %in% keep.col]))

saveRDS(df_dust_phylum,"df_dust_phylum.rds")
saveRDS(df_dust_class,"df_dust_class.rds")

saveRDS(df_nasal_phylum,"df_nasal_phylum.rds")
saveRDS(df_nasal_class,"df_nasal_class.rds")

saveRDS(df_sup_learn_dust_phylum,"df_sup_learn_dust_phylum.rds")
saveRDS(df_sup_learn_dust_class,"df_sup_learn_dust_class.rds")

saveRDS(df_sup_learn_nasal_phylum,"df_sup_learn_nasal_phylum.rds")
saveRDS(df_sup_learn_nasal_class,"df_sup_learn_nasal_class.rds")

saveRDS(df_all_phylum,"df_all_phylum.rds")
saveRDS(df_all_class,"df_all_class.rds")


