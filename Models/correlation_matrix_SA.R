temp=merge(sample_data_stacked[,1:3],rbind(nasal_OTU,dust_OTU),"sampleid")
rm.col=unique(c(which(colSums(temp[temp$data_source=="Nasal",4:10121])==0),
                which(colSums(temp[temp$data_source=="Dust",4:10121])==0)))+3
temp=temp[,-rm.col]

temper=merge(temp[temp$data_source=="Nasal",],temp[temp$data_source=="Dust",],"StudyID")
which(colnames(temper)=="sampleid.y")
which(colnames(temper)=="data_source.y")
tempest=cbind(temper[,1:3],temper[,1056:1057],temper[,4:1055],temper[,1058:2109])
head(colnames(tempest[,6:1057]))
head(colnames(tempest[,1058:2109]))

mat.cor=matrix(NA,1052,4)
for(i in 6:1057) {
  j=i+1052
  mat.cor[i-5,1] = cor.test(tempest[,i], tempest[,i+1052])$estimate
  mat.cor[i-5,2] = cor.test(tempest[,i], tempest[,i+1052])$p.value
  mat.cor[i-5,3] = cor.test(tempest[,i], tempest[,i+1052],method="spearman",exact=F)$estimate
  mat.cor[i-5,4] = cor.test(tempest[,i], tempest[,i+1052],method="spearman",exact=F)$p.value
}

cor_OTU=as.data.frame(mat.cor)
cor_OTU=data.frame(OTU_ID=colnames(temp)[4:1055],cor_OTU)
cor_OTU$pearson.sig=F
cor_OTU$pearson.sig[cor_OTU$V2<=0.01]=T
cor_OTU$spearman.sig=F
cor_OTU$spearman.sig[cor_OTU$V4<=0.01]=T
cor_OTU$both.sig=F
cor_OTU$both.sig[cor_OTU$V2<=0.01 & cor_OTU$V4<=0.01]=T

v.cor.OTU=cor_OTU$OTU_ID[cor_OTU$spearman.sig==T]

saveRDS(cor_OTU,"cor_OTU.rds")

v.cor.OTU[v.cor.OTU %in% volcano_OTU$OTU_ID]
taxonomy[taxonomy$OTU_ID %in% v.cor.OTU[v.cor.OTU %in% volcano_OTU$OTU_ID],]

volcano_OTU[volcano_OTU$OTU_ID %in% v.cor.OTU,]