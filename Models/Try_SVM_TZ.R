



library(dplyr)
library(e1071)   
library(mvtnorm)
library(elasticnet)
newdata<-readRDS("master_reduced.rds")



newx<-newdata[,16:2920]
newy<-newdata[,4]

out = svm(newx, newy, kernel="linear", type="C", cost = 1,data=newdata2)

psi=t(out$coefs)%*%as.matrix(newx[out$index,])
pripsi = eigen(t(psi)%*%psi)$vectors[,1:2]  
head((mat.sort(pripsi,1,decreasing=FALSE)))
View(pripsi)

head(pripsi)
pripsi[order(pripsi[,1],pripsi[,2],decreasing=TRUE),]

beta = var(newx)%*%pripsi

View(data.frame(beta))
a<-mat.sort(beta,1,decreasing=TRUE)



b<-mat.sort(beta,2,decreasing=TRUE)

saveRDS(a,"a.rds")

saveRDS(b,"b.rds")


















######################weezy##########




newy2<-newdata[,5]

out2 = svm(newx, newy2, kernel="linear", type="C", cost = 1)

psi2=t(out2$coefs)%*%as.matrix(newx[out2$index,])
pripsi2 = eigen(t(psi2)%*%psi2)$vectors[,1:2]  




beta2 = var(newx)%*%pripsi2


c<-mat.sort(beta2,1,decreasing=TRUE)



d<-mat.sort(beta2,2,decreasing=TRUE)

saveRDS(c,"c.rds")

saveRDS(d,"d.rds")

View(c)










