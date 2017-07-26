## This part of R.code is to do Elastic Net for Dust and Nasal.
## dataset:new_reduced_master
## Everything with subscript 6 is from dust
## Everything with subscript 5 is from nasal
suppressMessages(library(glmnet))
load("C:/Users/amberzhang/Desktop/SAMSI_2017/RHO/new_reduced.RData")
new_dust_data <- new_reduced[which(new_reduced$data_source=="Dust"),]
a <- which(new_dust_data$wheeze_atopy_grp=="neither")
new_dust_data$wheeze_atopy_grp[a]=rep(0,30)
b <-which(new_dust_data$wheeze_atopy_grp!=0)
new_dust_data$wheeze_atopy_grp[b]=rep(1,44)
new_dust_dm <- as.matrix(new_dust_data[,16:9431])

set.seed(111)
model_6 <- cv.glmnet (new_dust_dm,new_dust_data$wheeze_atopy_grp,
                      family = "binomial",
                      type.measure = "class",intercept=FALSE,
                      alpha = .5,nfolds = 5,standardize=TRUE)
model_6.fit<- model_6$glmnet.fit
model_62 <- cv.glmnet (new_dust_dm,new_dust_data$wheeze_atopy_grp,
                       family = "binomial",
                       type.measure = "class",intercept=FALSE,
                       alpha = 1,nfolds = 5,standardize=TRUE) # LASSO
par(mfrow=c(2,2))
plot(model_6)
plot(model_62)
plot(model_6$glmnet.fit,"dev")
plot(model_62$glmnet.fit,"dev")
mean(model_6$cvsd)
mean(model_62$cvsd)
min(model_6$cvm)
min(model_62$cvm)
par(mfrow=c(1,2))
plot(model_6$lambda,model_6$nzero,xlab = "lambda",ylab = "# of nonzero variables"
)
plot(model_62$lambda,model_62$nzero,xlab = "lambda",ylab = "# of nonzero variables",
     xlim =c(0,.4),ylim = c(0,200))

par(mfrow=c(1,1))
plot(log(model_6$lambda),model_6$cvm,pch=20,col="red",
     xlab="log(Lambda)",ylab=model_6$name,main = "alpha=.5--red alpha=1--blue")
points(log(model_62$lambda),model_62$cvm,pch=20,col="blue")
c_multi_6 <- as.matrix(coef(model_6,s="lambda.min"))
nz_c_6 <- cbind(rownames(c_multi_6)[which(c_multi_6!=0)],
                c_multi_6[which(c_multi_6!=0)])
nz_c_6 <-  nz_c_6[order(nz_c_6[,2],decreasing = TRUE),]

new_dust_valsel <- nz_c_6[c(1:5,36:39),]


new_nasal_data <- new_reduced[which(new_reduced$data_source=="Nasal"),]
a <- which(new_nasal_data$wheeze_atopy_grp=="neither")
new_nasal_data$wheeze_atopy_grp[a]=rep(0,30)
b <-which(new_nasal_data$wheeze_atopy_grp!=0)
new_nasal_data$wheeze_atopy_grp[b]=rep(1,44)
new_nasal_dm <- as.matrix(new_nasal_data[,16:9431])

set.seed(311)
model_5 <- cv.glmnet (new_nasal_dm,new_nasal_data$wheeze_atopy_grp,
                      family = "binomial",
                      type.measure = "class",intercept=FALSE,
                      alpha =0.5,nfolds = 5,standardize=TRUE)
model_52 <- cv.glmnet (new_nasal_dm,new_nasal_data$wheeze_atopy_grp,
                       family = "binomial",
                       type.measure = "class",intercept=FALSE,
                       alpha =1,nfolds = 5,standardize=TRUE)
par(mfrow=c(2,2))
plot(model_5)
plot(model_52)
plot(model_5$glmnet.fit,"dev")
plot(model_52$glmnet.fit,"dev")

mean(model_5$cvsd)
mean(model_52$cvsd)
min(model_5$cvm)
min(model_52$cvm)

par(mfrow=c(1,2))
plot(model_5$lambda,model_5$nzero,xlab = "lambda",ylab = "# of nonzero variables"
)
plot(model_52$lambda,model_52$nzero,xlab = "lambda",ylab = "# of nonzero variables",
     xlim = c(0,0.4),ylim = c(0,200)     )

suppressMessages(library(knitr))
c_multi_5 <- as.matrix(coef(model_5,s="lambda.min"))
nz_c_5 <- cbind(rownames(c_multi_5)[which(c_multi_5!=0)],
                c_multi_5[which(c_multi_5!=0)])
c_multi_52 <- as.matrix(coef(model_52,s="lambda.min"))
nz_c_52 <- cbind(rownames(c_multi_52)[which(c_multi_52!=0)],
                 c_multi_52[which(c_multi_52!=0)])
nz_c_5 <- as.data.frame( nz_c_5[order(nz_c_5[,2],decreasing = TRUE),])
nz_c_52 <- as.data.frame( nz_c_52[order(nz_c_52[,2],decreasing = TRUE),])
df <- merge(nz_c_5,nz_c_52,by = "V1",all=TRUE)
colnames(df) <- c("OTU","Elastic Net","LASSO")
df_1 <- df[1:40,]
df_2 <- df[41:83,]
#new_nasal_valsel <- nz_c_5[c(1:9,39:61),]
#kable(new_nasal_valsel)

df_not_lasso <- df[is.na(df$LASSO=="NA"),]
df_in_lasso <- df[which(df$LASSO!="NA"),]
compare_df <- df[which(df$LASSO!="NA"),]
compare_df_1 <- compare_df[,c(1,3)]
compare_df_2 <- compare_df[,c(1,2)]
names <-  as.character( df_not_lasso$OTU)
names_1 <- as.character(df_in_lasso$OTU)
not_lasso_otu <- new_nasal_data[,names]
in_lasso_otu <- new_nasal_data[,names_1]
dim(not_lasso_otu)
dim(in_lasso_otu)
corr_mat <- cor(not_lasso_otu,in_lasso_otu)
corr_mat[ corr_mat <0.5 & corr_mat > -.5]<- 0
index=as.matrix(which(corr_mat!=0,arr.ind = T))
length(unique(index[,1]))

