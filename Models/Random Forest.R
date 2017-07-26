#Random Forest
library(xgboost)
Y <- master[order(master$StudyID, master$data_source),]
Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_alorwh <- Y[which(Y$data_source=="Dust"),]
Y_alorwh <- ifelse(Y_alorwh$wheeze_atopy_grp=="neither",0,1)

xgdatax<-Y_dust[, 16:10136]
xgdatax<-as.matrix(xgdatax)
set.seed(1)
bst <- xgboost(data = xgdatax, label = Y_alorwh, max_depth = 20, nthread = 2, nrounds = 1, objective = "binary:logistic",                        num_parallel_tree = 500, subsample = 0.6, colsample_bytree = 0.1)
importance_matrix <- xgb.importance(colnames(xgdatax), model = bst)
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance", top_n = 10, main="Dust Sample")

xgdatax<-Y_nasal[, 16:10136]
xgdatax<-as.matrix(xgdatax)
set.seed(1)
bst <- xgboost(data = xgdatax, label = Y_alorwh, max_depth = 20, nthread = 2, nrounds = 1, objective = "binary:logistic",                        num_parallel_tree = 500, subsample = 0.6, colsample_bytree = 0.1)
importance_matrix <- xgb.importance(colnames(xgdatax), model = bst)
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance", top_n = 10, main="Nasal Sample")
