################pca#################3



###############################################################
library(e1071)   
library(mvtnorm)
newdata<-readRDS("master_reduced.rds")
View(newdata)

dim(newdata)


newx<-newdata[,16:2920]
newy<-newdata[,4]


library(pls)
library(MASS)
ir.pca <- prcomp(sqrt(newx))


biplot(ir.pca,scale=0)

summary(ir.pca)
b<-ir.pca$rotation[,1:2]

#compute standard deviation of each principal component
std_dev <-ir.pca$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var[1:10]


#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]





plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
head()

