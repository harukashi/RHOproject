#######################logistic regression and plot###############################
Y <- master[order(master$StudyID, master$data_source),]
Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust.Faith <- Y_dust$Faith
Y_dust.Pielou <-Y_dust$Pielou
Y_dust.Chao1 <-Y_dust$Chao1
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal.Faith <- Y_nasal$Faith
Y_nasal.Pielou <- Y_nasal$Pielou
Y_nasal.Chao1 <- Y_nasal$Chao1
Y_alorwh <- Y[which(Y$data_source=="Dust"),]
Y_alorwh <- ifelse(Y_alorwh$wheeze_atopy_grp=="neither",0,1)
Y_size <- Y[which(Y$data_source=="Dust"),]
Y_size <- Y_size$site  
summary(glm(Y_alorwh ~ Y_dust.Faith + Y_nasal.Faith + + Y_size, family="binomial"))

summary(glm(Y_alorwh ~ Y_dust.Pielou + Y_nasal.Pielou + Y_size, family="binomial"))

summary(glm(Y_alorwh ~ Y_dust.Chao1 + Y_nasal.Chao1 + Y_size, family="binomial"))

summary(glm(Y_alorwh ~ Y_dust.Faith + Y_dust.Pielou + Y_dust.Chao1 + Y_nasal.Faith + Y_nasal.Pielou + Y_nasal.Chao1 + Y_size, family="binomial"))

summary(glm(Y_alorwh ~ Y_dust.Faith + Y_dust.Pielou + Y_dust.Chao1 + Y_nasal.Faith + Y_nasal.Pielou + Y_nasal.Chao1, family="binomial"))

sum.glm.fit=summary(glm.fit)
head(glm.fit$coefficients)

names.vec=as.vector(names(glm.fit$coefficients))
names.vec
names(est.vec)=NA
est.vec=sum.glm.fit$coefficients[,1]
se.glm=sum.glm.fit$coefficients[,2]

#significant non-time variables
variables=names.vec[c(2:7)]
mean=as.numeric(est.vec[c(2:7)])
lower=as.numeric(est.vec[c(2:7)]-1.96*se.glm[c(2:7)])
upper=as.numeric(est.vec[c(2:7)]+1.96*se.glm[c(2:7)])

#bimonthly variables
variables.time=names.vec[c(2:7)]
mean.time=as.numeric(est.vec[c(2:7)])
lower.time=as.numeric(est.vec[c(2:7)]-1.96*se.glm[c(2:7)])
upper.time=as.numeric(est.vec[c(2:7)]+1.96*se.glm[c(2:7)])

#Ordering the variables
variables.time=as.character(variables.time)
variables.time=factor(variables.time,levels=variables.time)

d=data.frame(variables.time, mean.time, lower.time, upper.time)
range=d$upper.time-d$lower.time
d=cbind(d, range)
d<-d[order(d$range), ]
d<-d[,1:4]
rownames(d)<-NULL
library(ggplot2)
ggplot(d, aes(x = variables.time, y = mean.time, col=factor(variables.time))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = upper.time, ymin = lower.time), position=position_dodge(width=0.9))+
  coord_flip()+
  geom_hline(yintercept=0)+
  ggtitle("Confidence Intervals of Regression Coefficients")+
  ylab("Regression Coefficient")+
  xlab("Independent Variables")+
  theme(axis.text=element_text(size=16))

#######################DAG Model##############################
suppressMessages(library(MASS))
suppressMessages(library(stats))
suppressMessages(library(geepack))
suppressMessages(library(lme4))
suppressMessages(library(car))
suppressMessages(library(boot))

A <- matrix(rep(0, 3*3), nrow=3)
D <-diag(3)
Y <- master[order(master$StudyID, master$data_source),]
Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Faith
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Faith
Y_alorwh <- Y[which(Y$data_source=="Dust"),]
Y_alorwh <- ifelse(Y_alorwh$wheeze_atopy_grp=="neither",0,1)

D[1, 1] <- var(Y_dust)
LM1 <- lm(Y_nasal ~ Y_dust)
A[2, 1] <- coef(LM1)[2]
D[2, 2] <- sigma(LM1)^2
LM2 <- glm(Y_alorwh ~ Y_dust + Y_nasal, family = "binomial")
A[3, 1] <- coef(LM2)[2]
A[3, 2] <- coef(LM2)[3]
D[3, 3] <- sigma(LM2)^2

Lambda <- solve(diag(3)-A) %*% D %*% t(solve(diag(3)-A))
corm <- round(cov2cor(Lambda), 3)
colnames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
rownames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
corm   #correlation matrix fitted from DAG using Faith

Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Pielou
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Pielou

D[1, 1] <- var(Y_dust)
LM1 <- lm(Y_nasal ~ Y_dust)
A[2, 1] <- coef(LM1)[2]
D[2, 2] <- sigma(LM1)^2
LM2 <- glm(Y_alorwh ~ Y_dust + Y_nasal, family = "binomial")
A[3, 1] <- coef(LM2)[2]
A[3, 2] <- coef(LM2)[3]
D[3, 3] <- sigma(LM2)^2

Lambda <- solve(diag(3)-A) %*% D %*% t(solve(diag(3)-A))
corm <- round(cov2cor(Lambda), 3)
colnames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
rownames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
corm  #correlation matrix fitted from DAG using Pielou

Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Chao1
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Chao1

D[1, 1] <- var(Y_dust)
LM1 <- lm(Y_nasal ~ Y_dust)
A[2, 1] <- coef(LM1)[2]
D[2, 2] <- sigma(LM1)^2
LM2 <- glm(Y_alorwh ~ Y_dust + Y_nasal, family = "binomial")
A[3, 1] <- coef(LM2)[2]
A[3, 2] <- coef(LM2)[3]
D[3, 3] <- sigma(LM2)^2

Lambda <- solve(diag(3)-A) %*% D %*% t(solve(diag(3)-A))
corm <- round(cov2cor(Lambda), 3)
colnames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
rownames(corm) <- c("dust", "nasal", "allergy and/or wheeze")
corm #correlation matrix fitted from DAG using Chao1

#diagnostic plots for DAG
hist(Y_dust)
qqPlot(LM1, main="QQ Plot Faith_dust to Faith_nasal")


######################Mediation Analysis############################
suppressMessages(library(mediation))
Y <- master[order(master$StudyID, master$data_source),]
Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Faith
dustmean<-mean(Y_dust)
Y_dust <- ifelse(Y_dust>dustmean, 1, 0)
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Faith
Y_alorwh <- Y[which(Y$data_source=="Dust"),]
Y_alorwh <- ifelse(Y_alorwh$wheeze_atopy_grp=="neither",0,1)
med.fit <- lm(Y_nasal ~ Y_dust)
out.fit <- glm(Y_alorwh ~ Y_nasal + Y_dust, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "Y_dust", mediator = "Y_nasal", robustSE = TRUE, sims = 100)
summary(med.out)  #Mediation analysis for Faith


Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Pielou
dustmean<-mean(Y_dust)
Y_dust <- ifelse(Y_dust>dustmean, 1, 0)
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Pielou
med.fit <- lm(Y_nasal ~ Y_dust)
out.fit <- glm(Y_alorwh ~ Y_nasal + Y_dust, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "Y_dust", mediator = "Y_nasal", robustSE = TRUE, sims = 100)
summary(med.out)  #Mediation analysis for Pielou


Y_dust <- Y[which(Y$data_source=="Dust"),]
Y_dust <- Y_dust$Chao1
dustmean<-mean(Y_dust)
Y_dust <- ifelse(Y_dust>dustmean, 1, 0)
Y_nasal <- Y[which(Y$data_source=="Nasal"),]
Y_nasal <- Y_nasal$Chao1
med.fit <- lm(Y_nasal ~ Y_dust)
out.fit <- glm(Y_alorwh ~ Y_nasal + Y_dust, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "Y_dust", mediator = "Y_nasal", robustSE = TRUE, sims = 100)
summary(med.out)  #Mediation analysis for Chao1








