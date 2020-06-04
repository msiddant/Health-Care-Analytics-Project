########################################################################################
#Packages
########################################################################################
library(tidytext)
library(dplyr)
library(nnet)
library(neuralnet)
library(caret)
library(forecast)
library(gains)
library(GGally)
library(ggcorrplot)
library(reshape2)
library(forecast)
library(MASS)
library(randomForest)
library(naivebayes)
library(rpart)
library(caret)
library(e1071)
##########################################################################################

##HEALTH PROJECT _ LIVER

liver.df <- read.csv("LIVER_Final_2.csv",header = T,na.strings = c(""))

str(liver.df)
attach(liver.df)
#Creating a varaible
liver.df$gendermat<-ifelse(GENDER==GENDER_DON,1,0)

#Fixing the class for each variable 
liver.df$DIAB <- as.factor(liver.df$DIAB)
liver.df$ETHCAT <- as.factor(liver.df$ETHCAT)
liver.df$PSTATUS <- as.factor(liver.df$PSTATUS)
liver.df$REGION <- as.factor(liver.df$REGION)
liver.df$TX_PROCEDUR_TY <- as.factor(liver.df$TX_PROCEDUR_TY)
liver.df$MED_COND_TRR <- as.factor(liver.df$MED_COND_TRR)
liver.df$COD_CAD_DON <- as.factor(liver.df$COD_CAD_DON)
liver.df$ETHCAT_DON <- as.factor(liver.df$ETHCAT_DON)
liver.df$ABO_MAT <- as.factor(liver.df$ABO_MAT)
liver.df$TX_Year <- as.factor(liver.df$TX_Year)
liver.df$LISTYR <- as.factor(liver.df$LISTYR)
liver.df$GENDER <- as.factor(liver.df$GENDER)
liver.df$GENDER_DON <- as.factor(liver.df$GENDER_DON)
liver.df$gendermat<- as.factor(liver.df$gendermat)
str(liver.df)
#Only considering patients who survived the first year of surgery 
liver.df1 <- liver.df %>%
  filter(PTIME > 365) %>%
  mutate(PTIME=(PTIME-365))

Liver_Num<-liver.df1[,!unlist(lapply(liver.df1, is.factor))]
#########################################################################################################
#Correlation
##########################################################################################################
corr <- round(cor(Liver_Num,use="complete.obs"), 2)
ggcorrplot(corr)

###########################################################################################################
#LINEAR REGRESSION
###########################################################################################################
#DATA FILTERING AND NORMALIZATION
#Survival after one - three years of surgery
liver.df2<-liver.df1 %>%
  filter(PTIME<=1095)

#Survival after 1 - 5 years of surgery
liver.df3<-liver.df1 %>%
  filter(PTIME<=1825)
##Seperating Numerical and Categorical Varaibles for 
#liver.df2 ( 1-3 yrs)
Liver_Num1<-liver.df2[,!unlist(lapply(liver.df2, is.factor))]
Liver_Factor1<- liver.df2[,unlist(lapply(liver.df2, is.factor))]
#liver.df3 ( 1-5 yrs)
Liver_Num2<-liver.df3[,!unlist(lapply(liver.df3, is.factor))]
Liver_Factor2<- liver.df3[,unlist(lapply(liver.df3, is.factor))]


##Normalizong the the numerical varaibles and combining into a dataframe
#liver.df2 ( 1-3 yrs)
Liver_Num_Norm1<- scale(Liver_Num1)
Liver_Norm13 <- cbind(Liver_Num_Norm1,Liver_Factor1)
#liver.df3 ( 1-5 yrs)
Liver_Num_Norm2<- scale(Liver_Num2)
Liver_Norm15 <- cbind(Liver_Num_Norm2,Liver_Factor2)


#(1-3 yrs)
#Training and test data sets
set.seed(1)
training13 <- sample(1:nrow(Liver_Norm13), dim(Liver_Norm13)[1]*.6)
trainData13 <- Liver_Norm13[training13,]
testData13 <- Liver_Norm13[-training13,]

testData13<- testData13 %>%
  filter(TX_Year !=2002)
testData15<- testData15 %>%
  filter(TX_Year !=2002)

trainData13<- trainData13 %>%
  filter(TX_Year !=2002)
trainData15<- trainData15 %>%
  filter(TX_Year !=2002)

#(1-5 yrs)
#Training and test data sets
set.seed(1)
training15 <- sample(1:nrow(Liver_Norm15), dim(Liver_Norm15)[1]*.6)
trainData15 <- Liver_Norm15[training15,]
testData15 <- Liver_Norm15[-training15,]

##Predicting PTime for 1-3 yrs
lm.fit1 <- lm(PTIME~FINAL_MELD_SCORE+REGION+LiverSize+LiverSizeDon
              +ALCOHOL_HEAVY_DON+MALIG+TX_Year,data = trainData13)

summary(lm.fit1)
par(mfrow = c(2,2))
plot(lm.fit1)
#Identifying Outliers
#identify(trainData13)

#Prediction Accuracy
#training
lm.pred13 = predict(lm.fit1, newdata = trainData13)
accuracy(lm.pred13, trainData13$PTIME)
#test
lm.pred13 = predict(lm.fit1, newdata = testData13)
accuracy(lm.pred13, testData13$PTIME)

##Predicting PTime for 1-5 yrs
lm.fit2 <- lm(PTIME~FINAL_MELD_SCORE+REGION
              +LiverSize+LiverSizeDon
              +ALCOHOL_HEAVY_DON
              +MALIG+TX_Year,data = trainData15)
summary(lm.fit2)
par(mfrow = c(2,2))
plot(lm.fit2)

#Identifying Outliers
#identify(trainData15)

#Prediction Accuracy
#training
lm.pred15 = predict(lm.fit2, newdata = trainData15)
accuracy(lm.pred15, trainData15$PTIME)
#test
lm.pred15 = predict(lm.fit2, newdata = testData15)
accuracy(lm.pred15, testData15$PTIME)

#Distribution of the validation error
par(mfrow = c(1,1))
all.residuals13 <- testData13$PTIME-lm.pred13
hist(all.residuals13, xlab = "Residuals", main = "Years 1 to 3")


all.residuals15 <- testData15$PTIME-lm.pred15
hist(all.residuals15, xlab = "Residuals", main = "Years 1 to 5")


#On unscalled
#1-3
samplevaraibles <- c("PTIME","FINAL_MELD_SCORE","REGION",
                     "LiverSize","LiverSizeDon",
                     "ALCOHOL_HEAVY_DON",
                     "MALIG","TX_Year")

UnscalledLiver13 <- na.omit(liver.df2[,samplevaraibles])
UnscalledLiver15 <- na.omit(liver.df3[,samplevaraibles])

unscalledTest13<- UnscalledLiver13[-training13,] %>%
  filter(TX_Year!=2002)
unscalledTest15<- UnscalledLiver15[-training15,] %>%
  filter(TX_Year!=2002)

ht.pred = predict(lm.fit1, newdata = unscalledTest13)
View(ht.pred)
accuracy(ht.pred, unscalledTest13$PTIME)


#On unscalled
ht.pred2 = predict(lm.fit2, newdata = unscalledTest15)
View(ht.pred2)
accuracy(ht.pred2, unscalledTest15$PTIME)


##NN
Liver_Norm213 <- na.omit(Liver_Norm13[,samplevaraibles])
Liver_Norm215 <- na.omit(Liver_Norm15[,samplevaraibles])

#1-3
norm.var <-!unlist(lapply(Liver_Norm213, is.factor))
name.var<-names(Liver_Norm213[,!unlist(lapply(Liver_Norm213, is.factor))])
liver.norm.df213 = cbind(Liver_Norm213[,norm.var],
                         class.ind(Liver_Norm213$REGION),
                         class.ind(Liver_Norm213$ALCOHOL_HEAVY_DON),
                         class.ind(Liver_Norm213$MALIG),
                         class.ind(Liver_Norm213$TX_Year))

names(liver.norm.df213) = c(name.var,
                            paste("REGION",c(1:11),sep = ""),
                            paste("ALCOHOL_HEAVY_DON",c(1:3),sep = ""),
                            paste("MALIG",c(1:3),sep = ""),
                            paste("LISTYR",c(01:18),sep = ""))


training213 <- sample(1:nrow(liver.norm.df213), dim(liver.norm.df213)[1]*.60)
traindata213 <- liver.norm.df213[training213,]
Validdata213 <- liver.norm.df213[-training213,]

#2 hidden layer
nn213 <- neuralnet::neuralnet(PTIME ~., data = traindata213, linear.output = T,hidden = 2)
plot(nn213,col.entry.synapse = "red", col.entry = "brown",
     col.hidden = "green", col.hidden.synapse = "black",
     col.out = "yellow", col.out.synapse = "purple",
     col.intercept = "green", fontsize = 10,
     show.weights = TRUE ,rep="best")

#Prediction
#Training
train_pred213 <- compute(nn213, traindata213)
accuracy(unlist(train_pred213), traindata213$PTIME)
#Validation
Valid_pred213 <- compute(nn213, Validdata213)
accuracy(unlist(Valid_pred213), Validdata213$PTIME)

#1-5
norm.var <-!unlist(lapply(Liver_Norm215, is.factor))
name.var<-names(Liver_Norm215[,!unlist(lapply(Liver_Norm215, is.factor))])
liver.norm.df215 = cbind(Liver_Norm215[,norm.var],
                         class.ind(Liver_Norm215$REGION),
                         class.ind(Liver_Norm215$ALCOHOL_HEAVY_DON),
                         class.ind(Liver_Norm215$MALIG),
                         class.ind(Liver_Norm215$TX_Year))

names(liver.norm.df215) = c(name.var,
                            paste("REGION",c(1:11),sep = ""),
                            paste("ALCOHOL_HEAVY_DON",c(1:3),sep = ""),
                            paste("MALIG",c(1:3),sep = ""),
                            paste("LISTYR",c(01:18),sep = ""))



training215 <- sample(1:nrow(liver.norm.df215), dim(liver.norm.df215)[1]*.60)
traindata215 <- liver.norm.df215[training215,]
Validdata215 <- liver.norm.df215[-training215,]

#2 hidden layer
nn215 <- neuralnet::neuralnet(PTIME ~., data = traindata215, linear.output = T,hidden = 2)
plot(nn215,col.entry.synapse = "red", col.entry = "brown",
     col.hidden = "green", col.hidden.synapse = "black",
     col.out = "yellow", col.out.synapse = "purple",
     col.intercept = "green", fontsize = 10,
     show.weights = TRUE ,rep="best")

#prediction
#training
train_pred215 <- compute(nn215, traindata215)
accuracy(unlist(train_pred215), traindata215$PTIME)
#validation
Valid_pred215 <- compute(nn215, Validdata215)
accuracy(unlist(Valid_pred215), Validdata215$PTIME)
#####################################################################

#CLASSIFICATION MODELS FOR SURVIVAL AFTER ONE TO THREE YEARS AFTER SURGERY 
set.seed(1)
training <- sample(1:nrow(Liver_Norm13), dim(Liver_Norm13)[1]*.6)
trainData <- Liver_Norm13[training,]
testData <- Liver_Norm13[-training,]

###############
#LOGISTIC MODEL 
###############
log.mod = glm(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
summary(log.mod)
log.mod.pred = predict(log.mod,newdata = testData, type = "response")
log.mod.prob = ifelse(log.mod.pred >.5, "1", "0") 
table(log.mod.prob,testData$PSTATUS)
mean(na.omit(log.mod.prob == testData$PSTATUS)) #0.7143623
mean(na.omit(log.mod.prob != testData$PSTATUS)) #0.2860512

#same model with blood type 
log.mod = glm(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
summary(log.mod)
log.mod.pred = predict(log.mod,newdata = testData, type = "response")
log.mod.prob = ifelse(log.mod.pred >.5, "1", "0") 
table(log.mod.prob,testData$PSTATUS)
mean(na.omit(log.mod.prob == testData$PSTATUS)) #0.7145766 
mean(na.omit(log.mod.prob != testData$PSTATUS)) #0.2854234

###############
#LDA MODEL 
###############
lda.mod1 = lda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
lda.pred1=predict(lda.mod1, testData)
predictedClass=lda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) # 0.7151281
mean(na.omit(predictedClass!=testData$PSTATUS)) # 0.2848719

#same model but with blood type 
lda.mod1 = lda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
lda.pred1=predict(lda.mod1, testData)
predictedClass=lda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) # 0.7142695
mean(na.omit(predictedClass!=testData$PSTATUS)) # 0.2855305

###############
#QDA MODEL 
###############
qda.mod1 = qda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
qda.pred1=predict(qda.mod1, testData)
predictedClass=qda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) #0.7079314
mean(na.omit(predictedClass!=testData$PSTATUS)) #0.2920686

qda.mod1 = qda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
qda.pred1=predict(qda.mod1, testData)
predictedClass=qda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) #0.7056806
mean(na.omit(predictedClass!=testData$PSTATUS)) #0.2943194

######################
#NN FOR CLASSIFICATION  
######################
#####################################################
samplevaraibles <- c("FINAL_MELD_SCORE","REGION","HGT_CM_CALC"
                     ,"HGT_CM_DON_CALC"
                     ,"LiverSize","LiverSizeDon","DIAB"
                     ,"DIABETES_DON","MED_COND_TRR"
                     ,"gendermat","PSTATUS")

Liver_Norm213 <- na.omit(Liver_Norm13[,samplevaraibles])
Liver_Norm215 <- na.omit(Liver_Norm15[,samplevaraibles])

#1-3
norm.var <-!unlist(lapply(Liver_Norm213, is.factor))
name.var<-names(Liver_Norm213[,!unlist(lapply(Liver_Norm213, is.factor))])
str(Liver_Norm213)
liver.norm.df213 = cbind(Liver_Norm213[,norm.var],
                         class.ind(Liver_Norm213$REGION),
                         class.ind(Liver_Norm213$DIAB),
                         class.ind(Liver_Norm213$DIABETES_DON),
                         class.ind(Liver_Norm213$MED_COND_TRR),
                         class.ind(Liver_Norm213$gendermat),
                         class.ind(Liver_Norm213$PSTATUS))

names(liver.norm.df213) = c(name.var,
                            paste("REGION",c(1:11),sep = ""),
                            paste("DIAB", c(1:6),sep = ""),
                            paste("DIAB_DON", c(1:3),sep = ""),
                            paste("MED_COND_TRR",c(1:3),sep = ""),
                            paste("gendermat",c(0,1),sep = ""),
                            paste("PSTATUS",c(0,1),sep = ""))

training213 <- sample(1:nrow(liver.norm.df213), dim(liver.norm.df213)[1]*.60)
traindata213 <- liver.norm.df213[training213,]
Validdata213 <- liver.norm.df213[-training213,]
head(traindata213)
str(Liver.norm.df213)
attach(liver.norm.df213)
colnames(traindata213)
#1 hidden layer, 2 nodes
nn <- neuralnet::neuralnet(PSTATUS0+PSTATUS1~., data = traindata213,linear.output = F, hidden = 2, learningrate = 1) 
plot(nn)
train.pred = neuralnet::compute(nn, traindata213[,-c(31,32)])
train.class = apply(train.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(train.class), as.factor(traindata213$PSTATUS1)) #0.7214

valid.pred = neuralnet::compute(nn, Validdata213[,-c(31,32)])
valid.class = apply(valid.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(valid.class), as.factor(Validdata213$PSTATUS1)) #0.7140

#1 hidden layer, 3 nodes
nn <- neuralnet::neuralnet(PSTATUS0+PSTATUS1~., data = traindata213,linear.output = F, hidden = 3, learningrate = 1)
plot(nn)
train.pred = neuralnet::compute(nn, traindata213[,-c(31,32)])
train.class = apply(train.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(train.class), as.factor(traindata213$PSTATUS1)) #0.7201

valid.pred = neuralnet::compute(nn, Validdata213[,-c(31,32)])
valid.class = apply(valid.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(valid.class), as.factor(Validdata213$PSTATUS1)) #0.7111                  


###########################################################################################################################################################
#CLASSIFICATION MODELS FOR SURVIVAL AFTER ONE TO FIVE YEARS AFTER SURGERY 
set.seed(1)
training <- sample(1:nrow(Liver_Norm15), dim(Liver_Norm15)[1]*.6)
trainData <- Liver_Norm15[training,]
testData <- Liver_Norm15[-training,]

###############
#LOGISTIC MODEL 
log.mod = glm(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
summary(log.mod)
log.mod.pred = predict(log.mod,newdata = testData, type = "response")
log.mod.prob = ifelse(log.mod.pred >.5, "1", "0") 
table(log.mod.prob,testData$PSTATUS)
mean(na.omit(log.mod.prob == testData$PSTATUS)) #0.7153623  
mean(na.omit(log.mod.prob != testData$PSTATUS)) #0.2856377

#same model but with blood type
log.mod = glm(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
summary(log.mod)
log.mod.pred = predict(log.mod,newdata = testData, type = "response")
log.mod.prob = ifelse(log.mod.pred >.5, "1", "0") 
table(log.mod.prob,testData$PSTATUS)
mean(na.omit(log.mod.prob == testData$PSTATUS)) #0.7165766
mean(na.omit(log.mod.prob != testData$PSTATUS)) #0.2854234

###############
#LDA MODEL 
lda.mod1 = lda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
lda.pred1=predict(lda.mod1, testData)
predictedClass=lda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) # 0.7255828
mean(na.omit(predictedClass!=testData$PSTATUS)) # 0.2744172

#same model with blood type 
lda.mod1 = lda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
lda.pred1=predict(lda.mod1, testData)
predictedClass=lda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) # 0.7144695
mean(na.omit(predictedClass!=testData$PSTATUS)) # 0.2855305
###############
#QDA MODEL 
qda.mod1 = qda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC,data = trainData,family = binomial)
qda.pred1=predict(qda.mod1, testData)
predictedClass=qda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) #0.7045814
mean(na.omit(predictedClass!=testData$PSTATUS)) #0.2855305

qda.mod1 = qda(PSTATUS~gendermat+DIAB+DIABETES_DON+REGION+HGT_CM_TCR+
                 MED_COND_TRR+FINAL_MELD_SCORE+PTIME+LiverSize+LiverSizeDon+HGT_CM_DON_CALC+ABO_MAT,data = trainData,family = binomial)
qda.pred1=predict(qda.mod1, testData)
predictedClass=qda.pred1$class
table(predictedClass, testData$PSTATUS)
mean(na.omit(predictedClass==testData$PSTATUS)) #0.7046806
mean(na.omit(predictedClass!=testData$PSTATUS)) #0.2943194

######################
#NN FOR CLASSIFICATION  
samplevaraibles <- c("FINAL_MELD_SCORE","REGION","HGT_CM_CALC"
                     ,"HGT_CM_DON_CALC"
                     ,"LiverSize","LiverSizeDon","DIAB"
                     ,"DIABETES_DON","MED_COND_TRR"
                     ,"gendermat","PSTATUS")

Liver_Norm215 <- na.omit(Liver_Norm15[,samplevaraibles])

#1-5
norm.var <-!unlist(lapply(Liver_Norm215, is.factor))
name.var<-names(Liver_Norm215[,!unlist(lapply(Liver_Norm215, is.factor))])
str(Liver_Norm215)
liver.norm.df215 = cbind(Liver_Norm215[,norm.var],
                         class.ind(Liver_Norm215$REGION),
                         class.ind(Liver_Norm215$DIAB),
                         class.ind(Liver_Norm215$DIABETES_DON),
                         class.ind(Liver_Norm215$MED_COND_TRR),
                         class.ind(Liver_Norm215$gendermat),
                         class.ind(Liver_Norm215$PSTATUS))

names(liver.norm.df215) = c(name.var,
                            paste("REGION",c(1:11),sep = ""),
                            paste("DIAB", c(1:6),sep = ""),
                            paste("DIAB_DON", c(1:3),sep = ""),
                            paste("MED_COND_TRR",c(1:3),sep = ""),
                            paste("gendermat",c(0,1),sep = ""),
                            paste("PSTATUS",c(0,1),sep = ""))

training215 <- sample(1:nrow(liver.norm.df215), dim(liver.norm.df215)[1]*.60)
traindata215 <- liver.norm.df215[training215,]
Validdata215 <- liver.norm.df215[-training215,]
#1 hidden layer, 2 nodes
nn <- neuralnet::neuralnet(PSTATUS0+PSTATUS1~., data = traindata215,linear.output = F, hidden = 2, learningrate = 1) 
plot(nn)
train.pred = neuralnet::compute(nn, traindata215[,-c(31,32)])
train.class = apply(train.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(train.class), as.factor(traindata215$PSTATUS1)) #0.7201

valid.pred = neuralnet::compute(nn, Validdata215[,-c(31,32)])
valid.class = apply(valid.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(valid.class), as.factor(Validdata215$PSTATUS1)) #0.7111

#1 hidden layer, 3 nodes
nn <- neuralnet::neuralnet(PSTATUS0+PSTATUS1~., data = traindata215,linear.output = F, hidden = 3, learningrate = 1)
plot(nn)
train.pred = neuralnet::compute(nn, traindata215[,-c(31,32)])
train.class = apply(train.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(train.class), as.factor(traindata215$PSTATUS1)) #0.7253

valid.pred = neuralnet::compute(nn, Validdata213[,-c(31,32)])
valid.class = apply(valid.pred$net.result,1,which.max)-1
confusionMatrix(as.factor(valid.class), as.factor(Validdata213$PSTATUS1)) #0.7146

########################################################################################################################
# Naive Bayes Model
########################################################################################################################

############################# MODELS FOR SURVIVAL AFTER ONE TO THREE YEARS AFTER SURGERY ################################

liver.sample = liver.df2[,c("PSTATUS","gendermat","DIAB","DIABETES_DON","REGION","MED_COND_TRR"
                               ,"HGT_CM_DON_CALC","HGT_CM_CALC","LiverSize","LiverSizeDon","FINAL_MELD_SCORE")]

liver.sample=na.omit(liver.sample)

lq<-quantile(liver.sample$LiverSize)
liver.sample$LiverSize <- ifelse(liver.sample$LiverSize <lq[2],1,ifelse(liver.sample$LiverSize>lq[4],3,2))
liver.sample$LiverSize <- as.factor(liver.sample$LiverSize)

lq1<-quantile(liver.sample$LiverSizeDon)
liver.sample$LiverSizeDon <- ifelse(liver.sample$LiverSizeDon <lq1[2],1,ifelse(liver.sample$LiverSizeDon>lq1[4],3,2))
liver.sample$LiverSizeDon <- as.factor(liver.sample$LiverSizeDon)

hgt<-quantile(liver.sample$HGT_CM_CALC)
liver.sample$HGT_CM_CALC <- ifelse(liver.sample$HGT_CM_CALC <hgt[2],1,ifelse(liver.sample$HGT_CM_CALC>hgt[4],3,2))
liver.sample$HGT_CM_CALC <- as.factor(liver.sample$HGT_CM_CALC)

hgt1<-quantile(liver.sample$HGT_CM_DON_CALC)
liver.sample$HGT_CM_DON_CALC <- ifelse(liver.sample$HGT_CM_DON_CALC <hgt1[2],1,ifelse(liver.sample$HGT_CM_DON_CALC>hgt1[4],3,2))
liver.sample$HGT_CM_DON_CALC <- as.factor(liver.sample$HGT_CM_DON_CALC)

meld<-quantile(liver.sample$FINAL_MELD_SCORE)
liver.sample$FINAL_MELD_SCORE <- ifelse(liver.sample$FINAL_MELD_SCORE <meld[2],1,ifelse(liver.sample$FINAL_MELD_SCORE>meld[4],3,2))
liver.sample$FINAL_MELD_SCORE <- as.factor(liver.sample$FINAL_MELD_SCORE)

str(liver.sample)

set.seed(1)
training <- sample(1:nrow(liver.sample), dim(liver.sample)[1]*.6)
trainData <- liver.sample[training,]
testData <- liver.sample[-training,]

model <- naive_bayes(PSTATUS ~ ., data = trainData, usekernel = T)
model

plot(model)

# Confusion Matrix - train data
p <- predict(model, trainData)
confusionMatrix(as.factor(p),trainData$PSTATUS) # Accuracy : 0.7122  

# Confusion Matrix - test data  
p2 <- predict(model, testData)
confusionMatrix(as.factor(p2),testData$PSTATUS) # Accuracy : 0.7177  


######################################## MODELS FOR SURVIVAL AFTER ONE TO FIVE YEARS AFTER SURGERY ################################


liver.sample2 = liver.df3[,c("PSTATUS","gendermat","DIAB","DIABETES_DON","REGION","MED_COND_TRR"
                             ,"HGT_CM_DON_CALC","HGT_CM_CALC","LiverSize","LiverSizeDon","FINAL_MELD_SCORE")]

liver.sample2=na.omit(liver.sample2)

lq<-quantile(liver.sample2$LiverSize)
liver.sample2$LiverSize <- ifelse(liver.sample2$LiverSize <lq[2],1,ifelse(liver.sample2$LiverSize>lq[4],3,2))
liver.sample2$LiverSize <- as.factor(liver.sample2$LiverSize)

lq1<-quantile(liver.sample2$LiverSizeDon)
liver.sample2$LiverSizeDon <- ifelse(liver.sample2$LiverSizeDon <lq1[2],1,ifelse(liver.sample2$LiverSizeDon>lq1[4],3,2))
liver.sample2$LiverSizeDon <- as.factor(liver.sample2$LiverSizeDon)

hgt<-quantile(liver.sample2$HGT_CM_CALC)
liver.sample2$HGT_CM_CALC <- ifelse(liver.sample2$HGT_CM_CALC <hgt[2],1,ifelse(liver.sample2$HGT_CM_CALC>hgt[4],3,2))
liver.sample2$HGT_CM_CALC <- as.factor(liver.sample2$HGT_CM_CALC)

hgt1<-quantile(liver.sample2$HGT_CM_DON_CALC)
liver.sample2$HGT_CM_DON_CALC <- ifelse(liver.sample2$HGT_CM_DON_CALC <hgt1[2],1,ifelse(liver.sample2$HGT_CM_DON_CALC>hgt1[4],3,2))
liver.sample2$HGT_CM_DON_CALC <- as.factor(liver.sample2$HGT_CM_DON_CALC)

meld<-quantile(liver.sample2$FINAL_MELD_SCORE)
liver.sample2$FINAL_MELD_SCORE <- ifelse(liver.sample2$FINAL_MELD_SCORE <meld[2],1,ifelse(liver.sample2$FINAL_MELD_SCORE>meld[4],3,2))
liver.sample2$FINAL_MELD_SCORE <- as.factor(liver.sample2$FINAL_MELD_SCORE)

str(liver.sample2)

set.seed(1)
training2 <- sample(1:nrow(liver.sample2), dim(liver.sample2)[1]*.6)
trainData2 <- liver.sample2[training2,]
testData2 <- liver.sample2[-training2,]

model <- naive_bayes(PSTATUS ~ ., data = trainData2, usekernel = T)
model

plot(model)

# Confusion Matrix - train data
p3 <- predict(model, trainData2)
confusionMatrix(as.factor(p3),trainData2$PSTATUS) #  Accuracy : 0.6268

# Confusion Matrix - test data  
p4 <- predict(model, testData2)
confusionMatrix(as.factor(p4),testData2$PSTATUS) #  Accuracy : 0.7735 

####################################################################################################################
# Random Forest Model
####################################################################################################################

###################################### MODELS FOR SURVIVAL AFTER ONE TO THREE YEARS AFTER SURGERY ##################


liver.sample3 = Liver_Norm13[,c("PSTATUS","gendermat","DIAB","DIABETES_DON","REGION","MED_COND_TRR"
                               ,"HGT_CM_DON_CALC","HGT_CM_CALC","LiverSize","LiverSizeDon","FINAL_MELD_SCORE")]

liver.sample3=na.omit(liver.sample3)
sum(is.na(liver.sample3))

set.seed(1)
training3 <- sample(1:nrow(liver.sample3), dim(liver.sample3)[1]*.6)
trainData3 <- liver.sample3[training3,]
testData3 <- liver.sample3[-training3,]

# Creating a Random Forest Model
model2 <- randomForest(PSTATUS ~ ., data = trainData3, ntree = 500, mtry = 3, importance = TRUE)
model2

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model2$err.rate), times=3),
  Type=rep(c("OOB", "0", "1"), each=nrow(model2$err.rate)),
  Error=c(model2$err.rate[,"OOB"], 
          model2$err.rate[,"0"], 
          model2$err.rate[,"1"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

# Predicting on training set
predTrain3 <- predict(model2, trainData3, type = "class")

# Checking classification accuracy
table(predTrain3, trainData3$PSTATUS)  

# Predicting on Validation set
predValid3 <- predict(model2, testData3, type = "class")

# Checking classification accuracy
mean(predValid3 == testData3$PSTATUS)                    
table(predValid3,testData3$PSTATUS)

# checking important variables
importance(model2)        
varImpPlot(model2)

# Using For loop to identify the right mtry for model
a=c()
i=5
for (i in 3:8) {
  model3 <- randomForest(PSTATUS ~ ., data = trainData3, ntree = 500, mtry = i, importance = TRUE)
  predValid <- predict(model2, testData3, type = "class")
  a[i-2] = mean(predValid == testData3$PSTATUS)
}

a

plot(3:8,a)

# Comparing Model of Random Forest with Decision Tree model 

model_dt = train(PSTATUS ~ ., data = trainData3, method = "rpart")
model_dt_1 = predict(model_dt, data = trainData3)
table(model_dt_1, trainData3$PSTATUS)
mean(model_dt_1 == trainData3$PSTATUS)

# Running on Validation Set
model_dt_vs = predict(model_dt, newdata = testData3)
table(model_dt_vs, testData3$PSTATUS)
mean(model_dt_vs == testData3$PSTATUS)

####################################### MODELS FOR SURVIVAL AFTER ONE TO FIVE YEARS AFTER SURGERY #############################################

liver.sample4 = Liver_Norm15[,c("PSTATUS","gendermat","DIAB","DIABETES_DON","REGION","MED_COND_TRR"
                                ,"HGT_CM_DON_CALC","HGT_CM_CALC","LiverSize","LiverSizeDon","FINAL_MELD_SCORE")]

liver.sample4=na.omit(liver.sample4)
sum(is.na(liver.sample4))

set.seed(1)
training4 <- sample(1:nrow(liver.sample4), dim(liver.sample4)[1]*.6)
trainData4 <- liver.sample4[training4,]
testData4 <- liver.sample4[-training4,]

# Creating a Random Forest Model
model4 <- randomForest(PSTATUS ~ ., data = trainData4, ntree = 100, mtry = 7, importance = TRUE)
model4

oob.error.data4 <- data.frame(
  Trees=rep(1:nrow(model4$err.rate), times=3),
  Type=rep(c("OOB", "0", "1"), each=nrow(model4$err.rate)),
  Error=c(model4$err.rate[,"OOB"], 
          model4$err.rate[,"0"], 
          model4$err.rate[,"1"]))

ggplot(data=oob.error.data4, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))


# Predicting on training set
predTrain4 <- predict(model4, trainData4, type = "class")

# Checking classification accuracy
table(predTrain4, trainData4$PSTATUS)  

# Predicting on Validation set
predValid4 <- predict(model4, testData4, type = "class")

# Checking classification accuracy
mean(predValid4 == testData4$PSTATUS)                    
table(predValid4,testData4$PSTATUS)

# checking important variables
importance(model4)        
varImpPlot(model4)

# Using For loop to identify the right mtry for model
a=c()
i=5
for (i in 3:8) {
  model5 <- randomForest(PSTATUS ~ ., data = trainData4, ntree = 500, mtry = i, importance = TRUE)
  predValid4 <- predict(model4, testData4, type = "class")
  a[i-2] = mean(predValid4 == testData4$PSTATUS)
}

a

plot(3:8,a)

# Comparing Model of Random Forest with Decision Tree model 

model_dt4 = train(PSTATUS ~ ., data = trainData4, method = "rpart")
model_dt_14 = predict(model_dt4, data = trainData4)
table(model_dt_14, trainData4$PSTATUS)
mean(model_dt_14 == trainData4$PSTATUS)

# Running on Validation Set
model_dt_vs4 = predict(model_dt4, newdata = testData4)
table(model_dt_vs4, testData4$PSTATUS)
mean(model_dt_vs4 == testData4$PSTATUS)

