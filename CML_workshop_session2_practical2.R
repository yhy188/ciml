#' ---
#' title: "Turing Masterclass on Causal Machine learning: session 2"
#' author: "Stijn Vansteelandt"
#' date: "March 3rd, 2020"
#' ---

#' Install packages and load them  ##########
#install.packages("SuperLearner")
#install.packages("xgboost")
#install.packages("tmle")

library(boot)
library(MASS) 
library(SuperLearner)
library(survey)

#' The data: there are two sizes n=500 and n=10,000 and 
#' two covariate settings with dim W either 5 or 25
#' the data generating mechanism is so that only the first 5 covariates are confounders
#' In real life, we wouldn't know this, so we will always adjust for all available variables
#' as potential confounders
n=c(500,10000)
n_W=c(5,25)


#' let's start with small sample and n=500 and W=25
#' load the data

#setwd("~/Dropbox/CML_Masterclass_2020/R scripts/Practical 1")
data<-read.csv(paste("datasets/CausalML_Masterclass_n",n[1],"_",n_W[1],"W.csv",sep=""),row.names = 1)


#'take a look at the data
dim(data)
#'the covariates are named W (there are two different datasets with 25 Ws or 5 Ws)
#'the treatment is A
#' the (binary) outcome is Y
head(data)

#' we will consider all available covariates as potential controls for confounders
covariates<-subset(data, select=-Y)

#' Now, to predict potential outcomes
#' we fiirst copy the observed data
exp.data<- unexp.data <- data
#' and the set A=1 in the exp.data and A=0 in unexp.data
exp.data$A <-1
unexp.data$A <- 0

##### Super Learner 
#' load the Super learner package
library("SuperLearner")
#' Super Learner comes with a quite a few pre-specified learners
#' you can review the available leaners by 

listWrappers(what = "SL")

#' Specify the Super Learner library with the following algorithms:
#'  main terms logistic regression ("SL.glm"), 
#'  main terms logistic regression with all possible pairwise interactions ("SL.glm.interaction"),
#'  generalized additive models ("SL.gam") 
#'  lasso ("SL.glmnet")
#'  gradiant boosting ("SL.glmnet")
#'  
#'  In practice, we would want to use a larger library with a mixture of simple (e.g. parametric) and more aggressive libraries.
 #' note the learner names have to be in quotations

SL.library<- c("SL.glm", "SL.glm.interaction", "SL.xgboost", "SL.glmnet")
#SL.library<- c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.xgboost")

#' These should ideally be tested with multiple hyperparameter settings for each algorithm
#' which can be tuned using CV
#' now we only use the defaults. Check later you know which ones these are
#' you can do this by typing for example
# ?SL.xgboost

#### SL for the outcome regression ####
SL.outcome.regression<- SuperLearner(Y=data$Y, X=subset(data, select=-Y),
                                     SL.library=SL.library, family="binomial")
#' You can look at the Super learner object, to see how the alogorithms are weighted 
SL.outcome.regression

#' prediction for the actual exposure level received
SL.predict.outcome.obs<- predict(SL.outcome.regression, newdata=subset(data, select=-Y))$pred

#' predict the PO Y^1
SL.predict.outcome.exp<- predict(SL.outcome.regression, newdata=subset(exp.data, select=-Y))$pred
#' predict the PO Y^0
SL.predict.outcome.unexp<- predict(SL.outcome.regression, newdata=subset(unexp.data, select=-Y))$pred

##' SL g-computation
SL.plugin.gcomp<-mean(SL.predict.outcome.exp-SL.predict.outcome.unexp)

#' Warning:  no way of doing inference, bootstrap not valid when using ML

#' Note: we can also train the Super Learner separately in the exposed and unexposed 

#### SL for the prop score #####
SL.g<- SuperLearner(Y=data$A, X=subset(data, select=-c(A,Y)),
                         SL.library=SL.library, family="binomial")
#' You can look at the Super learner object, to see how the alogorithms are weighted 
SL.g


#' get the probability of getting the exposure
g1W <- SL.g$SL.predict
summary(g1W)
#' Look at the histogram of the weights
hist(1/SL.g$SL.predict)

#' noe the probability of being unexposed
g0W<- 1- g1W

#' We now save the SL fits, because we're going to use them in the second computer lab

Q=cbind(SL.predict.outcome.obs,SL.predict.outcome.unexp,SL.predict.outcome.exp)
colnames(Q)<-c("QAW","Q0W","Q1W")

save(Q, file = paste("Q0_",n[1],"_",n_W[2],".RData"))
save(g1W, file = paste("g1W_",n[1],"_",n_W[2],".RData"))

### TMLE

#' Plug-in AIPW
#' E(Y1)
mean((data$A/g1W)*(data$Y-Q[,"Q1W"])+Q[,"Q1W"])
#' E(Y0)
mean(((1-data$A)/g0W)*(data$Y-Q[,"Q0W"])+Q[,"Q0W"])
#' ATE = E(Y1)-E(Y0)
mean((data$A/g1W)*(data$Y-Q[,"Q1W"])+Q[,"Q1W"])-mean(((1-data$A)/g0W)*(data$Y-Q[,"Q0W"])+Q[,"Q0W"])

#' TMLE by hand 
#' E(Y1)
#' Constructing the clever covariate
H<-as.numeric(data$A/g1W)
#' Fitting a parametric extension model
model<-glm(data$Y~-1+H+offset(logit(Q[,"QAW"])),family=binomial)
summary(model)
#' Updating the predictions
Q1W.1<-plogis(logit(Q[,"Q1W"])+coef(model)[1]/g1W)
#' Estimating E(Y1)
mean(Q1W.1)

#' E(Y1)
#' Constructing the clever covariate
H<-as.numeric((1-data$A)/g0W)
#' Fitting a parametric extension model
model<-glm(data$Y~-1+H+offset(logit(Q[,"QAW"])),family=binomial)
summary(model)
#' Updating the predictions
Q0W.1<-plogis(logit(Q[,"Q0W"])+coef(model)[1]/g0W)
#' Estimating E(Y0)
mean(Q0W.1)

#' ATE = E(Y1)-E(Y0)
mean(Q1W.1)-mean(Q0W.1)

#' Alternative TMLE by hand 
#' E(Y1)
H<-as.numeric(data$A/g1W)
model<-glm(data$Y~offset(logit(Q[,"QAW"])),family=binomial,weight=H)
Q1W.1<-plogis(logit(Q[,"Q1W"])+coef(model)[1])
mean(Q1W.1)
#' E(Y0)
H<-as.numeric((1-data$A)/g0W)
model<-glm(data$Y~offset(logit(Q[,"QAW"])),family=binomial,weight=H)
Q0W.1<-plogis(logit(Q[,"Q0W"])+coef(model)[1])
mean(Q0W.1)
#' ATE = E(Y1)-E(Y0)
mean(Q1W.1)-mean(Q0W.1)

#' TMLE using the R package
result <- tmle(Y=data$Y,A=data$A,W=subset(data, select=-c(A,Y)), family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library)
summary(result)

#' Confidence interval for the plug-in AIPW estimator
IF<-((data$A/g1W)*(data$Y-Q[,"Q1W"])+Q[,"Q1W"])-(((1-data$A)/g0W)*(data$Y-Q[,"Q0W"])+Q[,"Q0W"])
mean(IF)
mean(IF)-qnorm(.975)*sd(IF)/sqrt(length(IF))
mean(IF)+qnorm(.975)*sd(IF)/sqrt(length(IF))

#'Repeat the above steps,now for the other settings, i.e.
#'n=10,000 and 25 covariates and 
#'n=500 and 10,000 with only the 5 covariates that are also the confounders

#' Make sure you change the name of the result object for each setting, or make a note
#' Make sure you change the values of n[j] and n_W[k] accordingly so you don't overwrite
#' Now we can look at all the results and compare them to the true
#' load the true ATE file 
trueATE<-read.csv("datasets/TrueATEs.csv",row.names = 1)
trueATE
#' you can now compare with the true ATE for the appropriate number of covariates
