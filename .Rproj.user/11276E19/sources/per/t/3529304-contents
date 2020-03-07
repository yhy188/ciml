#' ---
#' title: "Turing Masterclass on Causal Machine learning: session 1"
#' author: "Karla Diazordaz"
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

data<-read.csv(paste("datasets/CausalML_Masterclass_n",n[1],"_",n_W[2],"W.csv",sep=""),row.names = 1)


#'take a look at the data
dim(data)
#'the covariates are named W (there are two different datasets with 25 Ws or 5 Ws)
#'the treatment is A
#' the (binary) outcome is Y
head(data)

#' we will consider all available covariates as potential controls for confounders
covariates<-subset(data, select=-Y)

#' Parametric g computation 
#' using a logistic regression with main terms only 
y_form  = paste("Y","~",paste(names(covariates),collapse=" + "), collapse = "") 
Y.mod <- glm(as.formula(y_form), data=data, family=binomial())
summary(Y.mod) # show results


#' Now, to predict potential outcomes
#' we fiirst copy the observed data
exp.data<- unexp.data <- data
#' and the set A=1 in the exp.data and A=0 in unexp.data
exp.data$A <-1
unexp.data$A <- 0

#'predict the POs as if exposed
predict.outcome.exp<- predict(Y.mod, newdata=exp.data, type='response')
head(predict.outcome.exp)

#'predict the POs as if unexposed
predict.outcome.unexp<- predict(Y.mod, newdata=unexp.data, type='response')
head(predict.outcome.unexp)
 
 # Point eestimate of E[E(Y|A=1,W)]- E[E(Y|A=0,W)]
param.gcomp <- mean(predict.outcome.exp - predict.outcome.unexp)
 
# bootstrap to get CIs, the following fucntion does this
ACE_gcomp_boot <- function(y_form, data, indices){
  data_o = data[indices,]
  fit.out.reg <- glm(as.formula(y_form), data=data_o, family=binomial())
  Y1 = predict(fit.out.reg, newdata = data.frame(A = 1, data_o), type="response")
  Y0 = predict(fit.out.reg, newdata = data.frame(A = 0, data_o), type="response")
  ACE=mean(Y1)-mean(Y0)
  return(ACE=ACE) 
}

#### only run if you  have the time 
ACE.boot <- boot(y_form=y_form, data=data, ACE_gcomp_boot, R = 1999)
boot.ci.ace<-boot.ci(ACE.boot, conf = 0.95, type = "perc")

 
 ###IPW 
#' IPW models A on W, so 
#' remove A from the set of covariates
regresors<- subset(covariates, select=-A)
#' lets do again main terms logistic regression
g_form  = paste("A", "~", paste(names(regresors),collapse=" + "), collapse = "") 
g.mod<- glm(as.formula(g_form), data=data, family=binomial())
ps <- predict(g.mod, type="response")
wt <- data$A/ps + (1-data$A)/(1-ps)

#' feel free to play with the parametric PS model and change the form or vars
summary(wt)
hist(ps)
#' estimate of E[E(Y|A=1,W)]
mean(as.numeric(data$A==1)*wt*data$Y)
#' estimate of E[E(Y|A=0,W)]
mean(as.numeric(data$A==0)*wt*data$Y)
#' our point estimate for Psi(P)
IPW<- mean(as.numeric(data$A==1)*wt*data$Y) -mean(as.numeric(data$A==0)*wt*data$Y)
IPW

#' you can bootstrap the above to get CIs (code not provided)

#' the IPW can also be writen as 
H.AW.param<- (as.numeric(data$A==1)/ps) - (as.numeric(data$A==0)/(1-ps))
IPW.param<-mean(H.AW.param*data$Y)


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
#' the IPW estimator can also be coded as 
H.AW<- as.numeric(data$A==1)/g1W - as.numeric(data$A==0)/g0W
summary(H.AW)

SL.IPW<-mean(H.AW*data$Y)
#' Warning: no way of doing valid inference


#' We now save the SL fits, because we're going to use them in the second computer lab

Q=as.data.frame(cbind(SL.predict.outcome.obs,SL.predict.outcome.unexp,SL.predict.outcome.exp))
names(Q)<-cbind("QAW","Q0W","Q1W")

save(Q, file = paste("Q0_",n[1],"_",n_W[2],".RData"))
save(g1W, file = paste("g1W_",n[1],"_",n_W[2],".RData"))
#' Make sure you change the values of n[j] and n_W[k] accordingly so you don't overwrite

results<-cbind(param.gcomp,IPW,SL.plugin.gcomp,SL.IPW)
results


#'Repeat the above steps,now for the other settings, i.e.
#'n=10,000 and 25 covariates and 
#'n=500 and 10,000 with only the 5 covariates that are also the confounders

#' Make sure you change the name of the result object for each setting, or make a note

#' Now we can look at all the results and compare them to the true
#' load the true ATE file 
trueATE<-read.csv("datasets/TrueATEs.csv",row.names = 1)
trueATE
#' you can now compare with the true ATE for the appropriate number of covariates


#' In general, the plug-in bias is larger in small sample
#' you can see that the larger n has less plug-in bias
#' including variables that are not confounders also introduces bias
#' the parametric g-comp is biased, because of mis-specification of the Y model
#' the IPW is less biased, because the g model is not too badly mis-specified in this case

#


