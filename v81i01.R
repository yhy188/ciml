## ----loadPkgs, echo=FALSE, message=FALSE, results='hide'-----------------
library("methods")
library("ltmle")
library("gam")
library("nnet")
library("glmnet")
library("SuperLearner")
library("ggplot2")

## ----exData1, echo=FALSE-------------------------------------------------
set.seed(2345)
n.Ex1 <- 500
L0.a <- runif(n.Ex1)
L0.b <- runif(n.Ex1)
L0.c <- runif(n.Ex1)
A0 <- rbinom(n.Ex1, 1, plogis(-1 + L0.a + L0.b))
C0 <- factor(c("uncensored", "censored")[((runif(n.Ex1) + 0.2 * A0) > 0.9) + 1])
L1.a <- L0.a + (1 - L0.b) * runif(n.Ex1) + 0.3 * A0  #from 0 to 2.3
L1.b <- -1 + L0.b + (0.5 + A0) * runif(n.Ex1)  #-1 to 1.5
Y1 <- (L0.c + (L0.a * A0 + runif(n.Ex1))/2)/2
A1 <- rbinom(n.Ex1, 1, plogis(-1.5 + L0.a + L0.b + A0))
C1 <- factor(c("uncensored", "censored")[((runif(n.Ex1) + 0.2 * A1) > 0.9) + 1])
Y2 <- (Y1 + (L1.a/2.3 * A1 + runif(n.Ex1))/2)/2

exData1 <- data.frame(L0.a, L0.b, L0.c, A0, C0, L1.a, L1.b, Y1, A1, C1, Y2)
exData1[exData1$C1 == "censored", 11] <- NA
exData1[exData1$C0 == "censored", 6:11] <- NA

## ----Ex1_EY1, dependson=c('SetNodes', 'exData1')-------------------------
Lnodes <- c("L1.a", "L1.b")
Anodes <- c("A0", "A1")
Cnodes <- c("C0", "C1")
Ynodes <- c("Y1", "Y2")
EY.11 <- ltmle(exData1, Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes, Ynodes = Ynodes, 
  abar = c(1, 1), estimate.time = FALSE)
print(summary(EY.11))

## ----Ex1_EY0, dependson=c('SetNodes', 'exData1', 'Ex1_EY1')--------------
ATE <- ltmle(exData1, Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes, Ynodes = Ynodes, 
  abar = list(treament = c(1, 1), control = c(0, 0)), estimate.time = FALSE)
print(summary(ATE))

## ----Ex1_EYd, dependsOn=c('SetNodes', 'exData1')-------------------------
d <- function(row) c(1, ifelse(row["L1.b"] > 0, 1, 0))

SL.lib <- c("SL.glm", "SL.stepAIC", "SL.nnet", "SL.gam", "SL.glmnet")
EY.d <- ltmle(exData1, Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes, Ynodes = Ynodes, 
  rule = d, SL.library = list(Q = NULL, g = SL.lib), estimate.time = FALSE)
print(summary(EY.d))

## ----exData2, echo=FALSE, dependson='exData1'----------------------------
exData2 <- transform(exData1, Y1 = ifelse(Y1 > 0.7, 1, 0), Y2 = ifelse((Y1 > 0.7) | 
  (Y2 > 0.55), 1, 0))

## ----Ex2_EYd, dependsOn=c('SetNodes', 'exData2')-------------------------
ATE.survival <- ltmle(exData2, Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes, 
  Ynodes = Ynodes, survivalOutcome = TRUE, abar = list(treatment = c(1, 1), control = c(0, 
    0)), estimate.time = FALSE)
print(summary(ATE.survival))

## ----MSMData, echo=FALSE-------------------------------------------------
data("sampleDataForLtmleMSM")

## ----Ex3_MSM-------------------------------------------------------------
Lnodes <- c("CD4_1", "CD4_2")
Anodes <- c("A0", "A1", "A2")
Ynodes <- c("Y1", "Y2", "Y3")

D <- list(function(row) c(1, 1, 1), function(row) c(0, 1, 1), function(row) c(0, 
  0, 1), function(row) c(0, 0, 0))

summary.measures <- array(dim = c(4, 2, 3))
dimnames(summary.measures)[[2]] <- c("switch.time", "time")
summary.measures[, , 1] <- cbind(0:3, rep(1, 4))
summary.measures[, , 2] <- cbind(0:3, rep(2, 4))
summary.measures[, , 3] <- cbind(0:3, rep(3, 4))

MSM.estimates <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes = Anodes, Lnodes = Lnodes, 
  Ynodes = Ynodes, survivalOutcome = TRUE, regimes = D, summary.measures = summary.measures, 
  final.Ynodes = Ynodes, working.msm = "Y ~ time + pmax(time - switch.time, 0)", 
  estimate.time = FALSE)
print(summary(MSM.estimates))

## ----MSMplot, echo=FALSE, fig.width=6, fig.height=4----------------------
flat.sum.meas <- expand.grid(time = seq(1, 3, 1), switch.time = 0:3)
mm <- model.matrix(~time + pmax(time - switch.time, 0), data.frame(flat.sum.meas))
surv <- 1 - plogis(mm %*% summary(MSM.estimates)$cmat[, 1, drop = FALSE])

# offset survival so curves don't exactly overlap
df <- data.frame(surv = surv[, 1], surv.off = surv[, 1] - 0.001 * flat.sum.meas$switch.time, 
  flat.sum.meas)

df <- transform(df, plot.time = time + (switch.time - 1.5)/50)

p <- qplot(x = plot.time, y = surv.off, data = df, color = factor(switch.time), geom = "point", 
  xlab = "Time", ylab = "Estimated counterfactual survival curves\nunder different switch times") + 
  scale_color_discrete(name = "Switch time") + theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0.6, 1), breaks = c(0.6, 0.7, 0.8, 0.9, 1))

p

