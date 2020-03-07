## ----opts, eval = TRUE, echo = FALSE, message = FALSE------------------------------------
options(width = 60)
options(tinytex.verbose = TRUE)


## ----ltmle_data--------------------------------------------------------------------------
# set seed for reproducibility & set sample size of 500
set.seed(212); n <- 500
# baseline variables
L0 <- data.frame(L01 = rnorm(n), L02 = rbinom(n, 1, 0.5))
# first treatment
gA0 <- plogis(0.2 * L0$L01 - 0.2 * L0$L02)
A0 <- rbinom(n = n, size = 1, prob = gA0)
# intermediate variable at time 1
L1 <- rnorm(n = n, mean = -A0 + L0$L01 - L0$L02, sd = 1)
# second treatment decision
gA1 <- plogis(0.2 * A0 - L1 + L0$L01)
A1 <- rbinom(n = n, size = 1, prob = gA1)
# intermediate variable at time 2
L2 <- rnorm(n = n, mean = -A0*A1 + 2*A1 - L0$L01 + L1, sd = 2)
# third treatment decision
gA2 <- plogis(A0 - A1 + 2*A0*A1 - L0$L01 + 0.2 * L1*L0$L02)
A2 <- rbinom(n = n, size = 1, prob = gA2)
# outcome
Y <- rnorm(n = n, mean = L0$L01 * L0$L02 * L2 - A0 - A1 - A2*A0*L2, sd = 2)
# put into a data frame
full_data <- data.frame(L0, A0 = A0, L1 = L1, 
                        A1 = A1, L2 = L2, A2 = A2, Y = Y)


## ----head_data---------------------------------------------------------------------------
head(full_data)


## ----truth, echo = TRUE, eval = TRUE-----------------------------------------------------
compute_truth <- function(n = 1e6, a0 = 1, a1 = 1, a2 = 1){
	set.seed(212)
	L0 <- data.frame(L01 = rnorm(n), L02 = rbinom(n, 1, 0.5))
	A0 <- rep(a0, n)
	L1 <- rnorm(n = n, mean = -A0 + L0$L01 - L0$L02, sd = 1)
	A1 <- rep(a1, n)
	L2 <- rnorm(n = n, mean = -A0*A1 + 2*A1 - L0$L01 + L1, sd = 2)
	A2 <- rep(a2, n)
	# outcome
	Y <- rnorm(n = n, mean = L0$L01 * L0$L02 * L2 - A0 - A1 - A2*A0*L2, sd = 2)
	# put into a data frame
	return(mean(Y))
}

# E[Y(1,1,1)] = compute_truth()
# E[Y(0,0,0)] = compute_truth(a0 = 0, a1 = 0, a2 = 0)


## ----ltmle_cens_data---------------------------------------------------------------------
set.seed(12)
# censoring prior to time 1 (1 = censored)
gC1 <- plogis(-2 + 0.05 * L0$L01)
C1 <- rbinom(n = n, size = 1, prob = gC1)
# censoring prior to time 2 (1 = censored)
gC2 <- plogis(-3 + 0.05 * A0 + 0.025 * L1 - 0.025 * L0$L02)
C2 <- rbinom(n = n, size = 1, prob = gC2)
# censoring prior to time 3 (1 = censored)
gC3 <- plogis(-3.5 + 0.05*A0*A1 - 0.025*L2 + 0.025 * L1)
C3 <- rbinom(n = n, size = 1, prob = gC3)
# make a cumulative indicator of censoring
anyC1 <- C1 == 1; anyC2 <- C1 == 1 | C2 == 1 
anyC3 <- C1 == 1 | C2 == 1 | C3 == 1
# censored data set
cens_data <- data.frame(L0, A0 = A0, 
               C1 = ltmle::BinaryToCensoring(is.censored = C1),
               L1 = ifelse(anyC1, NA, L1), A1 = ifelse(anyC1, NA, A1), 
               C2 = ltmle::BinaryToCensoring(is.censored = ifelse(anyC1, NA, C2)),
               L2 = ifelse(anyC2, NA, L2), A2 = ifelse(anyC2, NA, A2), 
               C3 = ltmle::BinaryToCensoring(is.censored = ifelse(anyC2, NA, C3)),
               Y = ifelse(anyC3, NA, Y))


## ----look_ltmle_cens_data----------------------------------------------------------------
head(cens_data, 9)


## ----load_drtmle, eval = TRUE, echo = TRUE, message = FALSE------------------------------
library(ltmle)
library(SuperLearner)


## ----simple_call_to_ltmle, echo = TRUE, eval = TRUE, message = FALSE, warning = FALSE----
set.seed(123)
ltmle_fit1 <- ltmle(
    data = full_data, 
    Anodes = c("A0", "A1", "A2"),
    Lnodes = c("L01","L02","L1","L2"),
    Ynodes = "Y",
    SL.library = list(Q = c("SL.earth", "SL.glm", "SL.mean"),
                      g = c("SL.earth", "SL.glm", "SL.mean")),
    stratify = FALSE, abar = list(treatment = c(1,1,1),
                                  control = c(0,0,0))
    )


## ----ltmle_sum---------------------------------------------------------------------------
summary(ltmle_fit1)	


## ----look_at_sl_weights------------------------------------------------------------------
# weights for outcome regressions, because we set stratify = FALSE, the output in 
# ltmle_fit1$fit$Q[[1]] is the same as in ltmle_fit1$fit$Q[[2]]
ltmle_fit1$fit$Q[[1]]


## ----look_at_sl_weights2-----------------------------------------------------------------
# weights for propensity scores, because we set stratify = FALSE, the output in 
# ltmle_fit1$fit$g[[1]] is the same as in ltmle_fit1$fit$g[[2]]

ltmle_fit1$fit$g[[1]]


## ----ex-wrt1, echo = TRUE----------------------------------------------------------------
tmp <- summary(ltmle_fit1)
EY1 <- tmp$effect.measures$treatment$estimate
EY1_ci <- tmp$effect.measures$treatment$CI
EY0 <- tmp$effect.measures$control$estimate
EY0_ci <- tmp$effect.measures$control$CI


## ----ex-wrt2, echo = TRUE----------------------------------------------------------------
w1 <- formatC(ltmle_fit1$fit$Q[[1]][[1]][,2], digits = 2, format = "f")
w2 <- formatC(ltmle_fit1$fit$Q[[1]][[2]][,2], digits = 2, format = "f")
w3 <- formatC(ltmle_fit1$fit$Q[[1]][[3]][,2], digits = 2, format = "f")


## ----simple_call_to_ltmle2, echo=TRUE, eval=TRUE, results='hide', message=FALSE, warning=FALSE----
set.seed(123)
ltmle_fit2 <- ltmle(
    data = cens_data, 
    Anodes = c("A0", "A1", "A2"),
    Lnodes = c("L01","L02","L1","L2"),
    Cnodes = c("C1","C2","C3"),
    Ynodes = "Y",
    SL.library = list(Q = c("SL.earth", "SL.glm", "SL.mean"),
                      g = c("SL.earth", "SL.glm", "SL.mean")),
    stratify = FALSE, abar = list(treatment = c(1,1,1),
                                  control = c(0,0,0))
    )


## ----ltmle_sum2--------------------------------------------------------------------------
summary(ltmle_fit2)


## ----define_rule-------------------------------------------------------------------------
rule1 <- function(pt_data){
	# all patients start on control
	A0 <- 0
	# patients get treatment at time 1 if L1 > -1
	# set patients with missing L1 to NA
	if(!is.na(pt_data$L1)){
		A1 <- ifelse(pt_data$L1 > -1, 1, 0)
	}else{
		A1 <- NA
	}
	# patients get treatment at time 2 if L2 > -1
	# set patients with missing L2 to NA
	if(!is.na(pt_data$L1)){
		A2 <- ifelse(pt_data$L2 > -1, 1, 0)
	}else{
		A2 <- NA
	}
	return(c(A0,A1,A2))
}


## ----define_rule2------------------------------------------------------------------------
rule2 <- function(pt_data){
	# all patients start on control
	A0 <- 0
	# and stay on control unless censored
	A1 <- ifelse(is.na(pt_data$L1), NA, 0)
	A2 <- ifelse(is.na(pt_data$L2), NA, 0)
	return(c(A0,A1,A2))
}


## ----simple_call_to_ltmle3, echo=TRUE, eval=TRUE, results='hide', message=FALSE, warning=FALSE----
set.seed(123)
ltmle_fit3 <- ltmle(
    data = cens_data, 
    Anodes = c("A0", "A1", "A2"),
    Lnodes = c("L01","L02","L1","L2"),
    Cnodes = c("C1","C2","C3"),
    Ynodes = "Y", stratify = FALSE, 
    SL.library = list(Q = c("SL.earth", "SL.glm", "SL.mean"),
                      g = c("SL.earth", "SL.glm", "SL.mean")),
    rule = list(treatment = rule1, control = rule2)
    )


## ----summary_dr_ltmle--------------------------------------------------------------------
summary(ltmle_fit3)

