# ==============================================================================
# BAYES FACTORS FOR DEVELOPMENT OF MECHANISMS OF CHANGE
# ==============================================================================

rm(list=ls())

############################## LOAD PACKAGES ###################################

library(dplyr)
library(datawizard)
library(tidyr)
library(reshape2)
library(brms)
library(blavaan)
library(logspline)
library(nlme)

############################ LOAD MODEL FITS ###################################

load("../model_files/RICLPM_es_treat_fit.RData")
load("../model_files/RICLPM_cr_treat_fit.RData")
load("../model_files/RICLPM_av_treat_fit.RData")
load("../model_files/RICLPM_ax_treat_fit.RData")

############################ DATA WRANGLING ####################################

# load data
data <- read.csv("dat_TBB.csv", sep = ";")

# recode sex to binary variable
data$sex2[data$sex2 == "weiblich"] <- 0
data$sex2[data$sex2 == "m\xe4nnlich" ] <- 1

# remove IST condition
data <- data %>% subset(data$treatment != "IST")

# keep only variables needed for our analysis
data <- data %>% select(c("ID", "treatment","sex2", "age", "single_sum_total_n", "group_sum",
                          "t0_bdi", "t4_bdi", "t7_bdi", "t8_bdi",     # depressive symptoms
                          "t0_erq_re", "t4_erq_re", "t7_erq_re", "t8_erq_re",   # cognitive reappraisal
                          "t0_erq_su", "t4_erq_su", "t7_erq_su", "t8_erq_su",   # emotion suppression
                          "t0_rsqb_an", "t7_rsqb_an", "t8_rsqb_an",   # attachment avoidance
                          "t0_rsqb_at", "t7_rsqb_at", "t8_rsqb_at"))  # attachment anxiety

# rename variables
data <- data %>% rename(treat = treatment, sex = sex2, 
                        ds_1 = t0_bdi, ds_2 = t4_bdi, ds_3 = t7_bdi, ds_4 = t8_bdi,  # depressive symptoms
                        es_1 = t0_erq_su, es_2 = t4_erq_su, es_3 = t7_erq_su, es_4 = t8_erq_su,   # emotion suppression
                        cr_1 = t0_erq_re, cr_2 = t4_erq_re, cr_3 = t7_erq_re, cr_4 = t8_erq_re,   # cognitive reappraisal
                        av_1 = t0_rsqb_an, av_3 = t7_rsqb_an, av_4 = t8_rsqb_an,    # attachment avoidance
                        ax_1 = t0_rsqb_at, ax_3 = t7_rsqb_at, ax_4 = t8_rsqb_at) # attachment anxiety

## rescale all variables from 0 to 1
#data <- data %>% rescale(c(ds_1, ds_2, ds_3, ds_4), to = c(0, 1), range = c(0, 63))
#data <- data %>% rescale(c(es_1, es_2, es_3, es_4, cr_1, cr_2, cr_3, cr_4), to = c(0, 1), range = c(1, 7))
#data <- data %>% rescale(c(av_1, av_3, av_4, ax_1, ax_3, ax_4), to = c(0, 1), range = c(1, 5))

########################## Regressions #########################################

prior <- set_prior("normal(0, 2)", class = "b")

#### Emotion suppression ####

ESdata_long <- melt(data[, c("ID", "es_1", "es_2", "es_3", "es_4", "treat")], id.vars=c("ID", "treat"))
ESdata_long$variable <- as.factor(ESdata_long$variable)
ESdata_long$treat <- as.factor(ESdata_long$treat)
colnames(ESdata_long) <- c("ID", "treat", "time", "ESscore")

# Frequentist regression models
lm_es_noInt <- summary(lme(ESscore ~ treat + time, data = ESdata_long[!is.na(ESdata_long$ESscore),], random=~1|ID, method="REML"))
lm_es_withInt <- summary(lme(ESscore ~ treat * time, data = ESdata_long[!is.na(ESdata_long$ESscore),], random=~1|ID, method="REML"))

hyps_Interaction <- c("treatST = 0", "timees_2 = 0", "timees_3 = 0", "timees_4 = 0", "treatST:timees_2 = 0", "treatST:timees_3 = 0", "treatST:timees_4 = 0")
hyps_noInteraction <- c("treatST = 0", "timees_2 = 0", "timees_3 = 0", "timees_4 = 0")

# Bayesian regression models
bayesian_es_noInteraction <- summary(brm(ESscore ~ treat + time + (1|ID), data = na.omit(ESdata_long), iter=10000, sample_prior="yes", prior=prior))
bayesian_es_Interaction <- summary(brm(ESscore ~ treat * time + (1|ID), data = na.omit(ESdata_long), iter=10000, sample_prior="yes", prior=prior))

bf_es_noInt <- hypothesis(bayesian_es_noInteraction, hyps_noInteraction)
bf_es_withInt <- hypothesis(bayesian_es_Interaction, hyps_Interaction)

coef_tests_ES_noInt <- round(cbind(coef(lm_es_noInt), BF10 = c(NA, 1/bf_es_noInt$hypothesis$Evid.Ratio)), 3)
coef_tests_ES_Int <- round(cbind(coef(lm_es_withInt), BF10 = c(NA, 1/bf_es_withInt$hypothesis$Evid.Ratio)), 3)

save(coef_tests_ES_noInt, file="../exploratory/coef_tests_ES_noInt.RData")
save(coef_tests_ES_Int, file="../exploratory/coef_tests_ES_Int.RData")

# RI-CLPM Intercept comparison
draws_es <- as.matrix(blavInspect(RICLPM_es_treat_fit, "mcobj"))
T2minT1_es <- draws_es[,"Nu_free[2]"]-draws_es[,"Nu_free[1]"]
T3minT1_es <- draws_es[,"Nu_free[3]"]-draws_es[,"Nu_free[1]"]
T4minT1_es <- draws_es[,"Nu_free[4]"]-draws_es[,"Nu_free[1]"]
set.seed(12345)
interceptDifference_es_prior <- logspline(rnorm(10000, 0, 32)*rnorm(10000, 0, 32))
T2minT1_es_posterior <- logspline(T2minT1_es)
T3minT1_es_posterior <- logspline(T3minT1_es)
T4minT1_es_posterior <- logspline(T4minT1_es)
BF01_T2minT1_es <- dlogspline(0, T2minT1_es_posterior)/ dlogspline(0, interceptDifference_es_prior)
BF01_T3minT1_es <- dlogspline(0, T3minT1_es_posterior)/ dlogspline(0, interceptDifference_es_prior)
BF01_T4minT1_es <- dlogspline(0, T4minT1_es_posterior)/ dlogspline(0, interceptDifference_es_prior)

BF_RICLPM <- cbind(BF01_T2minT1_es, BF01_T3minT1_es, BF01_T4minT1_es)

#### Cognitive reappraisal ####

CRdata_long <- melt(data[, c("ID", "cr_1", "cr_2", "cr_3", "cr_4", "treat")], id.vars=c("ID", "treat"))
CRdata_long$variable <- as.factor(CRdata_long$variable)
CRdata_long$treat <- as.factor(CRdata_long$treat)
colnames(CRdata_long) <- c("ID", "treat", "time", "CRscore")

lm_cr_noInt <- summary(lme(CRscore ~ treat + time, data = CRdata_long[!is.na(CRdata_long$CRscore),], random=~1|ID, method="REML"))
lm_cr_withInt <- summary(lme(CRscore ~ treat * time, data = CRdata_long[!is.na(CRdata_long$CRscore),], random=~1|ID, method="REML"))

hyps_Interaction <- c("treatST = 0", "timecr_2 = 0", "timecr_3 = 0", "timecr_4 = 0", "treatST:timecr_2 = 0", "treatST:timecr_3 = 0", "treatST:timecr_4 = 0")
hyps_noInteraction <- c("treatST = 0", "timecr_2 = 0", "timecr_3 = 0", "timecr_4 = 0")

bayesian_cr_noInteraction <- brm(CRscore ~ treat + time + (1|ID), data = na.omit(CRdata_long), iter=10000, sample_prior="yes", prior=prior)
bayesian_cr_Interaction <- brm(CRscore ~ treat * time + (1|ID), data = na.omit(CRdata_long), iter=10000, sample_prior="yes", prior=prior)

bf_cr_noInt <- hypothesis(bayesian_cr_noInteraction, hyps_noInteraction)
bf_cr_withInt <- hypothesis(bayesian_cr_Interaction, hyps_Interaction)

coef_tests_CR_noInt <- round(cbind(coef(lm_cr_noInt), BF10 = c(NA, 1/bf_cr_noInt$hypothesis$Evid.Ratio)), 3)
coef_tests_CR_Int <- round(cbind(coefficients(lm_cr_withInt), BF10 = c(NA, 1/bf_cr_withInt$hypothesis$Evid.Ratio)), 3)

save(coef_tests_CR_noInt, file="../exploratory/coef_tests_CR_noInt.RData")
save(coef_tests_CR_Int, file="../exploratory/coef_tests_CR_Int.RData")

# RI-CLPM Intercept comparison
draws_cr <- as.matrix(blavInspect(RICLPM_cr_treat_fit, "mcobj"))
T2minT1_cr <- draws_cr[,"Nu_free[2]"]-draws_cr[,"Nu_free[1]"]
T3minT1_cr <- draws_cr[,"Nu_free[3]"]-draws_cr[,"Nu_free[1]"]
T4minT1_cr <- draws_cr[,"Nu_free[4]"]-draws_cr[,"Nu_free[1]"]
set.seed(12345)
interceptDifference_cr_prior <- logspline(rnorm(10000, 0, 32)*rnorm(10000, 0, 32))
T2minT1_cr_posterior <- logspline(T2minT1_cr)
T3minT1_cr_posterior <- logspline(T3minT1_cr)
T4minT1_cr_posterior <- logspline(T4minT1_cr)
BF01_T2minT1_cr <- dlogspline(0, T2minT1_cr_posterior)/ dlogspline(0, interceptDifference_cr_prior)
BF01_T3minT1_cr <- dlogspline(0, T3minT1_cr_posterior)/ dlogspline(0, interceptDifference_cr_prior)
BF01_T4minT1_cr <- dlogspline(0, T4minT1_cr_posterior)/ dlogspline(0, interceptDifference_cr_prior)

BF_RICLPM <- rbind(BF_RICLPM,cbind(BF01_T2minT1_cr, BF01_T3minT1_cr, BF01_T4minT1_cr))

# Anxious attachment

AXdata_long <- melt(data[, c("ID", "ax_1", "ax_3", "ax_4", "treat")], id.vars=c("ID", "treat"))
AXdata_long$variable <- as.factor(AXdata_long$variable)
AXdata_long$treat <- as.factor(AXdata_long$treat)
colnames(AXdata_long) <- c("ID", "treat", "time", "AXscore")

lm_ax_noInt <- summary(lme(AXscore ~ treat + time, data = AXdata_long[!is.na(AXdata_long$AXscore),], random=~1|ID, method="REML"))
lm_ax_withInt <- summary(lme(AXscore ~ treat * time, data = AXdata_long[!is.na(AXdata_long$AXscore),], random=~1|ID, method="REML"))

hyps_Interaction <- c("treatST = 0", "timeax_3 = 0", "timeax_4 = 0", "treatST:timeax_3 = 0", "treatST:timeax_4 = 0")
hyps_noInteraction <- c("treatST = 0", "timeax_3 = 0", "timeax_4 = 0")

bayesian_ax_noInteraction <- brm(AXscore ~ treat + time + (1|ID), data = na.omit(AXdata_long), iter=10000, sample_prior="yes", prior=prior)
bayesian_ax_Interaction <- brm(AXscore ~ treat * time + (1|ID), data = na.omit(AXdata_long), iter=10000, sample_prior="yes", prior=prior)

bf_ax_noInt <- hypothesis(bayesian_ax_noInteraction, hyps_noInteraction)
bf_ax_withInt <- hypothesis(bayesian_ax_Interaction, hyps_Interaction)

coef_tests_AX_noInt <- round(cbind(coef(lm_ax_noInt), BF10 = c(NA, 1/bf_ax_noInt$hypothesis$Evid.Ratio)), 3)
coef_tests_AX_Int <- round(cbind(coef(lm_ax_withInt), BF10 = c(NA, 1/bf_ax_withInt$hypothesis$Evid.Ratio)), 3)

save(coef_tests_AX_noInt, file="../exploratory/coef_tests_AX_noInt.RData")
save(coef_tests_AX_Int, file="../exploratory/coef_tests_AX_Int.RData")

# RI-CLPM Intercept comparison
draws_ax <- as.matrix(blavInspect(RICLPM_ax_treat_fit, "mcobj"))
T3minT1_ax <- draws_ax[,"Nu_free[2]"]-draws_ax[,"Nu_free[1]"]
T4minT1_ax <- draws_ax[,"Nu_free[3]"]-draws_ax[,"Nu_free[1]"]
set.seed(12345)
interceptDifference_ax_prior <- logspline(rnorm(10000, 0, 32)*rnorm(10000, 0, 32))
T3minT1_ax_posterior <- logspline(T3minT1_ax)
T4minT1_ax_posterior <- logspline(T4minT1_ax)
BF01_T3minT1_ax <- dlogspline(0, T3minT1_ax_posterior)/ dlogspline(0, interceptDifference_ax_prior)
BF01_T4minT1_ax <- dlogspline(0, T4minT1_ax_posterior)/ dlogspline(0, interceptDifference_ax_prior)

BF_RICLPM <- rbind(BF_RICLPM,cbind(NA, BF01_T3minT1_ax, BF01_T4minT1_ax))

# Avoidant attachment

AVdata_long <- melt(data[, c("ID", "av_1", "av_3", "av_4", "treat")], id.vars=c("ID", "treat"))
AVdata_long$variable <- as.factor(AVdata_long$variable)
AVdata_long$treat <- as.factor(AVdata_long$treat)
colnames(AVdata_long) <- c("ID", "treat", "time", "AVscore")

lm_av_noInt <- summary(lme(AVscore ~ treat + time, data = AVdata_long[!is.na(AVdata_long$AVscore),], random=~1|ID, method="REML"))
lm_av_withInt <- summary(lme(AVscore ~ treat * time, data = AVdata_long[!is.na(AVdata_long$AVscore),], random=~1|ID, method="REML"))

hyps_Interaction <- c("treatST = 0", "timeav_3 = 0", "timeav_4 = 0", "treatST:timeav_3 = 0", "treatST:timeav_4 = 0")
hyps_noInteraction <- c("treatST = 0", "timeav_3 = 0", "timeav_4 = 0")

bayesian_av_noInteraction <- brm(AVscore ~ treat + time + (1|ID), data = na.omit(AVdata_long), iter=10000, sample_prior="yes", prior=prior)
bayesian_av_Interaction <- brm(AVscore ~ treat * time + (1|ID), data = na.omit(AVdata_long), iter=10000, sample_prior="yes", prior=prior)

bf_av_noInt <- hypothesis(bayesian_av_noInteraction, hyps_noInteraction)
bf_av_withInt <- hypothesis(bayesian_av_Interaction, hyps_Interaction)

coef_tests_AV_noInt <- round(cbind(coef(lm_av_noInt), BF10 = c(NA, 1/bf_av_noInt$hypothesis$Evid.Ratio)), 3)
coef_tests_AV_Int <- round(cbind(coef(lm_av_withInt), BF10 = c(NA, 1/bf_av_withInt$hypothesis$Evid.Ratio)), 3)

save(coef_tests_AV_noInt, file="../exploratory/coef_tests_AV_noInt.RData")
save(coef_tests_AV_Int, file="../exploratory/coef_tests_AV_Int.RData")

# RI-CLPM Intercept comparison
draws_av <- as.matrix(blavInspect(RICLPM_av_treat_fit, "mcobj"))
T3minT1_av <- draws_av[,"Nu_free[2]"]-draws_av[,"Nu_free[1]"]
T4minT1_av <- draws_av[,"Nu_free[3]"]-draws_av[,"Nu_free[1]"]
set.seed(12345)
interceptDifference_av_prior <- logspline(rnorm(10000, 0, 32)*rnorm(10000, 0, 32))
T3minT1_av_posterior <- logspline(T3minT1_av)
T4minT1_av_posterior <- logspline(T4minT1_av)
BF01_T3minT1_av <- dlogspline(0, T3minT1_av_posterior)/ dlogspline(0, interceptDifference_av_prior)
BF01_T4minT1_av <- dlogspline(0, T4minT1_av_posterior)/ dlogspline(0, interceptDifference_av_prior)

BF_RICLPM <- rbind(BF_RICLPM,cbind(NA, BF01_T3minT1_av, BF01_T4minT1_av))

colnames(BF_RICLPM) <- c("BF01 T2-T1", "BF01 T3-T1", "BF01 T4-T1")
rownames(BF_RICLPM) <- c("Emotion Suppression", "Cognitive Reappraisal", "Anxious Attachment", "Avoidant Attachment")
BF_RICLPM <- round(BF_RICLPM, 6)

save(BF_RICLPM, file = "../exploratory/BF_RICLPM_Intercepts.RData")
