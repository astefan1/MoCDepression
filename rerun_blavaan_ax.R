# ==============================================================================
# BAYESIAN RANDOM INTERCEPT CROSS-LAGGED PANEL MODEL: ATTACHMENT ANXIETY
# ==============================================================================

rm(list=ls())

############################## LOAD PACKAGES ###################################

library(dplyr)
library(blavaan)
library(datawizard)
library(bayesplot)

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
data <- data %>% rescale(c(ds_1, ds_2, ds_3, ds_4), to = c(0, 1), range = c(0, 63))
data <- data %>% rescale(c(es_1, es_2, es_3, es_4, cr_1, cr_2, cr_3, cr_4), to = c(0, 1), range = c(1, 7))
data <- data %>% rescale(c(av_1, av_3, av_4, ax_1, ax_3, ax_4), to = c(0, 1), range = c(1, 5))

# check data
head(data)

################### MODEL 1: UNCONSTRAINED #####################################

RICLPM_ax_unconstrained <- '
  # Create between components (random intercepts)
  RI_ax =~ 1*ax_1 + 1*ax_3 + 1*ax_4  ## attachment anxiety
  RI_ds =~ 1*ds_1 + 1*ds_3 + 1*ds_4  ## depressive symptoms

  # Create within-person centered variables
  w_ax1 =~ 1*ax_1
  w_ax3 =~ 1*ax_3
  w_ax4 =~ 1*ax_4
  w_ds1 =~ 1*ds_1
  w_ds3 =~ 1*ds_3
  w_ds4 =~ 1*ds_4

  ### add autogregressive effects and cross-lagged effects between within centred variables
  w_ax3 ~ w_ax1 + w_ds1 
  w_ds3 ~ w_ax1 + w_ds1 

  w_ax4 ~ w_ax3 + w_ds3  
  w_ds4 ~ w_ax3 + w_ds3 

# Estimate covariance between within-person centered variables at first wave
  w_ax1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_ax3 ~~ w_ds3
  w_ax4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_ax ~~ RI_ax
  RI_ds ~~ RI_ds
  RI_ax ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_ax1 ~~ w_ax1 # Variances
  w_ds1 ~~ w_ds1 
  w_ax3 ~~ w_ax3 # Residual variances
  w_ds3 ~~ w_ds3
  w_ax4 ~~ w_ax4
  w_ds4 ~~ w_ds4
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 
RICLPM_ax_unconstrained_fit <- blavaan(RICLPM_ax_unconstrained, 
                                       data = data,
                                       target = "stan",
                                       bcontrol = list(cores = 4),
                                       n.chains = 4,
                                       sample = 10000,
                                       group = "treat",
                                       meanstructure = T, # includes means from observes vars into model
                                       int.ov.free = T,
                                       control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20, warmup = 6000),
                                       dp = mydp)

save(RICLPM_ax_unconstrained_fit, file="RICLPM_ax_unconstrained_fit.RData")

summary(RICLPM_ax_unconstrained_fit)                                       

# extract draws from posterior
draws_ax_unconstrained <- blavInspect(RICLPM_ax_unconstrained_fit, "mcmc")
draws_ax_unconstrained <- do.call("rbind", draws_ax_unconstrained)

# check traceplots 
mcmc_obj <- blavInspect(RICLPM_ax_unconstrained_fit, "mcobj")

pdf(file="ax_traceplots_unconstrained.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj, pars = colnames(draws_ax_unconstrained)[grep("bet_sign", colnames(draws_ax_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_ax_unconstrained)[grep("Psi_cov", colnames(draws_ax_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_ax_unconstrained)[grep("Psi_var", colnames(draws_ax_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_ax_unconstrained)[grep("Nu_free", colnames(draws_ax_unconstrained))])
dev.off()

# plot posterior distributions
par(mfrow=c(11,5), mar=c(2,0,2,0))
for(i in 1:52){
  plot(density(draws_ax_unconstrained[,i]), xlab="", ylab="", main=colnames(draws_ax_unconstrained)[i], yaxt="n")
}

##################### MODEL 2: CONSTRAINED #####################################

RICLPM_ax_constrained <- '
  # Create between components (random intercepts)
  RI_ax =~ 1*ax_1 + 1*ax_3 + 1*ax_4  ## attachment anxiety
  RI_ds =~ 1*ds_1 + 1*ds_3 + 1*ds_4  ## depressive symptoms

  # Create within-person centered variables
  w_ax1 =~ 1*ax_1
  w_ax3 =~ 1*ax_3
  w_ax4 =~ 1*ax_4
  w_ds1 =~ 1*ds_1
  w_ds3 =~ 1*ds_3
  w_ds4 =~ 1*ds_4

  ### add autogregressive effects and cross-lagged effects between within centred variables
  w_ax3 ~ c(a1, a1) * w_ax1 + c(a2, a2) * w_ds1 
  w_ds3 ~ c(a3, a3) * w_ax1 + c(a4, a4) * w_ds1 

  w_ax4 ~ c(b1, b1) * w_ax3 + c(b2, b2) * w_ds3  
  w_ds4 ~ c(b3, b3) * w_ax3 + c(b4, b4) * w_ds3 

# Estimate covariance between within-person centered variables at first wave
  w_ax1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_ax3 ~~ w_ds3
  w_ax4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_ax ~~ RI_ax
  RI_ds ~~ RI_ds
  RI_ax ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_ax1 ~~ w_ax1 # Variances
  w_ds1 ~~ w_ds1 
  w_ax3 ~~ w_ax3 # Residual variances
  w_ds3 ~~ w_ds3
  w_ax4 ~~ w_ax4
  w_ds4 ~~ w_ds4
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 

RICLPM_ax_constrained_fit <- blavaan(RICLPM_ax_constrained, 
                                     data = data,
                                     target = "stan",
                                     bcontrol = list(cores = 4),
                                     n.chains = 4,
                                     sample = 10000,
                                     group = "treat",
                                     meanstructure = T, # includes means from observes vars into model
                                     int.ov.free = T,
                                     control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20, warmup = 6000),
                                     dp = mydp)

save(RICLPM_ax_constrained_fit, file = "RICLPM_ax_constrained_fit.RData")

summary(RICLPM_ax_constrained_fit)

# extract draws from posterior
draws_ax_constrained <- blavInspect(RICLPM_ax_constrained_fit, "mcmc")
draws_ax_constrained <- do.call("rbind", draws_ax_constrained)

# check traceplots 
mcmc_obj2 <- blavInspect(RICLPM_ax_constrained_fit, "mcobj")

pdf(file="ax_traceplots_constrained.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj2, pars = colnames(draws_ax_constrained)[grep("bet_sign", colnames(draws_ax_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_ax_constrained)[grep("Psi_cov", colnames(draws_ax_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_ax_constrained)[grep("Psi_var", colnames(draws_ax_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_ax_constrained)[grep("Nu_free", colnames(draws_ax_constrained))])
dev.off()

par(mfrow=c(11,5), mar=c(2,0,2,0))
for(i in 1:52){
  plot(density(draws_ax_constrained[,i]), xlab="", ylab="", main=colnames(draws_ax_constrained)[i], yaxt="n")
}

################### MODEL 3: TREATMENT AS PREDICTOR ############################

RICLPM_ax_treat <- '
  # Create between components (random intercepts)
  RI_ax =~ 1*ax_1 + 1*ax_3 + 1*ax_4  ## attachment anxiety
  RI_ds =~ 1*ds_1 + 1*ds_3 + 1*ds_4  ## depressive symptoms

  # Create within-person centered variables
  w_ax1 =~ 1*ax_1
  w_ax3 =~ 1*ax_3
  w_ax4 =~ 1*ax_4
  w_ds1 =~ 1*ds_1
  w_ds3 =~ 1*ds_3
  w_ds4 =~ 1*ds_4

  ### add autogregressive effects and cross-lagged effects between within centred variables
  w_ax3 ~ w_ax1 + w_ds1 + treat
  w_ds3 ~ w_ax1 + w_ds1 + treat

  w_ax4 ~ w_ax3 + w_ds3 + treat
  w_ds4 ~ w_ax3 + w_ds3 + treat

# Estimate covariance between within-person centered variables at first wave
  w_ax1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_ax3 ~~ w_ds3
  w_ax4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_ax ~~ RI_ax
  RI_ds ~~ RI_ds
  RI_ax ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_ax1 ~~ w_ax1 # Variances
  w_ds1 ~~ w_ds1 
  w_ax3 ~~ w_ax3 # Residual variances
  w_ds3 ~~ w_ds3
  w_ax4 ~~ w_ax4
  w_ds4 ~~ w_ds4
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 
RICLPM_ax_treat_fit <- blavaan(RICLPM_ax_treat, 
                               data = data,
                               target = "stan",
                               bcontrol = list(cores = 4),
                               n.chains = 4,
                               sample = 10000,
                               meanstructure = T, 
                               int.ov.free = T,
                               control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20, warmup = 6000),
                               dp = mydp)

save(RICLPM_ax_treat_fit, file = "RICLPM_ax_treat_fit.RData")

summary(RICLPM_ax_treat_fit)

# extract draws from posterior
draws_ax_treat_fit <- blavInspect(RICLPM_ax_treat_fit, "mcmc")
draws_ax_treat_fit <- do.call("rbind", draws_ax_treat_fit)

# check traceplots 
mcmc_obj3 <- blavInspect(RICLPM_ax_treat_fit, "mcobj")

pdf(file="ax_traceplots_treat.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj3, pars = colnames(draws_ax_treat_fit)[grep("bet_sign", colnames(draws_ax_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_ax_treat_fit)[grep("Psi_cov", colnames(draws_ax_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_ax_treat_fit)[grep("Psi_var", colnames(draws_ax_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_ax_treat_fit)[grep("Nu_free", colnames(draws_ax_treat_fit))])
dev.off()

par(mfrow=c(5,6), mar=c(2,0,2,0))
for(i in 1:30){
  plot(density(draws_ax_treat_fit[,i]), xlab="", ylab="", main=colnames(draws_ax_treat_fit)[i], yaxt="n")
}


