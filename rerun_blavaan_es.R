# ==============================================================================
# BAYESIAN RANDOM INTERCEPT CROSS-LAGGED PANEL MODEL: EMOTION SUPPRESSION
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

RICLPM_es_unconstrained <- '
  # Create between components (random intercepts)
  RI_es =~ 1*es_1 + 1*es_2 + 1*es_3 + 1*es_4 ## emotion suppression
  RI_ds =~ 1*ds_1 + 1*ds_2 + 1*ds_3 + 1*ds_4 ## depressive symptoms

  # Create within-person centered variables
  w_es1 =~ 1*es_1
  w_es2 =~ 1*es_2
  w_es3 =~ 1*es_3
  w_es4 =~ 1*es_4
  w_ds1 =~ 1*ds_1
  w_ds2 =~ 1*ds_2
  w_ds3 =~ 1*ds_3
  w_ds4 =~ 1*ds_4

  ### add autogregressive effects and cross-lagged effects between within centred variables
  w_es2 ~ w_es1 + w_ds1
  w_ds2 ~ w_es1 + w_ds1

  w_es3 ~ w_es2 + w_ds2
  w_ds3 ~ w_es2 + w_ds2

  w_es4 ~ w_es3 + w_ds3
  w_ds4 ~ w_es3 + w_ds3

# Estimate covariance between within-person centered variables at first wave
  w_es1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_es2 ~~ w_ds2
  w_es3 ~~ w_ds3
  w_es4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_es ~~ RI_es
  RI_ds ~~ RI_ds
  RI_es ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_es1 ~~ w_es1 # Variances
  w_ds1 ~~ w_ds1 
  w_es2 ~~ w_es2 # Residual variances
  w_ds2 ~~ w_ds2 
  w_es3 ~~ w_es3 
  w_ds3 ~~ w_ds3 
  w_es4 ~~ w_es4 
  w_ds4 ~~ w_ds4 
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 
RICLPM_es_unconstrained_fit <- blavaan(RICLPM_es_unconstrained, 
                                       data = data,
                                       target = "stan",
                                       bcontrol = list(cores = 4),
                                       n.chains = 4,
                                       sample = 10000,
                                       group = "treat",
                                       meanstructure = T, # includes means from observes vars into model
                                       int.ov.free = T,
                                       control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20),
                                       dp = mydp)

save(RICLPM_es_unconstrained_fit, file = "RICLPM_es_unconstrained_fit.RData")

# extract draws from posterior
draws_es_unconstrained <- blavInspect(RICLPM_es_unconstrained_fit, "mcmc")
draws_es_unconstrained <- do.call("rbind", draws_es_unconstrained)

# check traceplots 
mcmc_obj <- blavInspect(RICLPM_es_unconstrained_fit, "mcobj")

pdf(file="es_traceplots_unconstrained.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj, pars = colnames(draws_es_unconstrained)[grep("bet_sign", colnames(draws_es_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_es_unconstrained)[grep("Psi_cov", colnames(draws_es_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_es_unconstrained)[grep("Psi_var", colnames(draws_es_unconstrained))])
mcmc_trace(mcmc_obj, pars = colnames(draws_es_unconstrained)[grep("Nu_free", colnames(draws_es_unconstrained))])
dev.off()

summary(RICLPM_es_unconstrained_fit)

# plot posterior distributions
par(mfrow=c(10,7), mar=c(2,0,2,0))
for(i in 1:70){
  plot(density(draws_es_unconstrained[,i]), xlab="", ylab="", main=colnames(draws_es_unconstrained)[i], yaxt="n")
}

########################### MODEL 2: CONSTRAINED ###############################

# Create constrained model
RICLPM_es_constrained <- '
# Create between components (random intercepts)
RI_es =~ 1*es_1 + 1*es_2 + 1*es_3 + 1*es_4 ## emotion suppression
RI_ds =~ 1*ds_1 + 1*ds_2 + 1*ds_3 + 1*ds_4 ## depressive symptoms

# Create within-person centered variables
w_es1 =~ 1*es_1
w_es2 =~ 1*es_2
w_es3 =~ 1*es_3
w_es4 =~ 1*es_4
w_ds1 =~ 1*ds_1
w_ds2 =~ 1*ds_2
w_ds3 =~ 1*ds_3
w_ds4 =~ 1*ds_4

### add autogregressive effects and cross-lagged effects between within centred variables
w_es2 ~ c(a,a) * w_es1 + c(b,b) * w_ds1 # specify parameter to be the same for both groups
w_ds2 ~ c(c,c) * w_es1 + c(d,d) * w_ds1 

w_es3 ~ c(e,e) * w_es2 + c(f,f) * w_ds2  
w_ds3 ~ c(g,g) * w_es2 + c(h,h) * w_ds2 

w_es4 ~ c(i,i) * w_es3 + c(j,j) * w_ds3  
w_ds4 ~ c(k,k) * w_es3 + c(l,l) * w_ds3 


# Estimate covariance between within-person centered variables at first wave
  w_es1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_es2 ~~ w_ds2
  w_es3 ~~ w_ds3
  w_es4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_es ~~ RI_es
  RI_ds ~~ RI_ds
  RI_es ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_es1 ~~ w_es1 # Variances
  w_ds1 ~~ w_ds1 
  w_es2 ~~ w_es2 # Residual variances
  w_ds2 ~~ w_ds2 
  w_es3 ~~ w_es3 
  w_ds3 ~~ w_ds3 
  w_es4 ~~ w_es4 
  w_ds4 ~~ w_ds4 
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 

RICLPM_es_constrained_fit <- blavaan(RICLPM_es_constrained, 
                                     data = data,
                                     target = "stan",
                                     bcontrol = list(cores = 4),
                                     n.chains = 4,
                                     sample = 10000,
                                     group = "treat", 
                                     meanstructure = T, 
                                     int.ov.free = T,
                                     control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20),
                                     dp = mydp
)

save(RICLPM_es_constrained_fit, file = "RICLPM_es_constrained_fit.RData")

# extract draws from posterior
draws_es_constrained <- blavInspect(RICLPM_es_constrained_fit, "mcmc")
draws_es_constrained <- do.call("rbind", draws_es_constrained)

# check traceplots 
mcmc_obj2 <- blavInspect(RICLPM_es_constrained_fit, "mcobj")

pdf(file="es_traceplots_constrained.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj2, pars = colnames(draws_es_constrained)[grep("bet_sign", colnames(draws_es_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_es_constrained)[grep("Psi_cov", colnames(draws_es_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_es_constrained)[grep("Psi_var", colnames(draws_es_constrained))])
mcmc_trace(mcmc_obj2, pars = colnames(draws_es_constrained)[grep("Nu_free", colnames(draws_es_constrained))])
dev.off()

summary(RICLPM_es_constrained_fit)

# plot posterior distributions
par(mfrow=c(10,7), mar=c(2,0,2,0))
for(i in 1:70){
  plot(density(draws_es_constrained[,i]), xlab="", ylab="", main=colnames(draws_es_constrained)[i], yaxt="n")
}

################### MODEL 3: TREATMENT AS PREDICTOR ############################

## model with treatment added as predcitor
RICLPM_es_treat <- '
  # Create between components (random intercepts)
  RI_es =~ 1*es_1 + 1*es_2 + 1*es_3 + 1*es_4 ## emotion suppression
  RI_ds =~ 1*ds_1 + 1*ds_2 + 1*ds_3 + 1*ds_4 ## depressive symptoms

  # Create within-person centered variables
  w_es1 =~ 1*es_1
  w_es2 =~ 1*es_2
  w_es3 =~ 1*es_3
  w_es4 =~ 1*es_4
  w_ds1 =~ 1*ds_1
  w_ds2 =~ 1*ds_2
  w_ds3 =~ 1*ds_3
  w_ds4 =~ 1*ds_4

  ### add autogregressive effects and cross-lagged effects between within centred variables
  w_es2 ~ w_es1 + w_ds1 + treat # here we add treatment as a predictor
  w_ds2 ~ w_es1 + w_ds1 + treat

  w_es3 ~ w_es2 + w_ds2 + treat
  w_ds3 ~ w_es2 + w_ds2 + treat

  w_es4 ~ w_es3 + w_ds3 + treat
  w_ds4 ~ w_es3 + w_ds3 + treat

  # Estimate covariance between within-person centered variables at first wave
  w_es1 ~~ w_ds1 
  
  # Estimate covariances between residuals of within-person centered variables 
  # (i.e., innovations)
  w_es2 ~~ w_ds2
  w_es3 ~~ w_ds3
  w_es4 ~~ w_ds4
  
  # Estimate variance and covariance of random intercepts 
  RI_es ~~ RI_es
  RI_ds ~~ RI_ds
  RI_es ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_es1 ~~ w_es1 # Variances
  w_ds1 ~~ w_ds1 
  w_es2 ~~ w_es2 # Residual variances
  w_ds2 ~~ w_ds2 
  w_es3 ~~ w_es3 
  w_ds3 ~~ w_ds3 
  w_es4 ~~ w_es4 
  w_ds4 ~~ w_ds4 
'

mydp <- dpriors(theta="gamma(.1,.1)[sd]", beta="normal(0,3)", psi="gamma(.1,.1)") 

# fit model with treatment as predictor
RICLPM_es_treat_fit <- blavaan(RICLPM_es_treat, 
                               data = data,
                               target = "stan",
                               bcontrol = list(cores = 4),
                               n.chains = 4,
                               sample = 10000,
                               meanstructure = T, 
                               int.ov.free = T,
                               control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20),
                               dp = mydp)

save(RICLPM_es_treat_fit, file = "RICLPM_es_treat_fit.RData")

# extract draws from posterior
draws_es_treat_fit <- blavInspect(RICLPM_es_treat_fit, "mcmc")
draws_es_treat_fit <- do.call("rbind", draws_es_treat_fit)

# check traceplots 
mcmc_obj3 <- blavInspect(RICLPM_es_treat_fit, "mcobj")

pdf(file="es_traceplots_treat.pdf", width = 17.72, height = 11.24)
mcmc_trace(mcmc_obj3, pars = colnames(draws_es_treat_fit)[grep("bet_sign", colnames(draws_es_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_es_treat_fit)[grep("Psi_cov", colnames(draws_es_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_es_treat_fit)[grep("Psi_var", colnames(draws_es_treat_fit))])
mcmc_trace(mcmc_obj3, pars = colnames(draws_es_treat_fit)[grep("Nu_free", colnames(draws_es_treat_fit))])
dev.off()

summary(RICLPM_es_constrained_fit)

# plot posterior distributions
par(mfrow=c(7,6), mar=c(2,0,2,0))
for(i in 1:41){
  plot(density(draws_es_treat_fit[,i]), xlab="", ylab="", main=colnames(draws_es_treat_fit)[i], yaxt="n")
}

