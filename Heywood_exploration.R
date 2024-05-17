################### HEYWOOD CASES EXPLORATION ##################################

# complete data

data_comp <- data[!is.na(rowMeans(data[,c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")])), ]
data_comp <- data_comp[, c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4", "treat") ]

data_comp_ST <- data_comp[data_comp$treat == "ST",]
data_comp_KVT <- data_comp[data_comp$treat == "KVT",]

# Create unconstrained model
RICLPM_unconstrained <- '
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
  RI_es ~~ c(v5, v6)*RI_es
  RI_ds ~~ c(v7, v8)*RI_ds
  RI_es ~~ RI_ds
  
  # Estimate (residual) variance of within-person centered variables
  w_es1 ~~ c(v1,v2)*w_es1 # Variances
  w_ds1 ~~ c(v3,v4)*w_ds1 
  w_es2 ~~ c(v9,v10)*w_es2 # Residual variances
  w_ds2 ~~ c(v11,v12)*w_ds2 
  w_es3 ~~ c(v13,v14)*w_es3 
  w_ds3 ~~ c(v15,v16)*w_ds3 
  w_es4 ~~ c(v17,v18)*w_es4 
  w_ds4 ~~ c(v19,v20)*w_ds4 
  
  # Constraints
  v1 > 0
  v2 > 0
  v3 > 0
  v4 > 0
  v5 > 0
  v6 > 0
  v7 > 0
  v8 > 0
  v9 > 0
  v10 > 0
  v11 > 0
  v12 > 0
  v13 > 0
  v14 > 0
  v15 > 0
  v16 > 0
  v17 > 0
  v18 > 0
  v19 > 0
  v20 > 0
'

# fit unconstrained model to data (or complete data)
RICLPM_unconstrained_fit <- lavaan(RICLPM_unconstrained, 
                                   data = data,
                                   estimator = "MLR", # robust ML estimation (can handle non-normally distributed data)
                                   missing = "ML", # maximum likelihood estimation
                                   group = "treat", # treatment here
                                   meanstructure = T, # includes means from observes vars into model
                                   int.ov.free = T,
                                   optim.method = "nlminb.constr")

summary(RICLPM_unconstrained_fit, fit.measures = T)

RICLPM_unconstrained_fit$summarize(selector = "raic")

cov_KVT <- cov(data[data$treat == "KVT",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")], use = "complete.obs")
cov_ST <- cov(data[data$treat == "ST",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")], use = "complete.obs")

is.positive.definite(cov_KVT)
is.positive.definite(cov_ST)

mean_KVT <- colMeans(data[data$treat == "KVT",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")], na.rm = TRUE)
mean_ST <- colMeans(data[data$treat == "ST",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")], na.rm = TRUE)

n_KVT <- sum(complete.cases(data[data$treat == "KVT",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")]))
n_ST <- sum(complete.cases(data[data$treat == "ST",c("ds_1", "ds_2", "ds_3", "ds_4", "es_1", "es_2", "es_3", "es_4")]))

# fit model again based on mean vectors and covariance matrices

RICLPM_unconstrained_fit <- lavaan(RICLPM_unconstrained,
                                   sample.cov = list(cov_KVT, cov_ST),
                                   sample.mean = list(mean_KVT, mean_ST),
                                   sample.nobs = list(n_KVT, n_ST),
                                   estimator = "ULS", # had to change to indicate that it's complete cases only and regular GLS algorithm didn't converge
                                   meanstructure = T, # includes means from observes vars into model
                                   int.ov.free = T,
                                   optim.method = "nlminb.constr")

summary(RICLPM_unconstrained_fit, fit.measures = T)

ycov(data_comp_KVT[,-9])

apply(data_comp_ST[,-9], 2, sd)
apply(data_comp_KVT[,-9], 2, sd)

apply(data_comp_ST[,-9], 2, mean)
apply(data_comp_KVT[,-9], 2, mean)
