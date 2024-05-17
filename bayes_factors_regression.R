# ==============================================================================
# BAYES FACTORS FOR REGRESSION COEFFICIENTS
# ==============================================================================

###################### EMOTION SUPPRESSION #####################################

regressions_es <- data.frame(pxnames = RICLPM_es_treat_fit@ParTable$pxnames,
                             est = RICLPM_es_treat_fit@ParTable$est,
                             path = paste(RICLPM_es_treat_fit@ParTable$lhs, 
                                          RICLPM_es_treat_fit@ParTable$op, 
                                          RICLPM_es_treat_fit@ParTable$rhs))[17:34,]
regressions_es$BF10 <- rep(NA, 18)

mcmc_obj_es <- blavInspect(RICLPM_es_treat_fit, "mcobj")
draws_es <- as.matrix(mcmc_obj_es)

for(i in 1:18){
  es_posterior <- logspline(draws_es[,regressions_es$pxnames[i]])
  regressions_es$BF10[i] <- dnorm(0, 0, 3) / dlogspline(0, es_posterior)
}

regressions_es$BF10 <- round(regressions_es$BF10, 3)

regressions_es$BF_interpretation <- ifelse(regressions_es$BF10 < 1/10,
                                           "strong evidence for H0",
                                           ifelse(regressions_es$BF10 > 10, 
                                                  "strong evidence for H1",
                                                  ifelse(regressions_es$BF10 < 1/3,
                                                         "moderate evidence for H0",
                                                         ifelse(regressions_es$BF10 > 3,
                                                                "moderate evidence for H1",
                                                                "inconclusive evidence"))))

###################### COGNITIVE REAPPRAISAL ###################################

regressions_cr <- data.frame(pxnames = RICLPM_cr_treat_fit@ParTable$pxnames,
                             est = RICLPM_cr_treat_fit@ParTable$est,
                             path = paste(RICLPM_cr_treat_fit@ParTable$lhs, 
                                          RICLPM_cr_treat_fit@ParTable$op, 
                                          RICLPM_cr_treat_fit@ParTable$rhs))[17:34,]
regressions_cr$BF10 <- rep(NA, 18)

mcmc_obj_cr <- blavInspect(RICLPM_cr_treat_fit, "mcobj")
draws_cr <- as.matrix(mcmc_obj_cr)

for(i in 1:18){
  cr_posterior <- logspline(draws_cr[,regressions_cr$pxnames[i]])
  regressions_cr$BF10[i] <- dnorm(0, 0, 3) / dlogspline(0, cr_posterior)
}

regressions_cr$BF10 <- round(regressions_cr$BF10, 3)

regressions_cr$BF_interpretation <- ifelse(regressions_cr$BF10 < 1/10,
                                           "strong evidence for H0",
                                           ifelse(regressions_cr$BF10 > 10, 
                                                  "strong evidence for H1",
                                                  ifelse(regressions_cr$BF10 < 1/3,
                                                         "moderate evidence for H0",
                                                         ifelse(regressions_cr$BF10 > 3,
                                                                "moderate evidence for H1",
                                                                "inconclusive evidence"))))

###################### ATTACHMENT AVOIDANCE ####################################

regressions_av <- data.frame(pxnames = RICLPM_av_treat_fit@ParTable$pxnames,
                             est = RICLPM_av_treat_fit@ParTable$est,
                             path = paste(RICLPM_av_treat_fit@ParTable$lhs, 
                                          RICLPM_av_treat_fit@ParTable$op, 
                                          RICLPM_av_treat_fit@ParTable$rhs))[13:24,]
regressions_av$BF10 <- rep(NA, 12)

mcmc_obj_av <- blavInspect(RICLPM_av_treat_fit, "mcobj")
draws_av <- as.matrix(mcmc_obj_av)

for(i in 1:12){
  av_posterior <- logspline(draws_av[,regressions_av$pxnames[i]])
  regressions_av$BF10[i] <- dnorm(0, 0, 3) / dlogspline(0, av_posterior)
}

regressions_av$BF10 <- round(regressions_av$BF10, 3)

regressions_av$BF_interpretation <- ifelse(regressions_av$BF10 < 1/10,
                                           "strong evidence for H0",
                                           ifelse(regressions_av$BF10 > 10, 
                                                  "strong evidence for H1",
                                                  ifelse(regressions_av$BF10 < 1/3,
                                                         "moderate evidence for H0",
                                                         ifelse(regressions_av$BF10 > 3,
                                                                "moderate evidence for H1",
                                                                "inconclusive evidence"))))

######################## ATTACHMENT ANXIETY ####################################

regressions_ax <- data.frame(pxnames = RICLPM_ax_treat_fit@ParTable$pxnames,
                             est = RICLPM_ax_treat_fit@ParTable$est,
                             path = paste(RICLPM_ax_treat_fit@ParTable$lhs, 
                                          RICLPM_ax_treat_fit@ParTable$op, 
                                          RICLPM_ax_treat_fit@ParTable$rhs))[13:24,]
regressions_ax$BF10 <- rep(NA, 12)

mcmc_obj_ax <- blavInspect(RICLPM_ax_treat_fit, "mcobj")
draws_ax <- as.matrix(mcmc_obj_ax)

for(i in 1:12){
  ax_posterior <- logspline(draws_ax[,regressions_ax$pxnames[i]])
  regressions_ax$BF10[i] <- dnorm(0, 0, 3) / dlogspline(0, ax_posterior)
}

regressions_ax$BF10 <- round(regressions_ax$BF10, 3)

regressions_ax$BF_interpretation <- ifelse(regressions_ax$BF10 < 1/10,
                                           "strong evidence for H0",
                                           ifelse(regressions_ax$BF10 > 10, 
                                                  "strong evidence for H1",
                                                  ifelse(regressions_ax$BF10 < 1/3,
                                                         "moderate evidence for H0",
                                                         ifelse(regressions_ax$BF10 > 3,
                                                                "moderate evidence for H1",
                                                                "inconclusive evidence"))))

