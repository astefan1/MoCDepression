# ==============================================================================
# ESTIMATING AND TESTING THE INDIRECT PATHS
# ==============================================================================

###################### EMOTION SUPPRESSION #####################################

# Extract Stan object

mcmc_obj_es <- blavInspect(RICLPM_es_treat_fit, "mcobj")

# Figure out which internal parameter belongs to which displayed connection
# 
# data.frame(pxnames = RICLPM_es_treat_fit@ParTable$pxnames,
#            est = RICLPM_es_treat_fit@ParTable$est,
#            connection = paste(RICLPM_es_treat_fit@ParTable$lhs, RICLPM_es_treat_fit@ParTable$op, RICLPM_es_treat_fit@ParTable$rhs))

# Extract parameter draws

draws_es <- as.matrix(mcmc_obj_es)

# Indirect effect treat - w_es2 - w_ds3
a1a2_es <- draws_es[,"bet_sign[13]"] * draws_es[,"bet_sign[4]"]

# Indirect effect treat - w_es3 - w_ds4
b1b2_es <- draws_es[,"bet_sign[14]"] * draws_es[,"bet_sign[6]"]

# Indirect effect treat - w_ds2 - w_es3
c1c2_es <- draws_es[,"bet_sign[16]"] * draws_es[,"bet_sign[9]"]

# Indirect effect treat - w_ds3 - w_es4
d1d2_es <- draws_es[,"bet_sign[17]"] * draws_es[,"bet_sign[11]"]

# Posterior mean, sd, 95%CI
post_indirect_es <- matrix(NA, nrow=4, ncol=4)
colnames(post_indirect_es) <- c("Estimate", "Post.SD", "pi.lower", "pi.upper")
rownames(post_indirect_es) <- c("a1a2", "b1b2", "c1c2", "d1d2")
post_indirect_es[,1] <- c(mean(a1a2_es), mean(b1b2_es), mean(c1c2_es), mean(d1d2_es))
post_indirect_es[,2] <- c(sd(a1a2_es), sd(b1b2_es), sd(c1c2_es), sd(d1d2_es))
post_indirect_es[,3:4] <- matrix(c(quantile(a1a2_es, probs = c(0.025, 0.975)), 
                                   quantile(b1b2_es, probs = c(0.025, 0.975)), 
                                   quantile(c1c2_es, probs = c(0.025, 0.975)), 
                                   quantile(d1d2_es, probs = c(0.025, 0.975))), 
                                 byrow=TRUE, nrow=4)

# Savage-Dickey-Density Ratio a1a2
set.seed(1234)
a1a2_es_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
a1a2_es_posterior <- logspline(a1a2_es)
BF01_a1a2_es <- dlogspline(0, a1a2_es_posterior)/ dlogspline(0, a1a2_es_prior)

# Savage-Dickey-Density Ratio b1b2
set.seed(1234)
b1b2_es_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
b1b2_es_posterior <- logspline(b1b2_es)
BF01_b1b2_es <- dlogspline(0, b1b2_es_posterior)/ dlogspline(0, b1b2_es_prior)

# Savage-Dickey-Density Ratio c1c2
set.seed(1234)
c1c2_es_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
c1c2_es_posterior <- logspline(c1c2_es)
BF01_c1c2_es <- dlogspline(0, c1c2_es_posterior)/ dlogspline(0, c1c2_es_prior)

# Savage-Dickey-Density Ratio d1d2
set.seed(1234)
d1d2_es_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
d1d2_es_posterior <- logspline(d1d2_es)
BF01_d1d2_es <- dlogspline(0, d1d2_es_posterior)/ dlogspline(0, d1d2_es_prior)

BF01_es <- cbind(BF01_a1a2_es, BF01_b1b2_es, BF01_c1c2_es, BF01_d1d2_es)

###################### COGNITIVE REAPPRAISAL ###################################

# Extract Stan object

mcmc_obj_cr <- blavInspect(RICLPM_cr_treat_fit, "mcobj")

# Figure out which internal parameter belongs to which displayed connection
# 
# data.frame(pxnames = RICLPM_cr_treat_fit@ParTable$pxnames,
#            est = RICLPM_cr_treat_fit@ParTable$est,
#            connection = paste(RICLPM_cr_treat_fit@ParTable$lhs, RICLPM_cr_treat_fit@ParTable$op, RICLPM_cr_treat_fit@ParTable$rhs))

# Extract parameter draws

draws_cr <- as.matrix(mcmc_obj_cr)

# Indirect effect treat - w_cr2 - w_ds3
a1a2_cr <- draws_cr[,"bet_sign[13]"] * draws_cr[,"bet_sign[4]"]

# Indirect effect treat - w_cr3 - w_ds4
b1b2_cr <- draws_cr[,"bet_sign[14]"] * draws_cr[,"bet_sign[6]"]

# Indirect effect treat - w_ds2 - w_cr3
c1c2_cr <- draws_cr[,"bet_sign[16]"] * draws_cr[,"bet_sign[9]"]

# Indirect effect treat - w_ds3 - w_cr4
d1d2_cr <- draws_cr[,"bet_sign[17]"] * draws_cr[,"bet_sign[11]"]

# Posterior mean, sd, 95%CI
post_indirect_cr <- matrix(NA, nrow=4, ncol=4)
colnames(post_indirect_cr) <- c("Estimate", "Post.SD", "pi.lower", "pi.upper")
rownames(post_indirect_cr) <- c("a1a2", "b1b2", "c1c2", "d1d2")
post_indirect_cr[,1] <- c(mean(a1a2_cr), mean(b1b2_cr), mean(c1c2_cr), mean(d1d2_cr))
post_indirect_cr[,2] <- c(sd(a1a2_cr), sd(b1b2_cr), sd(c1c2_cr), sd(d1d2_cr))
post_indirect_cr[,3:4] <- matrix(c(quantile(a1a2_cr, probs = c(0.025, 0.975)), 
                                   quantile(b1b2_cr, probs = c(0.025, 0.975)), 
                                   quantile(c1c2_cr, probs = c(0.025, 0.975)), 
                                   quantile(d1d2_cr, probs = c(0.025, 0.975))), 
                                 byrow=TRUE, nrow=4)

# Savage-Dickey-Density Ratio a1a2
set.seed(1234)
a1a2_cr_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
a1a2_cr_posterior <- logspline(a1a2_cr)
BF01_a1a2_cr <- dlogspline(0, a1a2_cr_posterior)/ dlogspline(0, a1a2_cr_prior)

# Savage-Dickey-Density Ratio b1b2
set.seed(1234)
b1b2_cr_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
b1b2_cr_posterior <- logspline(b1b2_cr)
BF01_b1b2_cr <- dlogspline(0, b1b2_cr_posterior)/ dlogspline(0, b1b2_cr_prior)

# Savage-Dickey-Density Ratio c1c2
set.seed(1234)
c1c2_cr_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
c1c2_cr_posterior <- logspline(c1c2_cr)
BF01_c1c2_cr <- dlogspline(0, c1c2_cr_posterior)/ dlogspline(0, c1c2_cr_prior)

# Savage-Dickey-Density Ratio d1d2
set.seed(1234)
d1d2_cr_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
d1d2_cr_posterior <- logspline(d1d2_cr)
BF01_d1d2_cr <- dlogspline(0, d1d2_cr_posterior)/ dlogspline(0, d1d2_cr_prior)

BF01_cr <- cbind(BF01_a1a2_cr, BF01_b1b2_cr, BF01_c1c2_cr, BF01_d1d2_cr)

###################### ATTACHMENT AVOIDANCE ####################################

# Extract Stan object

mcmc_obj_av <- blavInspect(RICLPM_av_treat_fit, "mcobj")

# Figure out which internal parameter belongs to which displayed connection
# 
# data.frame(pxnames = RICLPM_av_treat_fit@ParTable$pxnames,
#            est = RICLPM_av_treat_fit@ParTable$est,
#            connection = paste(RICLPM_av_treat_fit@ParTable$lhs, RICLPM_av_treat_fit@ParTable$op, RICLPM_av_treat_fit@ParTable$rhs))

# Extract parameter draws

draws_av <- as.matrix(mcmc_obj_av)

# Indirect effect treat - w_av3 - w_ds4
a1a2_av <- draws_av[,"bet_sign[9]"] * draws_av[,"bet_sign[4]"]

# Indirect effect treat - w_ds3 - w_av4
c1c2_av <- draws_av[,"bet_sign[11]"] * draws_av[,"bet_sign[7]"]

# Posterior mean, sd, 95%CI
post_indirect_av <- matrix(NA, nrow=2, ncol=4)
colnames(post_indirect_av) <- c("Estimate", "Post.SD", "pi.lower", "pi.upper")
rownames(post_indirect_av) <- c("a1a2", "c1c2")
post_indirect_av[,1] <- c(mean(a1a2_av), mean(c1c2_av))
post_indirect_av[,2] <- c(sd(a1a2_av), sd(c1c2_av))
post_indirect_av[,3:4] <- matrix(c(quantile(a1a2_av, probs = c(0.025, 0.975)), 
                                   quantile(c1c2_av, probs = c(0.025, 0.975))), 
                                 byrow=TRUE, nrow=4)

# Savage-Dickey-Density Ratio a1a2
set.seed(1234)
a1a2_av_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
a1a2_av_posterior <- logspline(a1a2_av)
BF01_a1a2_av <- dlogspline(0, a1a2_av_posterior)/ dlogspline(0, a1a2_av_prior)

# Savage-Dickey-Density Ratio c1c2
set.seed(1234)
c1c2_av_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
c1c2_av_posterior <- logspline(c1c2_av)
BF01_c1c2_av <- dlogspline(0, c1c2_av_posterior)/ dlogspline(0, c1c2_av_prior)

BF01_av <- cbind(BF01_a1a2_av, BF01_c1c2_av)

######################## ATTACHMENT ANXIETY ####################################

# Extract Stan object

mcmc_obj_ax <- blavInspect(RICLPM_ax_treat_fit, "mcobj")

# Figure out which internal parameter belongs to which displayed connection
# 
# data.frame(pxnames = RICLPM_ax_treat_fit@ParTable$pxnames,
#            est = RICLPM_ax_treat_fit@ParTable$est,
#            connection = paste(RICLPM_ax_treat_fit@ParTable$lhs, RICLPM_ax_treat_fit@ParTable$op, RICLPM_ax_treat_fit@ParTable$rhs))

# Extract parameter draws

draws_ax <- as.matrix(mcmc_obj_ax)

# Indirect effect treat - w_ax3 - w_ds4
a1a2_ax <- draws_ax[,"bet_sign[9]"] * draws_ax[,"bet_sign[4]"]

# Indirect effect treat - w_ds3 - w_ax4
c1c2_ax <- draws_ax[,"bet_sign[11]"] * draws_ax[,"bet_sign[7]"]

# Posterior mean, sd, 95%CI
post_indirect_ax <- matrix(NA, nrow=2, ncol=4)
colnames(post_indirect_ax) <- c("Estimate", "Post.SD", "pi.lower", "pi.upper")
rownames(post_indirect_ax) <- c("a1a2", "c1c2")
post_indirect_ax[,1] <- c(mean(a1a2_ax), mean(c1c2_ax))
post_indirect_ax[,2] <- c(sd(a1a2_ax), sd(c1c2_ax))
post_indirect_ax[,3:4] <- matrix(c(quantile(a1a2_ax, probs = c(0.025, 0.975)), 
                                   quantile(c1c2_ax, probs = c(0.025, 0.975))), 
                                 byrow=TRUE, nrow=4)

# Savage-Dickey-Density Ratio a1a2
set.seed(1234)
a1a2_ax_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
a1a2_ax_posterior <- logspline(a1a2_ax)
BF01_a1a2_ax <- dlogspline(0, a1a2_ax_posterior)/ dlogspline(0, a1a2_ax_prior)

# Savage-Dickey-Density Ratio c1c2
set.seed(1234)
c1c2_ax_prior <- logspline(rnorm(10000, 0, 3)*rnorm(10000, 0, 3))
c1c2_ax_posterior <- logspline(c1c2_ax)
BF01_c1c2_ax <- dlogspline(0, c1c2_ax_posterior)/ dlogspline(0, c1c2_ax_prior)

BF01_ax <- cbind(BF01_a1a2_ax, BF01_c1c2_ax)

