# ==============================================================================
# MODEL COMPARISON USING LOO-CV
# ==============================================================================

########################## EMOTION SUPPRESSION #################################

loo_es_unconstrained <- loo(blavInspect(RICLPM_es_unconstrained_fit, "mcobj"))
loo_es_constrained <- loo(blavInspect(RICLPM_es_constrained_fit, "mcobj"))
loo_es_treat <- loo(blavInspect(RICLPM_es_treat_fit, "mcobj"))

loo_compare_es <- loo_compare(loo_es_unconstrained, loo_es_constrained, loo_es_treat)

######################## COGNITIVE REAPPRAISAL #################################

loo_cr_unconstrained <- loo(blavInspect(RICLPM_cr_unconstrained_fit, "mcobj"))
loo_cr_constrained <- loo(blavInspect(RICLPM_cr_constrained_fit, "mcobj"))
loo_cr_treat <- loo(blavInspect(RICLPM_cr_treat_fit, "mcobj"))

loo_compare_cr <- loo_compare(loo_cr_unconstrained, loo_cr_constrained, loo_cr_treat)

######################### ATTACHMENT AVOIDANCE #################################

loo_av_unconstrained <- loo(blavInspect(RICLPM_av_unconstrained_fit, "mcobj"))
loo_av_constrained <- loo(blavInspect(RICLPM_av_constrained_fit, "mcobj"))
loo_av_treat <- loo(blavInspect(RICLPM_av_treat_fit, "mcobj"))

loo_compare_av <- loo_compare(loo_av_unconstrained, loo_av_constrained, loo_av_treat)

########################### ATTACHMENT ANXIETY #################################

loo_ax_unconstrained <- loo(blavInspect(RICLPM_ax_unconstrained_fit, "mcobj"))
loo_ax_constrained <- loo(blavInspect(RICLPM_ax_constrained_fit, "mcobj"))
loo_ax_treat <- loo(blavInspect(RICLPM_ax_treat_fit, "mcobj"))

loo_compare_ax <- loo_compare(loo_ax_unconstrained, loo_ax_constrained, loo_ax_treat)

