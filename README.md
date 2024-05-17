# Mechanisms of Change in Depression

This projects contains a collection of R scripts necessary to reproduce the Bayesian Random Intercept Cross-Lagged Panel Models in the manuscript "Emotion Regulation and Attachment as Mechanisms of Change in Schema Therapy and Cognitive Behaviour Therapy for Depression" (Egli, Badenbach, van Emmerik, Stefan, OPTIMA Study Group, & Kopf-Beck, 2024). To reproduce the results, follow the steps below.

1. Run the RI-CLPM models using the scripts `rerun_blavaan_av.R`, `rerun_blavaan_ax.R`, `rerun_blavaan_cr.R`, and `rerun_blavaan_es.R`. The model fits will be saved in RData files that will be re-used in the following scripts.
2. The script `model_comparison.R` contains the code to run the model comparisons between the unconstrained, constrained, and treatment model and computes the loo coefficients for each comparison.
3. The pre-registered hypothesis tests can be conducted with the scripts `bayes_factors_regression.R` and `indirect_paths.R`. They yield Bayes factors for tests on the individual model parameters and their combinations.
4. Exploratory analyses regarding the development of the mechanism of change variables across the time of the study can be found in the `exploratory_MoC_development.R` script.
5. The `Heywood_exploration.R` script contains some code that we used to come to the conclusions mentioned in the manuscript about the reasons for the Heywood cases found in the original frequentist analyses. 


 
