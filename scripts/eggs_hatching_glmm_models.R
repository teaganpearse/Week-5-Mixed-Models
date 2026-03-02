# Author: Teagan Pearse
# Date: 02/03/2026

# Purpose: Investigate hatching success across genetic crosses, addressing overdispersion using binomial generalised linear mixed models

# ---- Setup ----
library(tidyverse)
library(here)
library(lme4)
library(performance)
library(ggeffects)

# ---- Load and prepare data ----
egg_hatch_data <- read_csv(here("data", "egg_hatch_mixed.csv"))

egg_hatch_data <- egg_hatch_data |>
  mutate(unhatched_eggs = number_of_eggs - number_of_larvae)

# ---- Standard binomial GLM ----
binomial_model <- glm(
  cbind(number_of_larvae, unhatched_eggs) ~ cross,
  family = binomial,
  data = egg_hatch_data)

summary(binomial_model)
  # Fixed effects: cross A (intercept) = 1.05 (SE 0.05, p < 0.001); cross B = 0.39 (SE = 0.91, p < 0.001); cross C = 1.29 (SE = 0.12, p < 0.001)
  # Null deviance = 1307.6; residual deviance = 1157.3; AIC = 1363.6

# ---- Check overdispersion ----
performance::check_overdispersion(binomial_model)
  # overdispersion detected  (1.703; p < 0.001)

check_model(binomial_model)
  # Diagnostic plots suggest poor model fit and overdispersion
  # Residual patterns and influential points detected
  # GLM assumptions violated

# ---- Fit GLMMs to address overdispersion ----

# Random intercept for replicate only
binomial_mixed_rep <-glmer(
  cbind(number_of_larvae, unhatched_eggs) ~ cross + (1|replicate),
  family = binomial,
  data = egg_hatch_data)

summary(binomial_mixed_rep)
  # Fixed effects: cross A (intercept) = 1.09 (SE = 0.19; p < 0.001); cross B = 0.61 (SE = 0.10; p < 0.001); cross C = 1.04 (SE = 0.12; p < 0.001)
  # Random effects: replicate variance = 0.16 (SD = 0.40); 5 replicate groups
  # Model fit improved relative to GLM (AIC = 1263.2)

# Random intercepts for replicate and female (nested)
binomial_mixed_rep_fem <- glmer(
  cbind(number_of_larvae, unhatched_eggs) ~ cross + (1|replicate/female_num),
  family = binomial,
  data = egg_hatch_data)

summary(binomial_mixed_rep_fem)
  # Fixed effects: cross A (intercept) = 1.25 (SE = 0.24; p < 0.001); cross B = 0.93 (SE = 0.13; p < 0.001); cross C = 0.70 (SE = 0.15; p < 0.001)
  # Random effects: female_num:replicate variance = 1.72 (SD = 1.31); replicate variance = 0.01 (SD = 0.12)
  # Most variation occurs between females; minimal variation between replicates
# Model fit substantially improved (AIC = 644.1)

# ---- Compare models ----
AIC(binomial_model, binomial_mixed_rep, binomial_mixed_rep_fem)

  # GLM = 1363.6; GLMM (replicate) = 1263.2; GLMM (replicate/female) = 644.1
  # Nested random-effects model provides best fit

# ---- Re-check overdispersion ----
performance::check_overdispersion(binomial_mixed_rep_fem)
  # No overdispersion detected (1.153; p = 0.224)

# ---- Examine final model ----
summary(binomial_mixed_rep_fem)
  # Most variation is between females (SD = 1.31); little variation between replicates (SD = 0.12)
  # Nested random effects required because females are unique within replicates
  # Groups: female:replicate = 34; replicate = 5

# ---- Visualise results ----
egg_cross_population_pred <- ggpredict(binomial_mixed_rep_fem, terms = "cross") |>
  plot()

ggsave(here("figures", "egg_cross_population_pred.pdf"),
       plot = egg_cross_population_pred, width = 20, height = 15, units = "cm")
