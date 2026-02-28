# Author: Teagan Pearse
# Date: 27/02/2026

# Purpose: Assess the effect of toxin concentration on larval detoxification capacity, accounting for batch-level variation using a linear mixed-effects model.

# ---- Setup ----
library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(MuMIn)
library(performance)
library(sjPlot)

# ----Load and prepare data ----
benzo <- read_csv(here("data", "benzo.csv"))

benzo <- benzo %>%
  mutate(group = factor(group,
                        levels = 1:5,
                        labels = c("a", "b", "c", "d", "e")))

# ---- Exploratory plots ----

# Detoxification vs toxin concentration (pooled data, batch ignored)
benzo_pooled <- ggplot(benzo, aes(x = benzo_um, y = detox_exp)) +
  geom_point(alpha = 0.6) +
  labs(
    x = "Toxin concentration (µM)",
    y = "Detoxification capacity"
  ) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "benzo_pooled.pdf"),
       plot = benzo_pooled, width = 20, height = 15, units = "cm")

# Detoxification vs toxin concentration by batch
benzo_by_batch <- ggplot(benzo, aes(x = benzo_um, y = detox_exp, colour = group)) +
  geom_point(alpha = 0.6) +
  labs(
    x = "Toxin concentration (µM)",
    y = "Detoxification capacity"
  ) +
  theme(legend.position = "right") +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "benzo_by_batch.pdf"),
       plot = benzo_by_batch, width = 20, height = 15, units = "cm")

# Distribution of detoxification capacity across batches
benzo_by_batch_boxplot <- ggplot(benzo, aes(x = group, y = detox_exp, fill = group)) +
  geom_boxplot() +
  labs(
    x = "Experimental batch",
    y = "Detoxification capacity"
  ) +
  theme(legend.position = "none") +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "benzo_by_batch_boxplot.pdf"),
       plot = benzo_boxplot, width = 20, height = 15, units = "cm")

# ---- Fit models ----
  
# Pooled model (ignores grouping)
pooled_model <- lm(detox_exp ~ benzo_um, data = benzo)
summary(pooled_model)
  # intercept: estimate = 21.05, SE = 1.65, p < 0.001
  # benzo_um: estimate = 2.58, SE = 0.29, p < 0.001

# Fixed effects model (batch as a fixed factor)
fixed_model <- lm(detox_exp ~ benzo_um + group, data = benzo)
summary(fixed_model)
  # benzo_um: estimate = 2.02, SE = 0.17, p < 0.001
  # Significant differences between batches

# Mixed effects model (batch as random intercept)
mixed_model <- lmer(detox_exp ~ benzo_um + (1 | group), data = benzo)
summary(mixed_model)
  # Random effects: group variance = 205.0 (SD = 14.32), residual variance = 101.0 (SD = 10.05)
  # Fixed effects: benzo_um = 2.03, SE = 0.17, p < 0.001
  # Positive toxin effect shown (detox increases 2.03 units per µM toxin)

# ---- Predictions and visualisations ----
  
# Population predictions
emmeans(mixed_model,
          specs = "benzo_um",
          at = list(benzo_um = c(0, 2.5, 5, 7.5, 10)))
  # Detoxification increases ~5 units per 2.5 µM step (23.3 -> 28.3 -> 33.4 -> 38.5 -> 43.5)

# Population average prediction curve
benzo_population_pred <- emmeans(mixed_model,
                                   specs = "benzo_um",
                                   at = list(benzo_um = seq(0, 10, 0.5))) |>
  as.data.frame() |>
  ggplot(aes(x = benzo_um, y = emmean)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, alpha = 0.3)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Toxin concentration (µM)",
    y = "Detoxification capacity",
    title = "Population-average predictions") +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "benzo_population_pred.pdf"),
       plot = benzo_population_pred, width = 20, height = 15, units = "cm")

# Population-level predictions
benzo_population_points <- ggpredict(mixed_model,
                                       terms = "benzo_um",
                                       type = "fixed") |>
  plot(show_data = TRUE) +
  labs(x = "Toxin concentration (µM)",
       y = "Detoxification capacity") +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "benzo_population_points.pdf"),
       plot = benzo_population_points, width = 20, height = 15, units = "cm")

# Group-specific predictions
benzo_group_predictions <- ggpredict(mixed_model,
                                       terms = c("benzo_um", "group"),
                                       type = "random") |>
  plot(show_data = TRUE) +
  facet_wrap(~group) +
  labs(x = "Toxin concentration (µM)",
       y = "Detoxification capacity")

ggsave(here("figures", "benzo_group_predictions.pdf"),
       plot = benzo_group_predictions, width = 20, height = 15, units = "cm")

# Final group predictions
benzo_group_final <- ggpredict(mixed_model,
                                 terms = c("benzo_um", "group"),
                                 type = "random") |>
  
  plot(show_data = TRUE) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~group) +
  labs(
    x = "Toxin concentration (µM)",
    y = "Detoxification capacity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(here("figures", "benzo_group_final.pdf"),
       plot = benzo_group_final, width = 20, height = 15, units = "cm")

# ---- Model fit and diagnostics ----
  
# R² for mixed model
r.squaredGLMM(mixed_model)
  # R²m = 0.10 (toxin only); R²c = 0.70 (including batch random effects)
  # Fixed effects explain ~10% of variance; full mixed model explains ~70%
  # Random effects add considerable additional variance and improve overall model fit

# Model diagnostics
check_model(mixed_model, detrend = FALSE)
  # diagnostics show acceptable fit -> indicate no major violations of linearity, homoscedasticity, or normality assumptions

# ---- Reporting table ----
tab_model(mixed_model,
            show.re.var = TRUE,
            show.icc = TRUE,
            show.r2 = TRUE)
  # mixed model results table including fixed effects, random effects, and model fit

# ---- Write up ----
# A linear mixed-effects model was fitted using REML estimation in the lmerTest package, to investigate how dietary benzoxazinoid concentration influences larval detoxification capacity. Experimental batch was included as a random intercept. Detoxification capacity increased by 2.03 units per µM toxin (95% CI: 1.69-2.36, p < 0.001). Between-batch variance (σ² = 205) exceeded residual variance (σ² = 101). The model included 430 observations across 5 batches. Fixed effects explained 10% of variance (R²m = 0.10), while the full model explained 70% (R²c = 0.70). Diagnostic checks indicated no major violations of model assumptions.