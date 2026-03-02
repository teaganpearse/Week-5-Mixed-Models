# Author: Teagan Pearse
# Date: 01/03/2026

# Purpose: Test treatment effects on liver glycogen while accounting for nested sampling (liver samples nested within rats) using a linear mixed-effects model

# ---- Setup ----
library(tidyverse)
library(here)
library(lmerTest)
library(ggeffects)


# ---- Load and prepare data ----
rats <- readRDS(here("data", "rats.rds"))

rats <- rats %>%
  mutate(
    Rat = factor(Rat),
    Liver = factor(Liver),
    Treatment = factor(Treatment))

# ---- Examine data structure ----
rats |>
  aggregate(Glycogen ~ Rat + Treatment + Liver, data = _, mean) |>
  head(12)

# ---- Models ----

# Incorrect approach: treating effects as crossed
rats_wrong <- lmer(Glycogen ~ Treatment + (1|Rat) + (1|Liver), data = rats)
summary(rats_wrong)
  # treats liver sample "1" as comparable across all rats
  # only 3 liver groups detected (instead of 18), indicating crossed rather than nested design

# Correct approach: nested design
rats_correct <- lmer(Glycogen ~ Treatment + (1|Rat/Liver), data = rats)

summary(rats_correct)
  # fixed effects: treatment1 (intercept) = 140.500 (SE = 4.71, p <0.001); treatment2 = 10.50. (SE = 6.66, p = 0.213); treatment3 = -5.33 (SE = 6.66, p = 0.482)
  # random effects: rat variance = 36.06 (SD = 6.01); Rat:liver variance = 14.17 (SD = 3.76); residual variance = 21.17 (SD = 4.60)

# ---- Predictions and plots ----

# Treatment effect on liver glycogen levels
rats_treatment_pred <- ggpredict(rats_correct, terms = "Treatment") %>%
  plot() +
  labs(x = "Treatment", y = "Glycogen")

ggsave(here("figures", "rats_treatment_pred.pdf"),
       plot = rats_treatment_pred, width = 20, height = 15, units = "cm")

# Rat-specific predictions showing individual variation
rats_by_rat_pred <- ggpredict(rats_correct, 
                              terms = c("Treatment", "Rat"),
                              type = "random") |>
  plot(show_data = TRUE)

print(rats_by_rat_pred)

ggsave(here("figures", "rats_by_rat_pred.pdf"),
       plot = rats_by_rat_pred, width = 20, height = 15, units = "cm")
