# Author: Teagan Pearse
# Date: 01/03/2026

# Purpose: Investigate treatment effects on liver glycogen in rats, accounting for nested liver samples using a linear mixed-effects model.

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
rat_specific_glycogen_pred <- ggpredict(rats_correct, 
                              terms = c("Treatment", "Rat"),
                              type = "random") |>
  plot(show_data = TRUE) +
  labs(title = "Rat-specific predictions of liver glycogen across treatments",
       x = "Treatment",
       y = "Glycogen")

print(rat_specific_glycogen_pred)

ggsave(here("figures", "rat_specific_glycogen_pred.pdf"),
       plot = rat_specific_glycogen_pred, width = 20, height = 15, units = "cm")

# ---- Write up ----
# A linear mixed-effects model was fitted using REML in the lmerTest package to investigate treatment effects on liver glycogen, accounting for nested sampling of liver sections within rats. The model included 36 observations from 6 rats and 18 liver samples, with random intercepts for rats and liver samples nested within rats (1|Rat/Liver). Degrees of freedom were estimated using Satterthwaite's method. Treatment 2 (β = 10.50, 95% CI: -6.17-27.17, p = 0.213) and Treatment 3 (β = -5.33, 95% CI = -22.00-11.34, p = 0.482) did not significantly differ from Treatment 1. Variance components were 36.06 (rat), 14.17 (liver within rat), and 21.17 (residual). Diagnostics indicated no major assumption violations.

