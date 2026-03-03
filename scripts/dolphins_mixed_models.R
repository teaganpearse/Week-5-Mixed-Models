# Author: Teagan Pearse
# Date: 28/02/2026

# Purpose: Investigate how body mass and breath direction influence tidal volume in dolphins, accounting for repeated measures using a linear mixed-effects model.

# ---- Setup ----
library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(MuMIn)
library(performance)
library(sjPlot)

# ---- Load and prepare data ----
dolphins <- read_csv(here("data", "dolphins.csv"))

dolphins <- dolphins %>%
  mutate(direction = factor(direction)) |>
  drop_na()

# Fixed vs random effects:
  # bodymass = fixed; direction = fixed; animal = random (repeated measures)

# ---- Fit mixed model ----
dolphmod <- lmer(vt ~ bodymass + direction + (1|animal), data = dolphins)
summary(dolphmod)
  # Fixed effects: bodymass = 0.017 (SE = 0.003, p < 0.001); direction (inhalation vs exhalation) = 1.11 (SE = 0.20, p < 0.001)
  # Random effects: animal variance = 1.04 (SD = 1.02); residual variance = 1.16 (SD = 1.08); ICC ≈ 0.47 (≈ 47% of variance attributable to differences between dolphins)
  # Tidal volume increases with body mass and is higher during inhalation

# ---- Model fit and diagnostics ----
r.squaredGLMM(dolphmod)
  # R²m = 0.41; R²c = 0.69
  # Fixed effects explain ~41% of variance; full model explains ~ 69%
  # Random effects substantially increase explained variance

check_model(dolphmod, detrend = FALSE)
  # Diagnostics indicate no major assumption violations -> appropriate model fit

# ---- Predictions and plots ----

# Predicted tidal volume by direction at mean body mass
emmeans(dolphmod, specs = "direction")
  # mean tidal volume higher during inhalation (6.37 L, 95% CI: 5.90-6.84) than exhalation (5.26 L, 95% CI: 4.78-5.73)

# Population-average tidal volume vs bodymass
dolphins_bodymass_population_pred <- emmeans(
  dolphmod, specs = "bodymass",
  at = list(bodymass = seq(130, 260, 10))
) %>%
  as.data.frame() %>%
  ggplot(aes(x = bodymass, y = emmean)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
              fill = "steelblue",
              alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Body mass (kg)", y = "Tidal volume (L)") +
  theme_bw()+
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "dolphins_bodymass_population_pred.pdf"),
       plot = dolphins_bodymass_population_pred, width = 20, height = 15, units = "cm")

# Tidal volume vs bodymass by direction
dolphins_bodymass_direction_pred <- ggpredict(dolphmod, terms = c("bodymass", "direction")) %>%
  plot(show_data = TRUE) +
  labs(title = "Predicted values of tidal volume",
       x = "Body mass (kg)",
       y = "Tidal volume (L)") +
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here("figures", "dolphins_bodymass_direction_pred.pdf"),
       plot = dolphins_bodymass_direction_pred, width = 20, height = 15, units = "cm")

# ---- Reporting table ----
tab_model(dolphmod,
          show.re.var = TRUE, 
          show.icc = TRUE, 
          show.r2 = TRUE)
  # Formatted results table (fixed effects (CI), variance components, ICC, R²)

# ---- Write up ----
# A linear mixed-effects model was fitted using REML in the lmerTest package to investigate how body mass and breath direction influence tidal volume in dolphins. Dolphin ID was included as a random intercep to account for repeated measures, and degrees of freedom were estimated using Satterthwaite's method. The analysis included 112 observations from 31 dolphins (1-4 measurements per individual). Tidal volume increased by 0.017 L per kg body mass (95% CI: 0.008-0.024, p < 0.001) and was 1.11 L higher during inhalation than exhalation (95% CI: 0.71-1.52, p < 0.001). Between-dolphin variance (σ² = 1.04) was similar to residual variance (σ² = 1.16); ICC ≈ 0.47 (47%). Fixed effects explained 41% of variance (R²m = 0.41), while the full model explained 69% (R²c = 0.69). Diagnostic checks indicated no major assumption violations.
