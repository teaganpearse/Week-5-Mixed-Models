# Mixed Models - Week 5

## Project Overview

Analysis of grouped and nested experimental data using linear and generalised linear mixed models.

## Data

**Benzo**
- **Sample size**: 430 observations across 5 experimental batches
- **Response variable**: Detoxification capacity
- **Primary predictor**: Benzoxazinoid concentration (µM)
- **Random effects**: Experimental batch

**Dolphins**
- **Sample size**: 112 observations from 31 dolphins (1-4 measurements per individual)
- **Response variable**: Tidal volume (L)
- **Primary predictors**: Body mass (kg), breath direction
- **Random effects**: Dolphin ID

**Rats**
- **Sample size**: 36 observations from 6 rats (3 liver samples per rat)
- **Response variable**: Liver glycogen level
- **Primary predictor**: Treatment
- **Random effects**: Rat and liver samples nested within rat

**Egg hatching**
- **Sample size**: 55 egg batches from 5 replicates and 34 females
- **Response variable**: Hatched vs unhatched eggs
- **Primary predictor**: Genetic cross
- **Random effects**: Replicate and female neseted within replicate

## Files
**Raw data**
- `benzo.csv`: Benzo detoxification dataset
- `dolphins.csv`: Dolphin lung function dataset
- `rats.rds`: Rat glycogen nested dataset
- `egg_hatch_mixed.csv`: Egg hatching dataset

**Analysis scripts**
- benzo_mixed_models.R: Linear mixed-effects modelling workflow
- dolphins_mixed_models.R: Random intercept mixed model
- rats_rested_models.R: Nested random effects example
- eggs_hatching_glmm_models.R: Binomial GLMM and overdispersion analysis

**Outputs**
- `figures/`: Saved model visualisations
- `README.md`: Project documentation

## Analysis Workflow
1. **Exploratory analysis**: 
   - Visualised grouped data structure (batches, individuals, nested samples)
   - Compared pooled relationships vs grouped patterns
   - Identified non-independence requiring mixed-effects modelling
   
2. **Model specification**: 
   - Identified fixed vs random effects
   - Fitted pooled, fixed-effects, and mixed-effects models
   - Specified nested random effects where appropriate
   - Estimate marginal and conditional R² values
   
3. **Model comparison and diagnostics**: 
   - Compared alternative models (pooled vs mixed; crossed vs nested; GLM vs GLMM)
   - Used AIC to evaluate model fit
   - Tested for overdispersion in binomial models
   - Assessed assumptions using residual diagnostics

4. **Final models and interpretation**: 
   - Linear mixed-effects models for continuous outcomes
   - Nested random-effects models for hierarchical sampling
   - Binomial GLMM to resolve overdispersion
   - Interpreted fixed effects, variance components, ICC, and model fit statistics
   
5. **Key findings**: 
   **Benzo**
   - Detoxification capacity increased significantly with toxin concentration (+2.03 units per µM)
   - Experimental batch accounted for substantial variance ((R²m = 0.10 vs R²c = 0.70)
   - Ignoring batch would underestimate variablility and overstate precision
   
   **Dolphins**
   - Tidal volume increased with body mass (+0.017 L per kg)
   - Inhalation volumes were significantly higher than exhalation (+1.11 L)
   - Nearly half of total variance was attributable to differences between dolphins (ICC ≈ 0.47), justifying a random intercept
   
   **Rats**
   - No significant treatment effect on liver glycogen levels
   - Variation occurred both between rats and between liver samples within rats
   - Correcly specifying nested random effects was essential to reflect the sampling design
   
   **Egg hatching**
   - Standard binomial GLM showed clear overdispersion
   - Including replicate and female (nested) random effects significantly improved fit
   - Overdispersion was resolved in the final GLMM, with most variance occurring between females.

## Software
R version: 4.5.2

**Key packages**
- tidyverse v2.0.0
- here v1.0.2
- performance v0.16.0
- lmerTest. v3.2-0
- emmeans v2.0.1
- ggeffects v2.3.2
- MuMIn v1.48.11
- sjPlot v2.9.0
- lme4 v1.1-38


## Author
Teagan Pearse
03/03/2026
