
# ─────────────────────────────────────────────
# Title: Fishing for data: AI approaches to advance recreational fisheries monitoring.
# Author: Baker LR, Knott NA, Gorkin R, Aubin S, Brown C, Peters KJ 
# Date: 2025
# Description: Statistical modeling of fish classification confidence
# ─────────────────────────────────────────────

# Clean environment
rm(list = ls())

# ───────────────────────────────────────────────────────────────
# Load required packages
# ───────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(stringr)
library(performance)
library(lme4)
library(see)
library(patchwork)
library(glmm)
library(MASS)
library(glmmTMB)
library(AICcmodavg)
library(broom)
library(sjPlot)
library(gt)
library(installr)
library(reshape2)

# ───────────────────────────────────────────────────────────────
# Load and clean fish classification confidence data
# ───────────────────────────────────────────────────────────────

# -----  set your working directory -----

#load data
Rdata <- read.csv("Rdata.csv")
summary(Rdata)

# Remove "Groundtruth" from fish_angle
clean <- Rdata %>%
  filter(fish_angle != "Groundtruth")

# Remove rows where confidence is zero
withoutzero <- clean %>%
  filter(correct_classification_confidence > 0)

# Summarise average confidence by species
withoutzerospecies <- withoutzero %>%
  group_by(GTfish) %>%
  summarise(
    ave = mean(correct_classification_confidence),
    sd = sd(correct_classification_confidence),
    n = n(),
    se = sd / sqrt(n())
  )

# ───────────────────────────────────────────────────────────────
# Bar chart of classification confidence by species
# ───────────────────────────────────────────────────────────────

species.barchart <- ggplot(withoutzerospecies, aes(x = reorder(GTfish, -ave), y = ave, fill = GTfish)) +
  geom_col(colour = "black", position = 'dodge') +
  geom_errorbar(aes(ymin = ave - se, ymax = ave + se), width = 0.2) +
  scale_fill_manual(values = c(
    Snapper = "#FFA3A3", Flathead = "#996600",
    Bream = "#FFFFCC", Leatherjacket = "#82EEF2"
  )) +
  xlab("Species") +
  ylab("Mean Identification Confidence") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

species.barchart

# ───────────────────────────────────────────────────────────────
# GLM model: First pass using Poisson regression
# ───────────────────────────────────────────────────────────────

glm.confidence <- glm(
  correct_classification_confidence ~ fish_angle + fish_orientation + camera_height + GTfish,
  family = poisson(link = log),
  data = clean
)

summary(glm.confidence)

# NOTE:
# Residual deviance = 21.734 on 304 degrees of freedom
# → Suggests under-dispersion, likely from autocorrelation of repeated species
# → Consider GLMM with fishID as a random effect

# ───────────────────────────────────────────────────────────────
# Log-transform the response variable to improve normality
# ───────────────────────────────────────────────────────────────

logconfidence <- withoutzero %>%
  mutate(correct_classification_confidence = log(correct_classification_confidence + 1))

# ───────────────────────────────────────────────────────────────
# GLMM with beta distribution (logit link)
# ───────────────────────────────────────────────────────────────

# Model 1: Camera height + orientation + angle
tmb <- glmmTMB(
  correct_classification_confidence ~ camera_height + fish_orientation + fish_angle + (1 | fishID),
  family = beta_family(),
  data = logconfidence
)

# Model 2: Camera height + angle + species
tmb2 <- glmmTMB(
  correct_classification_confidence ~ camera_height + fish_angle + GTfish + (1 | fishID),
  family = beta_family(),
  data = logconfidence
)

# Model 3: Interaction between species and camera height
tmb3 <- glmmTMB(
  correct_classification_confidence ~ camera_height + GTfish * camera_height + (1 | fishID),
  family = beta_family(),
  data = logconfidence
)

# ───────────────────────────────────────────────────────────────
# Model selection and comparison
# ───────────────────────────────────────────────────────────────

# Store candidate models
Cand.set <- list(tmb, tmb2, tmb3)
Modnames <- c("tmb", "tmb2", "tmb3")

# Compute AICc model selection table
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)

# Evidence ratios
evidence(aic.table = aicctable.out)

# Model summaries and diagnostics
aicctable.out
summary(tmb)
summary(tmb2)
AIC(tmb2)

check_model(tmb)
check_model(tmb2)

# ───────────────────────────────────────────────────────────────
# CHAPTER 2: Length and Width Error Estimation
# ───────────────────────────────────────────────────────────────

# Load error datasets
errordata <- read.csv("Bellambi_error_direction.csv")
lengthwidth <- read.csv("Bellambi_length_width.csv")

# ───────────────────────────────────────────────────────────────
# Distribution check for error data
# ───────────────────────────────────────────────────────────────

ggplot(lengthwidth, aes(x = width_err_pc)) +
  geom_histogram()

# ───────────────────────────────────────────────────────────────
# GLMs for length and width errors
# ───────────────────────────────────────────────────────────────

# Length error model
glm1 <- glm(
  length_err_pc ~ Species + resolution + Species * resolution,
  family = quasipoisson(link = "log"),
  data = lengthwidth
)

# Width error model
glm2 <- glm(
  width_err_pc ~ Species + resolution + Species * resolution,
  family = quasipoisson(link = "log"),
  data = lengthwidth
)

# Model summaries and diagnostics
summary(glm1)
check_model(glm1)

summary(glm2)
check_model(glm2)

# Model tables
tab_model(glm1, glm2)

# ───────────────────────────────────────────────────────────────
# Visualise directional error by species
# ───────────────────────────────────────────────────────────────

sorted_error <- errordata %>%
  group_by(species, direction) %>%
  summarise(
    mean = mean(error),
    sd = sd(error),
    n = n(),
    se = sd / sqrt(n())
  )

bar <- ggplot(sorted_error, aes(x = direction, y = mean, fill = species)) +
  geom_col(colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  facet_grid(~species) +
  scale_fill_manual(values = c(
    "Snapper" = "#FFA3A3",
    "Dusky flathead" = "#996600",
    "Tiger flathead" = "orange",
    "Blue mackerel" = "skyblue"
  )) +
  labs(x = "Measurement", y = "Mean Error %", fill = "Species") +
  theme_classic() +
  theme(
    strip.text.x = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.8, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

bar

