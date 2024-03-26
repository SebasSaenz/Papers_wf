# Code by Johan S. Sáenz
# University of Hohenheim

# Load libraries ---------------------------------------------------------------
set.seed(150988)

library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)

# Load data --------------------------------------------------------------------

clean_metadata <- readRDS("conservation_protocol/output/clean_meta.rds") %>%
  filter(treatment != "F")

fr <- read_csv("conservation_protocol/rawdata/Record_FR_PCN.csv")

#------------------------------------------------------------------------------
fr_df <- fr %>%
  rename(sample = Sample) %>%
  mutate(
    sample = str_remove(sample, "Nr[0*]"),
    sample = str_remove(sample, "Nr"),
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_metadata, by = "sample")



fr_df %>%
  group_by(treatment, incubation, substrate) %>%
  identify_outliers(nFR)


model <- aov(nFR ~ treatment * incubation * substrate, data = fr_df)
summary(model)

# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


fr_df %>%
  group_by(treatment, incubation, substrate) %>%
  shapiro_test(nFR)



ggqqplot(fr_df, "nFR", ggtheme = theme_bw()) +
  facet_grid(treatment + incubation ~ substrate, labeller = "label_both")




model <- lm(nFR ~ incubation * substrate * treatment, data = fr_df)
summary(model)

treatment.effect <- fr_df %>%
  group_by(incubation, substrate) %>%
  anova_test(nFR ~ treatment, error = model)




fr_df %>%
  group_by(treatment) %>%
  anova_test(nFR ~ incubation * substrate, error = model)

fr_df %>%
  ggplot(aes(
    x = treatment,
    y = nFR
  )) +
  geom_boxplot(outlier.colour = "white", width = 0.4) +
  geom_jitter(aes(
    color = incubation,
    shape = substrate
  ), width = 0.3) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, 0.01)
  ) +
  scale_color_manual(values = c('#e66101', '#fdb863', '#b2abd2')) +
  labs(x = NULL) +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = "conservation_protocol/plots/FR.png", width = 5, height = 3)
