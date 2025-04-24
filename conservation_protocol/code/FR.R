# Code by Johan S. Sáenz
# University of Hohenheim

# Load libraries ---------------------------------------------------------------
set.seed(150988)

library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)
library(gt)

base_color <- c("#d8b365", "#67a9cf", "#2166ac", "#91cf60", "#1a9850")

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



res.aov <- aov(nFR ~ treatment*incubation*substrate, data = fr_df)
  sumary <- summary(res.aov)

tidy(res.aov) %>% 
  gt() %>% 
  gtsave(filename = "conservation_protocol/output/aov_fr.docx")



leveneTest(nFR ~ treatment, data = fr_df)

TukeyHSD(res.aov) 



incubation_plot <- fr_df %>%
  ggplot(aes(
    x = incubation,
    y = nFR
  )) +
  geom_boxplot(outlier.colour = "white", width = 0.4) +
  geom_jitter(width = 0.3,
              aes(color = treatment)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, 0.01)
  ) +
  scale_color_manual(values = base_color) +
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
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_x_discrete(labels=c("None" = "None", "Hay" = "Hay", "Rapeseed meal" = "Rapeseed\nmeal", 
                            "Wheat grain" = "Wheat\ngrain"))


compose_plot <- treatment_plot + incubation_plot + substrate_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(compose_plot,
       filename = "conservation_protocol/plots/FR.tif", width = 11, height = 4,
       dpi = 600)
