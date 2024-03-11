set.seed(150988)

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(ggtext)


# Load data --------------------------------------------------------------------
clean_metadata <- read_rds("conservation_protocol/output/clean_meta.rds") %>% 
  filter(!treatment == "F")

peptides <- read_rds("conservation_protocol/output/peptides.rds")

base_color <- c('#d8b365','#67a9cf', '#2166ac','#91cf60','#1a9850')

#
df <- peptides %>%
  pivot_longer(-Proteins, names_to = "sample", values_to = "Intensity") %>%
  group_by(sample) %>%
  summarise(sum_intensity = sum(Intensity)) %>%
  mutate(
    log2 = log10(sum_intensity),
    sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
    sample = str_remove(sample, "_WDH"),
    sample = str_remove(sample, "_STneu"),
    sample = str_remove(sample, "_wdh"),
    sample = str_remove(sample, "_neuST"),
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_metadata, by = "sample")

# Statistical test -------------------------------------------------------------
model <- aov(log2 ~ treatment * incubation * substrate, data = df)
x <- summary(model)

anova_results <- tidy(model)
write_tsv(anova_results, file = "conservation_protocol/output/anova_peptides.txt")

q# Plot ------------------------------------------------------------
plot_incubation <- df %>% 
  mutate(incubation = str_replace(incubation, "-", "<br>"),
         incubation = factor(incubation,
                             levels = c("Pre<br>incubation", "8", "24"))) %>% 
    ggplot(aes(
      x = incubation,
      y = log2
    )) +
    geom_boxplot(width = 0.5) +
    geom_jitter(aes(colour = treatment),
                width = 0.2,
                alpha = 0.8
    ) +
    scale_y_continuous(
      limits = c(10.90, 11.30),
      breaks = seq(10.90, 11.30, 0.05)
    ) +
    scale_color_manual(values = base_color) +
    labs(
      x = NULL,
      y = "Peptides intensity (Log2)"
    ) +
    theme(
      panel.background = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text.x = element_markdown()
    )

plot_substrate <- df %>%
  mutate(substrate = str_replace(substrate, "-", "<br>"),
         substrate = str_replace(substrate, " ", "<br>"),
         substrate = factor(substrate,
                             levels = c("None", "Hay", "Rapeseed<br>meal", "Wheat<br>grain"))) %>% 
  ggplot(aes(
    x = substrate,
    y = log2
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(aes(colour = treatment),
              width = 0.2,
              alpha = 0.8
  ) +
  scale_y_continuous(
    limits = c(10.90, 11.30),
    breaks = seq(10.90, 11.30, 0.05)
  ) +
  scale_color_manual(values = base_color) +
  labs(
    x = NULL,
    y = "Peptides intensity (Log2)"
  ) +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing  = unit(0, "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_markdown()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

compose_plot <- plot_incubation + plot_substrate + plot_layout(axes = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))

# Save plot --------------------------------------------------------------------
ggsave(compose_plot, file = "conservation_protocol/plots/peptides_plot.png",
       dpi = 450,
       width = 6.5,
       height = 4)
