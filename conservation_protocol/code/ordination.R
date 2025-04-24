### Code by Johan S. Sáenz
### University of Hohenheim


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(vegan)
library(here)
library(patchwork)
library(broom)
library(gt)

# Load data --------------------------------------------------------------------
taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>%
  rename(genome = Genome)

# Fresh samples were remove because they can not be compared with the other
# treatments

clean_metadata <- readRDS("conservation_protocol/output/clean_meta.rds") %>%
  filter(treatment != "F") %>%
  mutate(substrate = if_else(substrate == "None", "No substrate", substrate))

proteins <- read_tsv("conservation_protocol/rawdata/functions_report.tsv")


F <- "#8c510a"
C <- "#d8b365"
FDN <- "#2166ac"
FD <- "#67a9cf"
FRN <- "#1a9850"
FR <- "#91cf60"


base_color <- c("#d8b365", "#67a9cf", "#2166ac", "#91cf60", "#1a9850")

# Wrangle taxonomy-----------
df_taxonomy <- taxonomy %>%
  select(genome, starts_with("Intensity ")) %>%
  pivot_longer(-genome, names_to = "sample", values_to = "intensity") %>%
  mutate(
    sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
    sample = str_remove(sample, "_WDH"),
    sample = str_remove(sample, "_STneu"),
    sample = str_remove(sample, "_wdh"),
    sample = str_remove(sample, "_neuST")
  ) %>%
  mutate(
    rel_abun = 100 * (intensity / sum(intensity)),
    .by = sample,
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_metadata, by = "sample") %>%
  pivot_wider(id_cols = genome, names_from = sample, values_from = rel_abun)

# Create matrix-----------------------------------------------------------------
matrix_taxonomy <- df_taxonomy[c(2:ncol(df_taxonomy))] %>%
  t()

dist <- vegdist(matrix_taxonomy, method = "bray")

adonis_tax <- adonis2(as.dist(dist) ~ clean_metadata$treatment * clean_metadata$incubation * clean_metadata$substrate, permutations = 1000)


tidy(adonis_tax) %>% 
  gt() %>% 
  gtsave("conservation_protocol/output/adonis_tax.docx")

# Calculate nMDS
set.seed(1)
nmds1_taxonomy <- metaMDS(matrix_taxonomy, # perform nmds
  distance = "bray",
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3
)

nmds_best_taxonomy <- metaMDS(matrix_taxonomy,
  distance = "bray", # find best nmds
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3,
  previous.best = nmds1_taxonomy
)

nmds_best_taxonomy$stress

# Extract ndms data points -----------------------------------------------------
data_scores_tax <- as.data.frame(scores(nmds_best_taxonomy, display = c("sites"))) %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(sample = as.numeric(sample))


data_nmds_tax <- data_scores_tax %>%
  left_join(clean_metadata, by = "sample") %>%
  filter(treatment != "FD_noN3Wo")

centroid_tax <- data_nmds_tax %>%
  group_by(treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  )


tax_plot <- data_nmds_tax %>%
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = incubation
  )) +
  geom_point(size = 3) +
  geom_point(
    data = centroid_tax, size = 3,
    shape = 21,
    color = "black",
    aes(fill = treatment)
  ) +
  scale_color_manual(
    values = base_color,
    guide = guide_legend(override.aes = list(color = base_color))
  ) +
  scale_fill_manual(values = base_color) +
  scale_shape_manual(values = c(25, 24, 23)) +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

tax_plot_feed <- data_nmds_tax %>%
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = substrate
  )) +
  geom_point(size = 3) +
  geom_point(
    data = centroid_tax, size = 3,
    shape = 21,
    color = "black",
    aes(fill = treatment)
  ) +
  scale_color_manual(
    values = base_color,
    guide = guide_legend(override.aes = list(color = base_color))
  ) +
  scale_fill_manual(values = base_color) +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

ordination_tax <- ordination_plot <- tax_plot + tax_plot_feed + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(ordination_tax, file ="conservation_protocol/plots/ordination_tax.png",
       width = 10, height = 4, dpi = 300)

# Ordination core proteins -----------------------------------------------------

protein_names <- proteins %>%
  filter(Protein_ID == 1) %>%
  select(Name, starts_with("Intensity ")) %>%
  pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>%
  mutate(
    intensity_presence = if_else(intensity > 0, 1, 0),
    sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
    sample = str_remove(sample, "_WDH"),
    sample = str_remove(sample, "_STneu"),
    sample = str_remove(sample, "_wdh"),
    sample = str_remove(sample, "_neuST"),
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_metadata, by = "sample") %>%
  group_by(Name) %>%
  summarise(sum_presence = sum(intensity_presence)) %>%
  filter(sum_presence > 71) %>%
  pull(Name)


protein_filter <- proteins %>%
  filter(Protein_ID == 1) %>%
  select(Name, starts_with("Intensity ")) %>%
  pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>%
  filter(Name %in% protein_names) %>%
  mutate(
    sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
    sample = str_remove(sample, "_WDH"),
    sample = str_remove(sample, "_STneu"),
    sample = str_remove(sample, "_wdh"),
    sample = str_remove(sample, "_neuST"),
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_metadata, by = "sample") %>%
  mutate(
    rel_abun = 100 * (intensity / sum(intensity)),
    .by = sample
  ) %>%
  pivot_wider(id_cols = Name, names_from = sample, values_from = rel_abun)

matrix_proteins_filtered <- protein_filter[c(2:ncol(protein_filter))] %>%
  t()

# Calculate nMDS
set.seed(2)
nmds1_proteins_filtered <- metaMDS(matrix_proteins_filtered, # perform nmds
  distance = "bray",
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3
)

nmds_best_proteins_filtered <- metaMDS(matrix_proteins_filtered,
  distance = "bray", # find best nmds
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3,
  previous.best = nmds1_proteins_filtered
)

nmds_best_proteins_filtered$stress

# Extract ndms data points -----------------------------------------------------
data_scores_proteins_filtered <- as.data.frame(scores(nmds_best_proteins_filtered, display = c("sites"))) %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(sample = as.numeric(sample))

data_nmds_proteins_filtered <- data_scores_proteins_filtered %>%
  left_join(clean_metadata, by = "sample")


centroid <- data_nmds_proteins_filtered %>%
  group_by(treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  )

prot_70p_incubation <- data_nmds_proteins_filtered %>%
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = incubation
  )) +
  geom_point(size = 3) +
  geom_point(
    data = centroid, size = 3,
    shape = 21,
    color = "black",
    aes(fill = treatment)
  ) +
  scale_color_manual(
    values = base_color,
    guide = guide_legend(override.aes = list(color = base_color))
  ) +
  scale_shape_manual(values = c(25, 24, 23)) +
  scale_fill_manual(values = base_color) +
  theme(
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

prot_70p_substrate <- data_nmds_proteins_filtered %>%
  mutate(substrate = factor(substrate,
                            levels = c("No substrate", "Hay", "Rapeseed meal",
                                       "Wheat grain"))) %>% 
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = substrate
  )) +
  geom_point(size = 3) +
  geom_point(
    data = centroid, size = 3,
    shape = 21,
    color = "black",
    aes(fill = treatment)
  ) +
  scale_color_manual(
    values = base_color,
    guide = guide_legend(override.aes = list(color = base_color))
  ) +
  scale_fill_manual(values = base_color) +
  scale_shape_manual(values = c(18, 17, 15, 3)) + 
  theme(
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# save plot --------------------------------------------------------------------
ordination_protein <-  prot_70p_incubation + prot_70p_substrate + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(ordination_protein,
  filename = "conservation_protocol/plots/ordination_protein.tif",
  width = 10, height = 4, dpi = 600)

# Adonis test ----
dist <- vegdist(matrix_proteins_filtered, method = "bray")

adonis_prot <- adonis2(as.dist(dist) ~ clean_metadata$treatment * clean_metadata$incubation * clean_metadata$substrate, permutations = 1000)


tidy(adonis_prot) %>% 
  gt() %>% 
  gtsave(filename = "conservation_protocol/output/adonis_prot.docx")


#------------------------------------------------------------------------------

df_proteins <- proteins %>%
  filter(Protein_ID == 1) %>%
  select(Name, starts_with("Intensity ")) %>%
  pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>%
  mutate(
    sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
    sample = str_remove(sample, "_WDH"),
    sample = str_remove(sample, "_STneu"),
    sample = str_remove(sample, "_wdh"),
    sample = str_remove(sample, "_neuST")
  ) %>%
  mutate(
    rel_abun = 100 * (intensity / sum(intensity)),
    .by = sample
  ) %>%
  pivot_wider(id_cols = Name, names_from = sample, values_from = rel_abun)

matrix_proteins <- df_proteins[c(2:ncol(df_proteins))] %>%
  t()

# Calculate nMDS
nmds1_proteins <- metaMDS(matrix_proteins, # perform nmds
  distance = "bray",
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3
)

nmds_best_proteins <- metaMDS(matrix_proteins,
  distance = "bray", # find best nmds
  try = 20,
  trymax = 100,
  maxit = 1000,
  k = 3,
  previous.best = nmds1_proteins
)

nmds_best_proteins$stress

# Extract ndms data points -----------------------------------------------------
data_scores_proteins <- as.data.frame(scores(nmds_best_proteins, display = c("sites")))

data_scores_proteins$sample <- as.numeric(row.names(data_scores_proteins))

data_nmds_proteins <- data_scores_proteins %>%
  left_join(clean_meta, by = "sample") %>%
  filter(treatment != "FD_noN3Wo")



proteins_plot <- data_nmds_proteins %>%
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = incubation
  )) +
  geom_point(
    alpha = 0.7,
    size = 3
  ) +
  scale_color_manual(values = base_color) +
  labs(title = "Functions") +
  theme_classic() +
  theme(legend.title = element_blank())

proteins_plot_feed <- data_nmds_proteins %>%
  ggplot(aes(
    x = NMDS1,
    y = NMDS2,
    color = treatment,
    shape = substrate
  )) +
  geom_point(
    alpha = 0.7,
    size = 3
  ) +
  scale_color_manual(values = base_color) +
  labs(title = "Functions") +
  theme_classic() +
  theme(legend.title = element_blank())

ordination <- (tax_plot + proteins_plot) / (tax_plot_feed + proteins_plot_feed)

ggsave(ordination,
  file = "conservation_protocol/plots/ordination.png",
  width = 12, height = 8
)
