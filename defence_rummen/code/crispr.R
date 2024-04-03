# Johan S. Sáenz, University of Hohenheim
# Calculate the number of CRISPR-Cas arrays and their relative abundance

# Load libraries ---------------------------------------------------------------
pacman::p_load(
  tidyverse,
  data.table,
  patchwork
)

# Load data and colors ---------------------------------------------------------
hq_genomes <- read_tsv("defence_rummen/rawdata/hq_genomes.txt")

df_taxo <- read_tsv("defence_rummen/rawdata/new_gtdbtk.bac120.tsv") %>%
  select(genome = user_genome, classification) %>%
  separate(classification,
    into = c("domain", "phylum", "class", "order", "family", "genues", "species"),
    sep = ";"
  ) %>%
  select(genome, phylum) %>%
  mutate(phylum = str_remove(phylum, "p__")) %>%
  inner_join(hq_genomes)



colors <- c("#1f78b4", "#1b9e77")

# Concatenate files ------------------------------------------------------------
all_paths <- list.files(
  path = "defence_rummen/rawdata/crisparcas/crispr_cas_loci/",
  pattern = "*.tab", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
    header = TRUE,
    sep = "\t",
    encoding = "UTF-8"
  )

all_filenames <- all_paths %>% # Manipulate file paths
  basename() %>%
  as.list()

all_lists <- mapply(c,
  all_content,
  all_filenames,
  SIMPLIFY = FALSE
)

all_result <- rbindlist(all_lists, fill = T)
names(all_result)[9] <- "genome" # change column name

# Wrangle data -----------------------------------------------------------------
# Count hq genome numbers
total_phy <- df_taxo %>%
  mutate(phylum = str_remove(phylum, "_[A-Z]")) %>%
  count(phylum)

# Relative abundance  CRISPR by phylum
cas_phylum_genomes <- all_result %>%
  select(genome) %>%
  mutate(genome = str_remove(genome, "_CRISPR_Cas\\.tab")) %>%
  unique() %>%
  inner_join(df_taxo, by = "genome") %>%
  count(phylum) %>%
  mutate(phylum = str_remove(phylum, "_[A-Z]")) %>%
  group_by(phylum) %>%
  summarise(n = sum(n)) %>%
  inner_join(total_phy, by = "phylum") %>%
  mutate(
    percentage = 100 * (n.x / n.y),
    phylum = fct_reorder(phylum, percentage)
  ) %>%
  ggplot(aes(
    x = percentage,
    y = phylum
  )) +
  geom_col() +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 50),
    breaks = seq(0, 50, 10)
  ) +
  labs(
    x = "Genomes with CRISPR-Cas loci\n(%)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.text.y = element_text(size = 9)
  )

# Colors for barplot relative abundance
crispr_color <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "black", "lightgreen"
)

# Relative abundance of system subtypes
crispr_abundance <- all_result %>%
  select(genome, Prediction) %>%
  mutate(genome = str_remove(genome, "_CRISPR_Cas\\.tab")) %>%
  left_join(df_taxo, by = "genome") %>%
  mutate(phylum = str_remove(phylum, "_[A-Z]")) %>%
  count(phylum, Prediction) %>%
  mutate(
    rel_abun = 100 * (n / sum(n)),
    .by = phylum
  )

# Merge low abundant subsytems
pool_crispr <- crispr_abundance %>%
  summarise(
    max = max(rel_abun),
    mean = mean(rel_abun),
    .by = Prediction
  ) %>%
  mutate(pool = if_else(max > 7, FALSE, TRUE)) %>%
  select(Prediction, pool)

# Create a ffactor
phylum_factor <- c(
  "UBP6", "Actinomycetota", "Bacteroidota", "Verrucomicrobiota",
  "Eremiobacterota", "Desulfobacterota", "Fibrobacterota", "Bacillota",
  "Methanobacteriota", "Pseudomonadota", "Elusimicrobiota", "Spirochaetota",
  "Cyanobacteriota"
)


# Plot of subtypes relative abundance
plot_crispr_abundance <- crispr_abundance %>%
  inner_join(pool_crispr, by = "Prediction") %>%
  mutate(Prediction = if_else(pool, "Other", Prediction)) %>%
  summarise(
    sum_rel_abund = sum(rel_abun),
    .by = c(Prediction, phylum)
  ) %>%
  mutate(phylum = factor(phylum, levels = (phylum_factor))) %>%
  ggplot(aes(
    y = fct_rev(phylum),
    x = sum_rel_abund,
    fill = Prediction
  )) +
  geom_col() +
  scale_x_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = crispr_color) +
  labs(
    y = NULL,
    x = "Relative abundance (%)"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    text = element_text(size = 8)
  )

# Number of CRISPR per Domain
df_crispr <- data.frame(
  domain = c("Archaea", "Bacteria"),
  crispr = c(9, 659),
  total_genome = c(62, 2976)
)


# Plot for number of CRISPR
all_genomes_crispr <- df_crispr %>%
  mutate(percentage = 100 * (crispr / total_genome)) %>%
  ggplot(aes(
    x = domain,
    y = percentage,
    fill = domain
  )) +
  geom_col() +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_fill_manual(values = colors) +
  labs(
    x = NULL,
    y = "Genomes with CRISPR-Cas loci (%)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

# Compose plot --------------------------------------------------
crispr_fig <- all_genomes_crispr + cas_phylum_genomes + plot_crispr_abundance +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12))

ggsave(crispr_fig, file = "plots/crispr_fig.png", width = 8, height = 4, dpi = 450)


# Phylum factor
phylum_factor <- all_result %>%
  select(genome) %>%
  mutate(genome = str_remove(genome, "_CRISPR_Cas\\.tab")) %>%
  unique() %>%
  inner_join(df_taxo, by = "genome") %>%
  count(phylum) %>%
  mutate(phylum = str_remove(phylum, "_[A-Z]")) %>%
  group_by(phylum) %>%
  summarise(n = sum(n)) %>%
  inner_join(total_phy, by = "phylum") %>%
  mutate(
    percentage = 100 * (n.x / n.y),
    phylum = fct_reorder(phylum, percentage)
  ) %>%
  pull(phylum)
