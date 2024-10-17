# Johan S. Sáenz, University of Hohenheim
# Calculate the number of defence systems in the collection of genomes


# Load libraries ---------------------------------------------------------------
pacman::p_load(
  tidyverse,
  data.table,
  pheatmap,
  ggridges,
  ggtext,
  ggdist,
  gghalves,
  patchwork,
  viridis
)

# Concatenate files-----------------------------------------------------
all_paths <- list.files(
  path = "defence_rummen/rawdata/padloc/csv/",
  pattern = "*.csv", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
    header = TRUE,
    sep = ",",
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
names(all_result)[20] <- "genome" # change column name

#  Load other data frames ------------------------------------------------------
quality_genomes <- read_tsv("defence_rummen/rawdata/quality_report_genomes.tsv") %>%
  rename(genome = Name)

taxonomy <- read_tsv("defence_rummen/rawdata/new_gtdbtk.bac120.tsv") %>%
  rename(genome = user_genome)

hq_genomes <- read_tsv("defence_rummen/rawdata/hq_genomes.txt")

# Clean taxonomy data frame
taxonomy <- taxonomy %>%
  select(genome, classification) %>%
  rename_all(tolower) %>%
  mutate(classification = str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
    into = c(
      "kingdom", "phylum", "class", "order", "family", "genus",
      "species"
    ),
    sep = ";"
  ) %>%
  inner_join(hq_genomes)

# Size of genomes
quality_genomes %>%
  filter(Completeness >= 90) %>%
  summarise(
    mean = mean(Genome_Size),
    sd = sd(Genome_Size)
  )


# Create a file for DS families-------------------------------------------------
unique_system <- all_result %>%
  select(system) %>%
  unique()

write.table(unique_system, "rawdata/list_system.txt",
  sep = "\t",
  row.names = FALSE
)

list_system <- read_tsv("defence_rummen/rawdata/list_system.txt")

# Select colors
colors <- c("#1f78b4", "#1b9e77")

# Prevalence plot
kingdom_number <- taxonomy %>%
  select(kingdom) %>%
  count(kingdom, name = "total")

# Barplot of number of DS per genome (distribution)
bar_number_system_genomes <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  select(kingdom, n) %>%
  count(n, kingdom,
    name = "total"
  ) %>%
  group_by(kingdom) %>%
  # filter(kingdom == "Bacteria") %>%
  mutate(perc = 100 * (total / sum(total))) %>%
  ungroup() %>% 
  filter(n < 21) %>%
  ggplot(aes(
    x = n,
    y = perc,
    fill = kingdom
  )) +
  geom_col() +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 26),
    breaks = seq(0, 26, 3)
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 20, 2)
  ) +
  labs(
    y = "Number of genomes (%)",
    x = "Number of defence systems"
  ) +
  scale_fill_manual(values = colors) +
  facet_grid(rows = vars(kingdom)) +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position.inside = c(0.8, 0.95),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    panel.spacing = unit(2, "lines")
  )



# Plot of the correlation between genome size and DS
genome_size_ds <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  mutate(n = replace_na(n, 0)) %>%
  inner_join(taxonomy, by = "genome") %>%
  select(genome, kingdom, n) %>%
  left_join(quality_genomes, by = "genome") %>%
  select(kingdom, Genome_Size, n) %>%
  mutate(Genome_Size = Genome_Size / 1000000) %>%
  ggplot(aes(
    x = Genome_Size,
    y = n,
    fill = kingdom
  )) +
  geom_point(
    position = position_jitter(),
    shape = 21,
    alpha = 0.6,
    color = "black"
  ) +
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(0, 8, 1)) +
  scale_y_continuous(breaks = seq(0, 35, 5)) +
  scale_fill_manual(values = colors) +
  facet_grid(rows = vars(kingdom)) +
  labs(
    y = "Number of defence systems",
    x = "Genome size (Mbp)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  )

# Number DS in genomes with proviruses
number_system_provirus <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  select(genome, kingdom, n)


# Proviruses

# Load the viral prediction file
viral_qc <- read_tsv("defence_rummen/rawdata/viral_quality_summary.tsv")

# Number of proviruses in the hg genomes
provirus <- viral_qc %>%
  filter(
    provirus == "Yes",
    viral_genes >= 1,
    proviral_length > 10000
  ) %>%
  select(contig_id) %>%
  mutate(
    contig_id = str_replace_all(contig_id, "_sc.*", ""),
    contig_id = str_replace_all(contig_id, "_\\d+.*", ""),
    contig_id = str_replace_all(contig_id, "_deg.*", "")
  ) %>%
  group_by(contig_id) %>%
  summarise(proviruses = n(), .groups = "drop") %>%
  rename(genome = contig_id)

# Plot of number of proviruses and DS pero genome
ds_provirus <- left_join(taxonomy, provirus, by = "genome") %>%
  select(genome, proviruses) %>%
  left_join(number_system_provirus, provirus, by = "genome") %>%
  mutate(
    proviruses = if_else(is.na(proviruses), 0, proviruses),
    proviruses = factor(proviruses, levels = c(0, 1, 2, 3, 4, 5))
  ) %>%
  ggplot(aes(
    x = proviruses,
    y = n,
    fill = kingdom
  )) +
  geom_violin(
    aes(fill = kingdom, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  # geom_point(
  #   shape = 21,
  #   position = position_jitter(
  #     width = 0.2,
  #     seed = 0
  #   ),
  #   size = 3,
  #   alpha = 0.8,
  #   color = "black"
  # ) +
  # stat_summary(
  #   fun = median,
  #   show.legend = FALSE,
  #   geom = "crossbar",
  #   color = "black",
  #   linewidth = 0.5
  # ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 35),
    breaks = seq(0, 35, 5),
    oob = scales::squish
  ) +
  labs(
    y = "Defence system per genome",
    x = "Number of proviruses"
  ) +
  facet_grid(rows = vars(kingdom)) +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    panel.spacing = unit(2, "lines")
  )


# Make compose plot ------------------------------------------------------------
figure2 <-  ((bar_number_system_genomes / genome_size_ds) | (ds_provirus)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave(figure2, file = "defence_rummen/new_plots/Figure_2.png", width = 10, height = 8, dpi = 400)


# Prevalence of system in all genomes ------------------------------------------

# Select more abundant DS
system_list_filter <- all_result %>%
  filter(!grepl("_other", system)) %>%
  select(genome, system) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  left_join(taxonomy, by = "genome") %>%
  select(genome, system, kingdom) %>%
  group_by(system, kingdom) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(system, kingdom, fill = list(n = 0)) %>%
  mutate(total_genomes = case_when(
    kingdom == "Bacteria" ~ 2976,
    kingdom == "Archaea" ~ 62
  )) %>%
  mutate(rel_abund = 100 * (n / total_genomes)) %>%
  group_by(system) %>%
  summarise(max = max(rel_abund), .groups = "drop") %>%
  filter(max > 3.5) %>%
  pull(system)

# Relative abundance of DS
system_prevalence_genomes <- all_result %>%
  filter(!grepl("_other", system)) %>%
  select(genome, system) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  left_join(taxonomy, by = "genome") %>%
  select(genome, system, kingdom) %>%
  group_by(system, kingdom) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(system, kingdom, fill = list(n = 0)) %>%
  mutate(total_genomes = case_when(
    kingdom == "Bacteria" ~ 2976,
    kingdom == "Archaea" ~ 62
  )) %>%
  mutate(rel_abund = 100 * (n / total_genomes)) %>%
  filter(system %in% system_list_filter) %>%
  mutate(system = str_replace_all(system, "_", " "))

# Order by abundance
order <- system_prevalence_genomes %>%
  arrange(rel_abund) %>%
  pull(system) %>%
  unique()

# Barplot of prevalence most abundant DS type
barplot_prevalence <- system_prevalence_genomes %>%
  mutate(system = factor(system, levels = order)) %>%
  ggplot(aes(
    y = system,
    x = rel_abund,
    fill = kingdom
  )) +
  geom_bar(
    stat = "identity",
    position = "dodge"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 65),
    breaks = seq(0, 65, 5)
  ) +
  scale_fill_manual(values = colors) +
  labs(
    y = NULL,
    x = "Number of genomes (%)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.2)
  )


# DS relative abundance
rel_abund_system_kingdom <- all_result %>% # replace status or kingdom
  filter(!grepl("_other", system)) %>%
  select(system.number, genome, system) %>%
  unique() %>%
  count(genome, system) %>%
  mutate(genome = str_replace(genome, ".fasta_padloc.csv", "")) %>%
  left_join(taxonomy, by = "genome") %>%
  left_join(list_system, by = "system") %>%
  select(genome, system_family, kingdom, n) %>%
  group_by(kingdom, system_family) %>%
  summarise(
    total = sum(n),
    .groups = "drop"
  ) %>%
  group_by(kingdom) %>%
  mutate(rel_abun = 100 * total / sum(total)) %>% 
  ungroup()

# Merge less abundant systems
taxon_pool <- rel_abund_system_kingdom %>%
  group_by(system_family) %>%
  summarise(
    pool = mean(rel_abun) < 1.5,
    mean = mean(rel_abun),
    .groups = "drop"
  )

# Set color for barplot
col_list <- c(
  "lightblue", "#8dd3c7", "#ffffb3", "#fb8072", "#bebada", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "black"
)


# Barplot of relative abundance of DS families
rel_abun_ds_family <- inner_join(rel_abund_system_kingdom, taxon_pool, by = "system_family") %>%
  mutate(system_family = if_else(pool, "Other", system_family)) %>%
  group_by(kingdom, system_family) %>%
  summarise(
    rel_abun = sum(rel_abun),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = kingdom,
    y = rel_abun,
    fill = system_family
  )) +
  geom_bar(
    stat = "identity",
    position = "stack"
  ) +
  scale_fill_manual(values = col_list) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    y = "Relative abundance (%)",
    x = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_markdown(size = 18),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -25)
  )

# Make compose plot ------------------------------------------------------------
figure2 <- rel_abun_ds_family + barplot_prevalence +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18))

ggsave(figure2,
  file = "defence_rummen/plots/figure2.png", width = 12, height = 8
)



# Heatmap phylum defence system ------------------------------------------------



# Number of phylum per genome
phylum_count <- taxonomy %>%
  count(phylum, name = "phylum_number")

# Find the best quantile for color of hearmap
quantile(heat_map_df$percent, 0.95)

# Wrangel df for heatmap
heat_map_df <- all_result %>%
  select(genome, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome, system) %>%
  right_join(taxonomy, by = "genome") %>%
  select(phylum, system, n) %>%
  group_by(phylum, system) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  complete(phylum, system, fill = list(total = 0)) %>%
  filter(!is.na(system)) %>%
  inner_join(phylum_count, by = "phylum") %>%
  group_by(phylum) %>%
  mutate(percent = 100 * (total / phylum_number)) %>%
  ungroup() %>% 
  mutate(
    percent2 = case_when(
      percent >= 20 ~ 20,
      T ~ percent
    ),
    phylum = str_replace_all(phylum, "_", " "),
    system = str_replace_all(system, "_", " ")
  ) %>%
  filter(phylum_number >= 5)

# Df for selecting most abundant DS
remove_system <- heat_map_df %>%
  group_by(system) %>%
  summarise(remove = max(percent) >= 10)

# Heatmap plot
heat_map_plot <- inner_join(remove_system, heat_map_df, by = "system") %>%
  filter(remove == TRUE) %>%
  ggplot(aes(
    x = system,
    y = phylum,
    fill = percent2
  )) +
  geom_tile() +
  geom_text(aes(label = round(percent)),
    color = "white",
    size = 4, angle = 90
  ) +
  scale_fill_gradientn(
    colours = c("white", "blue"),
    breaks = c(0, 5, 10, 15, 20),
    labels = c("0%", "5%", "10%", "15%", ">20%")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.title = element_blank()
  )



# Save heatmap plot
ggsave(heat_map_plot,
  file = "plots/figure3.png", width = 16, height = 8
)

# Defence sytems at genus level ------------------------------------------------

# Select genera with more than 10 genomes
more_5_genomes <- taxonomy %>%
  count(genus, name = "number_genomes") %>%
  filter(number_genomes >= 10)

# Make a list of those genera
select_genus <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  group_by(genome) %>%
  summarise(total = mean(n), .groups = "drop") %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(total = replace_na(total, 0)) %>%
  select(family, genus, total) %>%
  group_by(family, genus) %>%
  summarise(
    mean = mean(total),
    sd = sd(total), .groups = "drop"
  ) %>%
  filter(
    genus %in% more_5_genomes$genus,
    !genus == "",
    mean > 6
  ) %>%
  select(genus)

# Plot DS per genome
genus_ds_plot <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  group_by(genome) %>%
  summarise(total = mean(n), .groups = "drop") %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(total = replace_na(total, 0)) %>%
  select(family, genus, total) %>%
  filter(genus %in% select_genus$genus) %>%
  mutate(
    genus = str_replace(genus, "_", " "),
    new_name = str_c(family, ": ", genus),
    new_name = str_replace(new_name, " ", "<br>")
  ) %>%
  ggplot(aes(
    x = total, digits = 1,
    y = reorder(new_name, total)
  )) +
  geom_point(
    shape = 21,
    position = position_jitter(
      height = 0.1,
      seed = 0
    ),
    size = 2,
    alpha = 0.9,
    fill = "grey"
  ) +
  stat_summary(
    fun = mean,
    show.legend = FALSE,
    geom = "crossbar",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 2)
  ) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Number of defence system",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_markdown()
  )

# Number of defence system families by genus
genus_ds_fam_plot <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  inner_join(list_system, by = "system") %>%
  select(genome, system_family) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  group_by(genome) %>%
  summarise(total = mean(n), .groups = "drop") %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(total = replace_na(total, 0)) %>%
  select(family, genus, total) %>%
  filter(genus %in% select_genus$genus) %>%
  mutate(
    genus = str_replace(genus, "_", " "),
    new_name = str_c(family, ": ", genus),
    new_name = str_replace(new_name, " ", "<br>")
  ) %>%
  ggplot(aes(
    x = total, digits = 1,
    y = reorder(new_name, total)
  )) +
  geom_point(
    shape = 21,
    position = position_jitter(
      height = 0.1,
      seed = 0
    ),
    size = 2,
    alpha = 0.9,
    fill = "grey"
  ) +
  stat_summary(
    fun = mean,
    show.legend = FALSE,
    geom = "crossbar",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 14),
    breaks = seq(0, 14, 2)
  ) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Number of defence system families",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_markdown()
  )

# Make compose plot
figure4 <- genus_ds_fam_plot + genus_ds_plot +
  plot_annotation(tag_levels = "A")

# Save plot
ggsave(figure4,
  file = "plots/figure4.png",
  width = 12, height = 6
)


# Correlation ds and ds families -----------------------------------------------

# Number of ds
genus_ds_plot <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  group_by(genome) %>%
  summarise(total = mean(n), .groups = "drop") %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(total = replace_na(total, 0)) %>%
  select(family, genus, total) %>%
  filter(genus %in% select_genus$genus) %>%
  group_by(genus) %>%
  summarise(mean_total = mean(total), .groups = "drop")

# Number of ds families
genus_ds_fam <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  inner_join(list_system, by = "system") %>%
  select(genome, system_family) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  group_by(genome) %>%
  summarise(total = mean(n), .groups = "drop") %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(total = replace_na(total, 0)) %>%
  select(family, genus, total) %>%
  filter(genus %in% select_genus$genus) %>%
  group_by(genus) %>%
  summarise(mean_total = mean(total), .groups = "drop")

# Correlation plot
key_genus <- inner_join(genus_ds_fam_plot, genus_ds_plot, by = "genus") %>%
  rename(system = mean_total.y, family = mean_total.x) %>%
  ggplot(aes(
    x = family,
    y = system
  )) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth"
  ) +
  geom_point() +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 14),
    breaks = seq(0, 14, 1)
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 6),
    breaks = seq(0, 6, 1)
  ) +
  labs(
    y = "Number of defence system",
    x = "Number of defence system families"
  ) +
  theme_classic()

# Save correlation plot
ggsave(key_genus, file = "plots/key_member.png", dpi = 400, width = 6, height = 4)
